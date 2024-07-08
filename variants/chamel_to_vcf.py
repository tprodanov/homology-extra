#!/usr/bin/env python3

import sys
import argparse
import pysam
import collections
from intervaltree import IntervalTree

import parascopy.inner.common as common


def create_header(ref):
    header = pysam.VariantHeader()
    for chrom, length in zip(ref.references, ref.lengths):
        header.add_line(f'##contig=<ID={chrom},length={length}>')
    header.add_line('##FORMAT=<ID=GT,Type=String,Number=1,Description="Sample genotype.">')
    header.add_sample('SAMPLE')
    return header


def load_region_cn(filename):
    trees = collections.defaultdict(IntervalTree)
    with common.open_possible_gzip(filename) as inp:
        for line in inp:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            if line[3] != 'PASS':
                continue
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            cn = int(line[4])
            trees[chrom].addi(start, end, cn)
    return trees


    if filename is None:
        return None

    with common.open_possible_gzip(filename) as inp:
        intervals = []
        for line in inp:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            intervals.append(Interval(genome.chrom_id(line[0]), int(line[1]), int(line[2])))
    intervals.sort()
    intervals = Interval.combine_overlapping(intervals)
    return itree.MultiNonOverlTree(intervals)


def combine_overlaps(var1, var2, header):
    assert var1.chrom == var2.chrom
    ref1 = var1.ref
    ref2, alt2 = var2.alleles
    start1 = var1.start
    start2 = var2.start
    assert start1 <= start2
    len1 = len(ref1)
    len2 = len(ref2)
    end1 = start1 + len1
    end2 = start2 + len2

    prefix2 = ref1[: start2 - start1]
    suffix2 = ref1[min(len1, end2 - start1):]
    suffix1 = ref2[min(len2, end1 - start2):]

    new_ref = ref1 + suffix1
    new_ref2 = prefix2 + ref2 + suffix2
    if new_ref != new_ref2:
        print(str(var1).strip())
        print(str(var2).strip())
        print('...{},   {}...{}'.format(suffix1, prefix2, suffix2))
        print(f'{ref1 + suffix1}, {prefix2 + ref2 + suffix2}')
        exit(1)

    var = header.new_record()
    var.chrom = var1.chrom
    alleles = [new_ref]
    for alt1 in var1.alts:
        alleles.append(alt1 + suffix1)
    alleles.append(prefix2 + alt2 + suffix2)

    var.start = start1
    var.alleles = alleles
    var.filter.add('PASS')
    var.qual = 100
    var.samples[0]['GT'] = var1.samples[0]['GT'] + (len(alleles) - 1, )
    return var


def process(inp, header, ref, var_col):
    chrom_order = { chrom: i for i, chrom in enumerate(ref.references) }
    skipped = 0
    incorr_ref = 0

    variants = []
    for line in inp:
        line = line.strip().split('\t')
        chrom = line[0]
        if chrom not in chrom_order:
            skipped += 1
            continue
        start = int(line[1]) - 1
        ref_allele0, alt_allele = line[var_col].rsplit('~', 1)[1].split('/')

        if ref_allele0 == '-' or alt_allele == '-':
            start -= 1
            first_nt = ref.fetch(chrom, start, start + 1)
            ref_allele0 = first_nt + ref_allele0.replace('-', '')
            alt_allele = first_nt + alt_allele.replace('-', '')

        end = start + len(ref_allele0)
        ref_allele = ref.fetch(chrom, start, end)
        if ref_allele0 != ref_allele:
            incorr_ref += 1
        if ref_allele == alt_allele:
            continue

        variant = header.new_record()
        variant.chrom = chrom
        variant.start = start
        variant.alleles = (ref_allele, alt_allele)
        variant.filter.add('PASS')
        variant.qual = 100
        variant.samples[0]['GT'] = (1,)
        variants.append(variant)

    sys.stderr.write(f'Skipped {skipped} variants\n')
    sys.stderr.write(f'{incorr_ref} variants with incorrect reference\n')

    variants.sort(key=lambda var: (chrom_order[var.chrom], var.start))
    merged_vars = []
    last_chrom = None
    last_end = 0
    for var in variants:
        if var.chrom == last_chrom and var.start < last_end:
            merged = combine_overlaps(merged_vars[-1], var, header)
            merged_vars[-1] = merged
            last_end = merged.start + len(merged.ref)
        else:
            assert last_chrom is None or chrom_order[last_chrom] <= chrom_order[var.chrom]
            last_chrom = var.chrom
            last_end = var.start + len(var.ref)
            merged_vars.append(var)
    return merged_vars


def set_cn(vars, trees):
    for var in vars:
        start = var.start
        end = var.start + len(var.ref)
        cns = trees[var.chrom].overlap(start, end)
        n = len(cns)
        if n == 0:
            cn = 4
        elif n == 1:
            cn = next(iter(cns)).data
        else:
            cn = max(obj.data for obj in cns)

        gt = var.samples[0]['GT']
        gt = sorted(gt + (0,) * max(cn - len(gt), 0))
        var.samples[0]['GT'] = gt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='Input Chameleolyzer file.')
    parser.add_argument('-r', '--reference', metavar='FILE', required=True,
        help='Indexed FASTA reference file.')
    parser.add_argument('-c', '--target-cn', metavar='FILE', required=True,
        help='Copy numbers of the target regions (parascopy examine).')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output VCF file.')
    parser.add_argument('--raw', action='store_true',
        help='Input file contains raw variants.')
    args = parser.parse_args()

    ref = pysam.FastaFile(args.reference)
    header = create_header(ref)
    trees = load_region_cn(args.target_cn)

    with open(args.input) as inp, pysam.VariantFile(args.output, 'wz', header=header) as vcf:
        if not args.raw:
            first_line = next(inp)
            assert first_line.startswith('Chromosome\t')
        var_col = 3 if args.raw else 2

        variants = process(inp, header, ref, var_col)
        set_cn(variants, trees)
        for var in variants:
            vcf.write(var)
    pysam.tabix_index(args.output, preset='vcf', force=True)


if __name__ == '__main__':
    main()
