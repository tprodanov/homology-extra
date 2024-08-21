#!/usr/bin/env python3

import sys
import pysam
import argparse
import warnings
import itertools
import numpy as np
from scipy.stats import pearsonr, spearmanr, ConstantInputWarning, NearConstantInputWarning
from collections import defaultdict
from intervaltree import IntervalTree

import parascopy.inner.common as common


def load_bed(f):
    genes = defaultdict(IntervalTree)
    exons = defaultdict(IntervalTree)

    for line in f:
        line = line.strip().split('\t')
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        if end <= start:
            continue
        gene = line[3]
        gene_ty = line[4]
        if gene_ty != 'protein_coding':
            gene += '*'
        ty = line[5]
        if ty == 'gene':
            genes[chrom].addi(start, end, gene)
        elif ty == 'exon':
            exons[chrom].addi(start, end, gene)
    return genes, exons


def reconcile_vars(var1, var2):
    """
    Reconciles allele sets of two variants.
    Returns combined set of alleles, as well index schemes:
        alleles[i] = var1.alleles[c1[i]]
    """
    assert var1.ref == var2.ref
    alleles = [var1.ref]
    c1 = [0]
    c2 = [0]
    for allele in itertools.chain(var1.alts, var2.alts):
        if allele in alleles:
            continue
        alleles.append(allele)
        try:
            c1.append(var1.alleles.index(allele))
        except ValueError:
            c1.append(None)
        try:
            c2.append(var2.alleles.index(allele))
        except ValueError:
            c2.append(None)
    return alleles, c1, c2


def fmt_regions(regions, var):
    curr_regions = set(reg.data for reg in regions[var.chrom].overlap(var.start, var.start + len(var.ref)))
    return ';'.join(curr_regions) if curr_regions else '*'


def gt_str(gt):
    return '/'.join(map(str, gt))


def format_af(var, i):
    if i is None:
        return '0'
    elif i == 0:
        return 'NA'
    else:
        return f'{var.info["AF"][i - 1]:.5f}'


def process_var(wgs_var, wes_var, samples, genes, exons, out, depth, quality):
    chrom = wgs_var.chrom
    start = wgs_var.start
    assert chrom == wes_var.chrom and start == wgs_var.start
    if wgs_var.ref != wes_var.ref or len(wgs_var.alleles) != len(wes_var.alleles):
        return False

    alleles, c1, c2 = reconcile_vars(wgs_var, wes_var)
    n_alleles = len(alleles)
    # Genotype ploidies for each applicable sample.
    ploidies = []
    # Allele counts for each allele and each applicable sample.
    wgs = [[] for _ in range(n_alleles)]
    wes = [[] for _ in range(n_alleles)]

    for sample in samples:
        wgs_call = wgs_var.samples[sample]
        wes_call = wes_var.samples[sample]
        if (wgs_call.get('GQ') or 0) < quality or (wes_call.get('GQ') or 0) < quality:
            continue
        wgs_gt = wgs_call.get('GT')
        wes_gt = wes_call.get('GT')
        if wgs_gt is None or wgs_gt[0] is None or wes_gt is None or wes_gt[0] is None:
            continue
        ploidy = len(wgs_gt)
        if len(wes_gt) != ploidy:
            sys.stderr.write('Genotype lengths do not match at '
                f'{chrom}:{start+1} for {sample}: {gt_str(wgs_gt)} and {gt_str(wes_gt)}\n')
            continue
        if wgs_call.get('DP') < depth * ploidy or wes_call.get('DP') < depth * ploidy:
            continue
        ploidies.append(ploidy)
        for i in range(n_alleles):
            wgs[i].append(wgs_gt.count(c1[i]))
            wes[i].append(wes_gt.count(c2[i]))

    ploidies = np.array(ploidies)
    wgs = np.array(wgs)
    wes = np.array(wes)

    count = len(ploidies)
    mean_ploidy = np.mean(ploidies) if count else np.nan
    prefix = f'{chrom}\t{start+1}\t{alleles[0]}\t{",".join(alleles[1:])}\t'
    prefix += f'{fmt_regions(genes, wgs_var)}\t{fmt_regions(exons, wgs_var)}\t'
    prefix += f'{wgs_var.info["overlPSV"]}\t{mean_ploidy:.5f}'
    if not count:
        out.write(prefix)
        out.write('\n')
        return True

    sum_ploidies = np.sum(ploidies)
    for i in range(n_alleles):
        out.write(f'{prefix}\t{i}\t{format_af(wgs_var, c1[i])}\t{format_af(wes_var, c2[i])}\t')
        a = wgs[i]
        b = wes[i]
        differences = np.bincount(np.abs(a - b))
        out.write('{}\t{}\t'.format(count, ','.join(map(str, differences))))
        out.write('{}\t{}\t{}\t{}\t'.format(
            np.sum((a == 0) & (b == 0)), np.sum((a == 0) & (b > 0)),
            np.sum((a > 0) & (b == 0)), np.sum((a > 0) & (b > 0)),
        ))
        out.write('{}\t{}\t'.format(np.sum(ploidies),
            np.sum(ploidies) - np.sum(np.maximum(a, b)) + np.sum(np.minimum(a, b))))
        if count >= 5:
            out.write('{:.6f}\t{:.6f}\n'.format(pearsonr(a, b).statistic, spearmanr(a, b).statistic))
        else:
            out.write('NA\tNA\n')
    return True


def process_vcfs(wgs_vcf, wes_vcf, genes, exons, out, depth, quality):
    samples = sorted(set(wgs_vcf.header.samples) & set(wes_vcf.header.samples))
    total1 = 0
    total2 = 0
    processed = 0
    try:
        wgs_var = next(wgs_vcf)
        total1 += 1
        wes_var = next(wes_vcf)
        total2 += 1

        while True:
            if (wgs_var.rid, wgs_var.start) < (wes_var.rid, wes_var.start):
                wgs_var = next(wgs_vcf)
                total1 += 1
            elif (wes_var.rid, wes_var.start) < (wgs_var.rid, wgs_var.start):
                wes_var = next(wes_vcf)
                total2 += 1
            else:
                processed += process_var(wgs_var, wes_var, samples, genes, exons, out, depth, quality)
                wgs_var = next(wgs_vcf)
                wes_var = next(wes_vcf)
                total1 += 1
                total2 += 1
    except StopIteration:
        pass

    sys.stderr.write(f'Processed {processed:,} variants. Total variants: {total1:,} (WGS), {total2:,} (WES)\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--wgs', required=True, metavar='FILE',
        help='WGS VCF file.')
    parser.add_argument('--wes', required=True, metavar='FILE',
        help='WES VCF file.')
    parser.add_argument('-e', '--exons', required=True, metavar='FILE',
        help='BED file with exon information. 4-6 columns: gene, gene type, record type (gene/exon).')
    parser.add_argument('-o', '--output', required=True, metavar='FILE',
        help='Output CSV file.')
    parser.add_argument('-d', '--depth', metavar='FLOAT', type=float, default=5,
        help='Only compare variant calls with read depth over FLOAT * ploidy [%(default)s].')
    parser.add_argument('-q', '--quality', metavar='FLOAT', type=float, default=10,
        help='Quality threshold for comparison [%(default)s].')
    args = parser.parse_args()

    sys.stderr.write('Loading genes and exons\n')
    with common.open_possible_gzip(args.exons) as f:
        genes, exons = load_bed(f)

    warnings.simplefilter('ignore', (ConstantInputWarning, NearConstantInputWarning))

    with common.open_possible_gzip(args.output, 'w') as out, \
            pysam.VariantFile(args.wgs) as wgs_vcf, pysam.VariantFile(args.wes) as wes_vcf:
        out.write('# {}\n'.format(' '.join(sys.argv)))
        out.write(f'# depth threshold = {args.depth} * ploidy\n')
        out.write(f'# quality threshold = {args.quality}\n')
        out.write('chrom\tstart\tref\talts\tgenes\texons\tpsv\tmean_ploidy\tallele_ix\t')
        out.write('wgs_AF\twes_AF\tcount\tdifferences\tev00\tev01\tev10\tev11\tsum_length\tsum_match'
            '\tpearson\tspearman\n')
        sys.stderr.write('Processing VCF files\n')
        process_vcfs(wgs_vcf, wes_vcf, genes, exons, out, args.depth, args.quality)


if __name__ == '__main__':
    main()
