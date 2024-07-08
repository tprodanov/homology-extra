#!/usr/bin/env python3

import sys
import pysam
import argparse
import warnings
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


def allele_correspondence(alleles1, alleles2):
    """
    Retuns correspondence from the first to the second set of alleles:
    alleles1[i] == alleles2[corresp[i]]
    """
    try:
        return [alleles2.index(a) for a in alleles1]
    except ValueError:
        return None


def fmt_regions(regions, var):
    curr_regions = set(reg.data for reg in regions[var.chrom].overlap(var.start, var.start + len(var.ref)))
    return ';'.join(curr_regions) if curr_regions else '*'


def gt_str(gt):
    return '/'.join(map(str, gt))


def af_str(var):
    af = var.info.get('AF')
    if af is None:
        return '*'
    return ','.join(map('{:.6f}'.format, af))


def process_var(wgs_var, wes_var, samples, genes, exons, out):
    chrom = wgs_var.chrom
    start = wgs_var.start
    assert chrom == wes_var.chrom and start == wgs_var.start
    if wgs_var.ref != wes_var.ref or len(wgs_var.alleles) != len(wes_var.alleles):
        return False

    corresp = allele_correspondence(wgs_var.alleles, wes_var.alleles)
    if corresp is None:
        return False
    assert corresp[0] == 0

    count = 0
    matches = 0
    mism1 = 0

    sum_match = 0
    sum_size = 0
    n_zeros1 = []
    n_zeros2 = []

    for sample in samples:
        wgs_call = wgs_var.samples[sample]
        wes_call = wes_var.samples[sample]
        wgs_gt = wgs_call.get('GT')
        wes_gt = wes_call.get('GT')
        if wgs_gt is None or wgs_gt[0] is None or wes_gt is None or wes_gt[0] is None:
            continue

        wgs_gt = sorted(wgs_gt)
        wes_gt = sorted(corresp[j] for j in wes_gt)
        # except IndexError:
        #     sys.stderr.write(f'{chrom}\t{start+1}\t{wgs_var.alleles}\t{wes_var.alleles}\t{corresp}\t'
        #         f'{gt_str(wgs_gt)}\t{gt_str(wes_gt)}\n')
        #     raise
        gt_size = len(wgs_gt)
        if len(wes_gt) != gt_size:
            sys.stderr.write('Genotype lengths do not match at '
                f'{chrom}:{start+1} for {sample}: {gt_str(wgs_gt)} and {gt_str(wes_gt)}\n')
            continue
        count += 1
        match_size = sum(i == j for i, j in zip(wgs_gt, wes_gt))
        matches += match_size == gt_size
        mism1 += match_size == gt_size - 1
        sum_size += gt_size
        sum_match += match_size
        n_zeros1.append(sum(i == 0 for i in wgs_gt))
        n_zeros2.append(sum(j == 0 for j in wes_gt))

    biallelic = len(wgs_var.alts) == 1
    out.write(f'{chrom}\t{start+1}\t{wgs_var.ref}\t{",".join(wgs_var.alts)}\t{"FT"[biallelic]}\t')
    out.write(f'{fmt_regions(genes, wgs_var)}\t{fmt_regions(exons, wgs_var)}\t')
    out.write('{}\t{}\t{}\t'.format(wgs_var.info['overlPSV'], af_str(wgs_var), af_str(wes_var)))
    out.write(f'{count}\t{matches}\t{mism1}\t{sum_size}\t{sum_match}\t')
    if count >= 5:
        out.write(f'{pearsonr(n_zeros1, n_zeros2).statistic:.5f}\t{spearmanr(n_zeros1, n_zeros2).statistic:.5f}\n')
    else:
        out.write('nan\tnan\n')
    return True


def process_vcfs(wgs_vcf, wes_vcf, genes, exons, out):
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
                processed += process_var(wgs_var, wes_var, samples, genes, exons, out)
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
    args = parser.parse_args()

    sys.stderr.write('Loading genes and exons\n')
    with common.open_possible_gzip(args.exons) as f:
        genes, exons = load_bed(f)

    warnings.simplefilter('ignore', (ConstantInputWarning, NearConstantInputWarning))

    with common.open_possible_gzip(args.output, 'w') as out, \
            pysam.VariantFile(args.wgs) as wgs_vcf, pysam.VariantFile(args.wes) as wes_vcf:
        out.write('# {}\n'.format(' '.join(sys.argv)))
        out.write('chrom\tstart\tref\talts\tbiallelic\tgenes\texons\tpsv\twgs_AF\twes_AF\t')
        out.write('count\tfull_match\tmism1\tsum_size\tsum_match\tpearson\tspearman\n')
        sys.stderr.write('Processing VCF files\n')
        process_vcfs(wgs_vcf, wes_vcf, genes, exons, out)


if __name__ == '__main__':
    main()
