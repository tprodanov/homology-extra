#!/usr/bin/env python3

import sys
import pysam
import argparse
import warnings
import itertools
import numpy as np
from scipy.stats import pearsonr, spearmanr, ConstantInputWarning, NearConstantInputWarning
from collections import defaultdict, Counter
from intervaltree import IntervalTree
import gzip

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
        shared_alleles[c1[allele_ix1]] = var1.alleles[allele_ix1]
    """
    assert var1.ref == var2.ref
    shared_alleles = [var1.ref]
    alleles_set = set()
    for allele in var1.alts:
        if allele not in alleles_set:
            shared_alleles.append(allele)
            alleles_set.add(allele)
    for allele in var2.alts:
        if allele not in alleles_set:
            shared_alleles.append(allele)
            alleles_set.add(allele)

    c1 = [shared_alleles.index(allele) for allele in var1.alleles]
    c2 = [shared_alleles.index(allele) for allele in var2.alleles]
    assert c1[0] == c2[0] == 0
    return shared_alleles, c1, c2


def fmt_regions(regions, var):
    curr_regions = set(reg.data for reg in regions[var.chrom].overlap(var.start, var.start + len(var.ref)))
    return ';'.join(curr_regions) if curr_regions else '*'


def gt_str(gt):
    return '/'.join(map(str, gt))


def process_var(wgs_var, wes_var, samples, genes, exons, out1, out2, depth, quality):
    chrom = wgs_var.chrom
    start = wgs_var.start
    assert chrom == wes_var.chrom and start == wgs_var.start
    if wgs_var.ref != wes_var.ref or len(wgs_var.alleles) != len(wes_var.alleles):
        return False

    alleles, c1, c2 = reconcile_vars(wgs_var, wes_var)
    n_alleles = len(alleles)
    # Genotype cns for each applicable sample.
    cns = []
    # Allele counts for each allele and each applicable sample.
    wgs = [[] for _ in range(n_alleles)]
    wes = [[] for _ in range(n_alleles)]

    events = Counter()
    for sample in samples:
        wgs_call = wgs_var.samples[sample]
        wes_call = wes_var.samples[sample]
        if (wgs_call.get('GQ') or 0) < quality or (wes_call.get('GQ') or 0) < quality:
            continue
        wgs_gt = wgs_call.get('GT')
        wes_gt = wes_call.get('GT')
        if wgs_gt is None or wgs_gt[0] is None or wes_gt is None or wes_gt[0] is None:
            continue
        cn = len(wgs_gt)
        if len(wes_gt) != cn:
            sys.stderr.write('Genotype lengths do not match at '
                f'{chrom}:{start+1} for {sample}: {gt_str(wgs_gt)} and {gt_str(wes_gt)}\n')
            continue
        if wgs_call.get('DP') < depth * cn or wes_call.get('DP') < depth * cn:
            continue

        cns.append(cn)
        wgs_gt = tuple(sorted(c1[allele_ix] for allele_ix in wgs_gt))
        wes_gt = tuple(sorted(c2[allele_ix] for allele_ix in wes_gt))
        events[(wgs_gt, wes_gt)] += 1
        for i in range(n_alleles):
            wgs[i].append(wgs_gt.count(i))
            wes[i].append(wes_gt.count(i))

    cns = np.array(cns)
    wgs = np.array(wgs)
    wes = np.array(wes)

    count = len(cns)
    ref_cn = 2 + 2 * len(wgs_var.info['pos2'])
    mean_cn = np.mean(cns) if count else np.nan
    prefix = f'{chrom}:{start+1}\t{alleles[0]}\t{",".join(alleles[1:])}\t'
    prefix += f'{fmt_regions(genes, wgs_var)}\t{fmt_regions(exons, wgs_var)}\t'
    prefix += f'{wgs_var.info["overlPSV"]}\t{ref_cn}\t{mean_cn:.5f}'
    if not count:
        out1.write(prefix)
        out1.write('\n')
        return True

    af1 = wgs_var.info['AF']
    af2 = wes_var.info['AF']

    af1 = [np.nan] + [0] * (n_alleles - 1)
    for i, val in enumerate(wgs_var.info['AF'], 1):
        af1[c1[i]] = val
    af2 = [np.nan] + [0] * (n_alleles - 1)
    for i, val in enumerate(wes_var.info['AF'], 1):
        af2[c2[i]] = val

    sum_cns = np.sum(cns)
    for i in range(n_alleles):
        out1.write(f'{prefix}\t{i}\t{af1[i]:.5f}\t{af2[i]:.5f}\t')
        a = wgs[i]
        b = wes[i]
        differences = np.bincount(np.abs(a - b))
        out1.write('{}\t{}\t'.format(count, ','.join(map(str, differences))))
        out1.write('{}\t{}\t'.format(np.sum(cns),
            np.sum(cns) - np.sum(np.maximum(a, b)) + np.sum(np.minimum(a, b))))
        if count >= 5:
            out1.write('{:.6f}\t{:.6f}\n'.format(pearsonr(a, b).statistic, spearmanr(a, b).statistic))
        else:
            out1.write('NA\tNA\n')

    for (wgs_gt, wes_gt), count in events.most_common():
        out2.write('{}:{}\t{}\t{}\t{}\n'
            .format(chrom, start + 1, '/'.join(map(str, wgs_gt)), '/'.join(map(str, wes_gt)), count))
    return True


def process_vcfs(wgs_vcf, wes_vcf, genes, exons, out1, out2, depth, quality):
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
                processed += process_var(wgs_var, wes_var, samples, genes, exons, out1, out2, depth, quality)
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
    parser.add_argument('-o', '--output', required=True, metavar='STR',
        help='Output prefix.')
    parser.add_argument('-d', '--depth', metavar='FLOAT', type=float, default=5,
        help='Only compare variant calls with read depth over FLOAT * cn [%(default)s].')
    parser.add_argument('-q', '--quality', metavar='FLOAT', type=float, default=10,
        help='Quality threshold for comparison [%(default)s].')
    args = parser.parse_args()

    sys.stderr.write('Loading genes and exons\n')
    with common.open_possible_gzip(args.exons) as f:
        genes, exons = load_bed(f)

    warnings.simplefilter('ignore', (ConstantInputWarning, NearConstantInputWarning))

    o_fname1 = f'{args.output}.csv.gz'
    o_fname2 = f'{args.output}.events.csv.gz'
    with gzip.open(o_fname1, 'wt') as out1, gzip.open(o_fname2, 'wt') as out2, \
            pysam.VariantFile(args.wgs) as wgs_vcf, pysam.VariantFile(args.wes) as wes_vcf:
        out1.write('# {}\n'.format(' '.join(sys.argv)))
        out1.write(f'# depth threshold = {args.depth} * cn\n')
        out1.write(f'# quality threshold = {args.quality}\n')
        out1.write('pos\tref\talts\tgenes\texons\tpsv\tref_cn\tmean_cn\tallele_ix\t')
        out1.write('wgs_AF\twes_AF\tcount\tdifferences\tsum_length\tsum_match'
            '\tpearson\tspearman\n')
        out2.write('pos\twgs_gt\twes_gt\tcount\n')

        sys.stderr.write('Processing VCF files\n')
        process_vcfs(wgs_vcf, wes_vcf, genes, exons, out1, out2, args.depth, args.quality)


if __name__ == '__main__':
    main()
