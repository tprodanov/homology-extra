#!/usr/bin/env python3

import os
import sys
import itertools
import collections
import operator
import argparse
import pysam
import gzip
from intervaltree import IntervalTree

import parascopy.inner.common as common


def load_regions(filename):
    trees = collections.defaultdict(IntervalTree)
    loci = set()
    with common.open_possible_gzip(filename) as inp:
        for line in inp:
            line = line.strip().split('\t')
            locus = line[3] if len(line) > 3 else '+'
            trees[line[0]].addi(int(line[1]), int(line[2]), locus)
            loci.add(locus)
    return trees, loci


def count(vcf_filename, trees):
    """
    For each locus, returns pair (total number, number of SNPs, number of indels) for any and high qualities.
    """
    counts = collections.defaultdict(lambda: [0] * 6)
    with pysam.VariantFile(vcf_filename) as vcf:
        for record in vcf:
            ref_len = len(record.ref)
            overlap = trees[record.chrom].overlap(record.start, record.start + ref_len)
            overlap = set(map(operator.attrgetter('data'), overlap))
            is_snp = all(len(alt) == ref_len for alt in record.alts)
            high_qual = record.samples[0].get('GQ', 10000) >= 10

            for locus in overlap or ('*',):
                counts[locus][0] += 1
                counts[locus][2 - int(is_snp)] += 1
                if high_qual:
                    counts[locus][3] += 1
                    counts[locus][5 - int(is_snp)] += 1
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--eval', metavar='DIR', required=True,
        help='RTG evaluation directory.')
    parser.add_argument('-R', '--regions', metavar='FILE', required=True,
        help='BED file, where fourth column declares region type.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    parser.add_argument('-a', '--all', action='store_true',
        help='Output counts for all entries, even empty.')
    args = parser.parse_args()

    trees, loci = load_regions(args.regions)
    counts_tpb = count(os.path.join(args.eval, 'tp-baseline.vcf.gz'), trees)
    counts_tpc = count(os.path.join(args.eval, 'tp.vcf.gz'), trees)
    counts_fp = count(os.path.join(args.eval, 'fp.vcf.gz'), trees)
    counts_fn = count(os.path.join(args.eval, 'fn.vcf.gz'), trees)

    types = 'any\tall any\tsnps any\tindels high\tall high\tsnps high\tindels'.split(' ')
    with common.open_possible_gzip(args.output, 'w') as out:
        out.write('# {}\n'.format(' '.join(sys.argv)))
        out.write('region\tqual\tvar_type\ttp_base\ttp_call\tfp\tfn\n')
        for locus in itertools.chain(sorted(loci), ('*',)):
            for i, ty in enumerate(types):
                tpb = counts_tpb[locus][i]
                tpc = counts_tpc[locus][i]
                fp = counts_fp[locus][i]
                fn = counts_fn[locus][i]
                if i == 0 and not args.all and tpb + tpc + fp + fn == 0:
                    break
                out.write(f'{locus}\t{ty}\t{tpb}\t{tpc}\t{fp}\t{fn}\n')


if __name__ == '__main__':
    main()
