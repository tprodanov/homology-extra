#!/usr/bin/env python3

import argparse
import gzip
import csv
from parascopy.inner import common
from parascopy.inner.genome import Genome, Interval


def process(inp, genome, out, args):
    for line in inp:
        if line.startswith('##'):
            continue
        assert line.startswith('#')
        fields = line.lstrip('#').strip().split('\t')
        break

    reader = csv.DictReader(inp, delimiter='\t', fieldnames=fields)
    res_regions = []
    for row in reader:
        agcn_qual = float(row['agCN_qual'])
        if agcn_qual < args.qual or (not args.all_filters and row['agCN_filter'] != 'PASS'):
            continue
        if not row['agCN'].isdigit():
            continue

        region1 = Interval(genome.chrom_id(row['chrom']), int(row['start']), int(row['end']))
        regions = [region1]
        regions2_str = row['homologous_regions']
        if regions2_str != '*':
            for region_str in regions2_str.split(','):
                assert region_str.endswith(':+') or region_str.endswith(':-')
                regions.append(Interval.parse(region_str[:-2], genome))
        ref_cn = 2 * len(regions)
        if ref_cn < args.cn_bounds[0] or ref_cn > args.cn_bounds[1] or ref_cn != int(row['agCN']):
            continue
        if args.pooled:
            res_regions.append(region1)
            continue

        if not args.all_filters and row['psCN_filter'] != 'PASS':
            continue
        pscns = row['psCN'].split(',')
        pscn_quals = list(map(float, row['psCN_qual'].split(',')))
        assert len(pscns) == len(pscn_quals) == len(regions)

        for pscn, pscn_qual, region in zip(pscns, pscn_quals, regions):
            if pscn.isdigit() and int(pscn) == 2 and pscn_qual >= args.qual:
                res_regions.append(region)

    res_regions.sort()
    res_regions = Interval.combine_overlapping(res_regions)
    return res_regions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, metavar='<file>',
        help='Input res.samples.bed.gz file.')
    parser.add_argument('-f', '--fasta-ref', required=True, metavar='<file>',
        help='Fasta reference file.')
    parser.add_argument('-o', '--output', required=True, metavar='<file>',
        help='Output BED file.')
    parser.add_argument('-q', '--qual', default=20, type=float, metavar='<float>',
        help='Quality threshold [default: %(default)s].')
    parser.add_argument('-c', '--cn-bounds', default=(4, 10), nargs=2, type=int, metavar='<int>',
        help='Minimal and maximal reference copy numbers [default: 4 10].')
    parser.add_argument('-p', '--pooled', action='store_true',
        help='Create pooled calling regions.')
    parser.add_argument('-a', '--all-filters', action='store_true',
        help='Ignore filter values.')
    args = parser.parse_args()

    with gzip.open(args.input, 'rt') as inp, open(args.output, 'w') as out, Genome(args.fasta_ref) as genome:
        regions = process(inp, genome, out, args)
        out.write('# {}\n'.format(common.command_to_str()))
        for region in regions:
            out.write(region.to_bed(genome) + '\n')


if __name__ == '__main__':
    main()
