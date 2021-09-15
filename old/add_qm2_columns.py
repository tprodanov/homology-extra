#!/usr/bin/env python3

import argparse
import pysam
from collections import defaultdict
import glob
import numpy as np
import os
import operator
import sys


def load_summary(inp, out):
    by_samples = defaultdict(list)
    for line in inp:
        if line.startswith('#'):
            if line.startswith('##'):
                out.write(line)
            else:
                out.write('{}\tqm2_ploidy\tqm2_paralog\tqm2_entries\tpooled_err\tparalog_err\tparalog_int_err\n'
                    .format(line.strip()))
            continue
        line = line.strip().split('\t')
        line[1] = int(line[1])
        line[2] = int(line[2])
        sample = line[4]
        by_samples[sample].append(line)
    return by_samples


def find_bed_files(qm2_dir):
    res = {}
    for filename in glob.iglob(os.path.join(qm2_dir, '**/*.bed.gz'), recursive=True):
        sample = os.path.basename(filename).split('.', 1)[0]
        res[sample] = filename
    return res


def single_sample_add_columns(summary_entries, qm2_bed, out):
    for line in summary_entries:
        subregions = [tuple(line[:3])]
        if line[5] != '*':
            for region in line[5].split(','):
                chrom, pos, _ = region.split(':')
                start, end = map(int, pos.split('-'))
                subregions.append((chrom, start - 1, end))

        qm2_pooled = 0
        qm2_paralog = []
        qm2_entries_str = []
        for chrom, start, end in subregions:
            qm2_entries = list(qm2_bed.fetch(chrom, start, end)) if qm2_bed is not None else []
            qm2_entries_str.append(','.join(map(operator.itemgetter(3), qm2_entries)))
            sum_ploidy = 0.0
            sum_weight = 0.0
            for qm2_chrom, qm2_start, qm2_end, qm2_ploidy in qm2_entries:
                qm2_start = int(qm2_start)
                qm2_end = int(qm2_end)
                qm2_ploidy = float(qm2_ploidy)
                assert qm2_chrom == chrom
                weight = (min(qm2_end, end) - max(qm2_start, start)) / (end - start)
                sum_weight += weight
                sum_ploidy += qm2_ploidy * weight
            total_qm2_ploidy = sum_ploidy / sum_weight if sum_weight else np.nan
            qm2_pooled += total_qm2_ploidy
            qm2_paralog.append(total_qm2_ploidy)

        out.write('\t'.join(map(str, line)))
        out.write('\t{:.3f}\t{}\t{}\t'.format(qm2_pooled, ','.join(map('{:.3f}'.format, qm2_paralog)),
            '; '.join(qm2_entries_str)))

        pooled_ploidy = int(line[7]) if line[7] != '*' else np.nan
        out.write('{:.3f}\t'.format(pooled_ploidy - qm2_pooled))
        paralog_ploidy = np.array([int(s) if s != '?' else np.nan for s in line[10].split(',')])
        paralog_err = np.sqrt(np.sum(np.power(
            np.where(np.isnan(paralog_ploidy), 0.0, paralog_ploidy - qm2_paralog), 2)))
        paralog_int_err = np.sqrt(np.sum(np.power(
            np.where(np.isnan(paralog_ploidy), 0.0, paralog_ploidy - np.round(qm2_paralog)), 2)))
        out.write('{:.3f}\t{:.0f}\n'.format(paralog_err, paralog_int_err))


def add_columns(by_sample_summary, bed_files, out):
    for i, sample in enumerate(sorted(by_sample_summary), 1):
        sys.stderr.write('{:4}: {}\n'.format(i, sample))
        if sample in bed_files:
            bed_file = pysam.TabixFile(bed_files[sample], parser=pysam.asTuple())
        else:
            bed_file = None
        single_sample_add_columns(by_sample_summary[sample], bed_file, out)
        if bed_file is not None:
            bed_file.close()


def main():
    parser = argparse.ArgumentParser(
        description='Add columns with Quick-Mer2 results.',
        formatter_class=argparse.RawTextHelpFormatter, add_help=False,
        usage='%(prog)s -s <file> -q <dir> -o <file>')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-s', '--summary', type=argparse.FileType(), metavar='<file>', required=True,
        help='Homology tools summary.')
    io_args.add_argument('-q', '--quick-mer2', metavar='<dir>', required=True,
        help='Directory with Quick-Mer2 results.')
    io_args.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<file>', required=True,
        help='Output extended file.')

    oth_args = parser.add_argument_group('Other arguments')
    oth_args.add_argument('-h', '--help', action='help', help='Show this message and exit.')
    args = parser.parse_args()

    bed_files = find_bed_files(args.quick_mer2)
    if not bed_files:
        sys.stderr.write('Cannot find any BED files\n')
        exit(1)
    by_sample_summary = load_summary(args.summary, args.output)
    add_columns(by_sample_summary, bed_files, args.output)


if __name__ == '__main__':
    main()
