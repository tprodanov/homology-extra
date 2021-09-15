#!/usr/bin/env python3

import pysam
import gzip
import sys
import argparse


def load_sample_depth(depth, sample):
    res = [None] * 101
    for line in depth:
        if line.startswith('#'):
            continue
        curr_sample, gc_count, depth1, depth2 = line.split('\t')
        if curr_sample != sample:
            continue
        res[int(gc_count)] = float(depth1)
    return res


def calc_gc_count(seq):
    return seq.count('G') + seq.count('C')


def analyze_window(line, fasta, bg_depth, out):
    line = line.split('\t')
    if line[0] != 'window':
        return
    positions = [line[1]]
    if line[2] != '*':
        positions.extend(reg[:-2] for reg in line[2].split(' '))
    seqs = [fasta.fetch(region=pos) for pos in positions]
    info = line[3].split(';')
    key, window_ix = info[0].split('=')
    assert key == 'window_ix'

    copy_num = len(seqs)
    out.write('{}\t{}\t{}\t{}\t'.format(window_ix, copy_num * 2, len(seqs[0]), sum(map(len, seqs))))
    gc_counts = list(map(calc_gc_count, seqs))
    sum_gc_counts = sum(gc_counts)
    out.write('{}\t{}\t{:.2f}\t'.format(gc_counts[0], sum_gc_counts, sum_gc_counts / gc_counts[0] / copy_num))

    bg_depth1 = bg_depth[gc_counts[0]]
    sum_bg_depth = sum(bg_depth[gc_count] for gc_count in gc_counts)
    out.write('{:.3f}\t{:.3f}\t{:.3f}\n'.format(bg_depth1, sum_bg_depth, sum_bg_depth / bg_depth1 / copy_num))


def main():
    parser = argparse.ArgumentParser(
        description='Count GC content of the windows and their copies.',
        formatter_class=argparse.RawTextHelpFormatter, add_help=False,
        usage='%(prog)s -i <joint_data> -f <fasta> -o <csv>')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', metavar='<file>', required=True,
        help='Input CSV.GZ joint data file.')
    io_args.add_argument('-f', '--fasta-file', metavar='<file>', required=True,
        help='Input FASTA file with reference sequence.')
    io_args.add_argument('-d', '--depth', metavar='<file>', type=argparse.FileType(), required=True,
        help='Input CSV file with background read depth.')
    io_args.add_argument('-s', '--sample', metavar='<str>', required=True,
        help='Sample name.')
    io_args.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<file>', required=False, default='-',
        help='Outpu CSV file.')

    oth_args = parser.add_argument_group('Other arguments')
    oth_args.add_argument('-h', '--help', action='help', help='Show this message and exit.')
    args = parser.parse_args()

    depth = load_sample_depth(args.depth, args.sample)

    out = args.output
    out.write('# {}\n'.format(' '.join(sys.argv)))
    out.write('window_ix\tploidy\tlen1\tsum_len\tgc_count1\tsum_gc_count\tgc_count_cons\tbg_depth1\tsum_bg_depth\tbg_depth_cons\n')
    fasta = pysam.FastaFile(args.fasta_file)
    with gzip.open(args.input, 'rt') as inp:
        for line in inp:
            analyze_window(line, fasta, depth, out)


if __name__ == '__main__':
    main()
