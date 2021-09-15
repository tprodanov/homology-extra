#!/usr/bin/env python3

import argparse
import pysam
import sys
import re
import numpy as np
import csv


class FastaWrapper:
    def __init__(self, filename):
        self.fasta = pysam.FastaFile(filename)
        self.chrom = None
        self.start = None
        self.end = None
        self.seq = None

    def fetch(self, chrom, start, end):
        if chrom != self.chrom or start < self.start or end > self.end:
            self.seq = self.fasta.fetch(chrom, start, start + 100000).upper()
            self.chrom = chrom
            self.start = start
            self.end = self.start + len(self.seq)
        return self.seq[start - self.start : end - self.start]


def get_read_groups(bam_file):
    """
    Returns list of pairs (group_id, sample).
    """
    read_groups = []
    for line in str(bam_file.header).splitlines():
        if line.startswith('@RG'):
            has_rg = True

            id_m = re.search(r'ID:([ -~]+)', line)
            sample_m = re.search(r'SM:([ -~]+)', line)
            if id_m is None or sample_m is None:
                common.log('ERROR: Cannot find ID or SM field in the header line: "%s"' % line)
                exit(1)
            read_groups.append((id_m.group(1), sample_m.group(1)))
    return read_groups


class Windows:
    def __init__(self, inp):
        reader = csv.DictReader(inp, delimiter='\t')
        starts = []
        ends = []

        for row in reader:
            start = int(row['start'])
            end = int(row['end'])
            window_ix = int(row['window_ix'])
            assert window_ix == len(starts)
            starts.append(start)
            ends.append(end)
        self.starts = np.array(starts)
        self.ends = np.array(ends)

    def find_window(self, center):
        ix = np.searchsorted(self.starts, center, 'right') - 1
        if ix == -1:
            return -0.5
        assert self.starts[ix] <= center
        if self.ends[ix] <= center:
            return ix + 0.5
        return ix



def write_read(read_num, read, read_groups, fasta, windows, out_short, out_long):
    read_seq = read.query_sequence
    ref_chrom = read.reference_name
    ref_start = read.reference_start
    ref_end = read.reference_end
    ref_seq = fasta.fetch(ref_chrom, ref_start, ref_end)

    mism_len = 0
    ins_len = 0
    del_len = 0
    clip_left = 0
    clip_right = 0

    read_pos = 0
    ref_pos = ref_start
    for (op, length) in read.cigartuples:
        if op == 0 or op == 7 or op == 8:
            # M, X or =
            for i in range(length):
                read_nt = read_seq[read_pos]
                ref_nt = ref_seq[ref_pos - ref_start]
                if read_nt != ref_nt and read_nt != 'N' and ref_nt != 'N':
                    out_long.write('{}\t{}\tm\t1\n'.format(read_num, ref_pos + 1))
                    mism_len += 1
                ref_pos += 1
                read_pos += 1

        elif op == 1:
            # Insertion
            out_long.write('{}\t{}\ti\t{}\n'.format(read_num, ref_pos, length))
            ins_len += length
            read_pos += length
        elif op == 2:
            # Deletion
            out_long.write('{}\t{}\td\t{}\n'.format(read_num, ref_pos, length))
            del_len += length
            ref_pos += length
        elif op == 4:
            # Soft clipping
            out_long.write('{}\t{}\t{}\t{}\n'.format(read_num, ref_pos + 1, 's1' if read_pos == 0 else 's2', length))
            if read_pos == 0:
                clip_left += length
            else:
                clip_right += length
            read_pos += length
        else:
            raise RuntimeError('Unexpected operation {} in read {}'.format(op, read.query_name))

    orig_aln = read.get_tag('OA')
    orig_chrom, orig_pos = orig_aln.split(':')
    dist_to_orig = abs(int(orig_pos) - 1 - ref_start) if orig_chrom == ref_chrom else 'inf'
    sample = read_groups[read.get_tag('RG')]
    center = (ref_start + ref_end) // 2
    window_ix = windows.find_window(center)

    out_short.write('{}\t{}\t{}\t{:.1f}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(read_num, read.query_name, sample, window_ix,
        ref_chrom, ref_start, ref_end, orig_aln, dist_to_orig, 'F' if read.is_read2 else 'T'))
    out_short.write('{}\t{}\t{}\t{}\t{}\n'.format(mism_len, ins_len, del_len, clip_left, clip_right))


def main():
    parser = argparse.ArgumentParser(
        description='Extract all errors from pooled reads.',
        formatter_class=argparse.RawTextHelpFormatter, add_help=False,
        usage='%(prog)s -i <bam> -f <fasta> -w <csv> -o <csv>')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', metavar='<file>', required=True,
        help='Input pooled reads BAM file.')
    io_args.add_argument('-f', '--fasta-file', metavar='<file>', required=True,
        help='Input FASTA file with reference sequence.')
    io_args.add_argument('-w', '--windows', metavar='<file>', required=True, type=argparse.FileType(),
        help='Input CSV file with windows.')
    io_args.add_argument('-o', '--output', nargs=2, type=argparse.FileType('w'), metavar='<file>', required=True,
        help='Two output CSV files: first has one line for each read and second has one line for each error.')

    oth_args = parser.add_argument_group('Other arguments')
    oth_args.add_argument('-h', '--help', action='help', help='Show this message and exit.')
    args = parser.parse_args()

    out_short, out_long = args.output
    out_short.write('# {}\n'.format(' '.join(sys.argv)))
    out_short.write('read_num\tread_name\tsample\twindow_ix\tchr\tstart\tend\torig_aln\tdist_to_orig\tis_first\t'
        'mism_len\tins_len\tdel_len\tleft_clip\tright_clip\n')
    out_long.write('# {}\n'.format(' '.join(sys.argv)))
    out_long.write('read_num\tpos\ttype\tlen\n')

    fasta = FastaWrapper(args.fasta_file)
    windows = Windows(args.windows)
    with pysam.AlignmentFile(args.input) as inp:
        read_groups = dict(get_read_groups(inp))
        for read_num, read in enumerate(inp):
            write_read(read_num, read, read_groups, fasta, windows, out_short, out_long)


if __name__ == '__main__':
    main()
