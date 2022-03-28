#!/usr/bin/env python3

import os
import pysam
import argparse
import re
import subprocess


def parse_rg(line):
    rg_id = re.search(r'\bID:([a-zA-Z0-9_-]+)\b', line).group(1)
    sample = re.search(r'\bSM:([a-zA-Z0-9_-]+)\b', line).group(1)
    return rg_id, sample


def create_header(alns, out_chrom):
    sq_entry = None
    rgs_header = ''
    rgs = {}
    for aln_i, aln in enumerate(alns):
        for line in str(aln.header).strip().split('\n'):
            if line.startswith('@SQ') and sq_entry is None:
                if f'SN:{out_chrom}\t' in line:
                    sq_entry = line.strip() + '\n'
            elif line.startswith('@RG'):
                rg_id, sample = parse_rg(line)
                new_rg_id = f'{rg_id}-{aln_i}'
                rgs_header += f'@RG\tID:{new_rg_id}\tSM:{sample}\n'
                rgs[rg_id] = new_rg_id
    new_header = sq_entry + rgs_header
    return new_header, rgs


def load_bed(bed):
    for line in bed:
        line = line.strip().split('\t')
        yield line[0], int(line[1]), int(line[2])


_rev_comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C','a':'T', 't':'A', 'c':'G', 'g':'C', 'N':'N', 'n':'N' }
def rev_comp(seq): # reverse complement of string
    return ''.join(_rev_comp.get(nt, 'X') for nt in seq[::-1])


def cond_rev_comp(seq, *, strand):
    if strand:
        return seq
    return rev_comp(seq)


def cond_reverse(qual, *, strand):
    return qual if strand else qual[::-1]


PADDING = 5500
_MAX_INSERT_SIZE = 50000


def create_record(old_record, header, rgs, new_chrom, new_start, new_cigar, new_strand, status):
    record = pysam.AlignedSegment(header)
    record.query_name = old_record.query_name
    record.query_sequence = cond_rev_comp(old_record.query_sequence, strand=new_strand)
    record.query_qualities = cond_reverse(old_record.query_qualities, strand=new_strand)
    record.mapping_quality = 60
    record.reference_name = new_chrom
    record.reference_start = new_start
    record.cigarstring = new_cigar

    if not new_cigar:
        record.is_unmapped = True

    if old_record.is_paired:
        record.is_paired = True
        if old_record.is_read1:
            record.is_read1 = True
        else:
            record.is_read2 = True

        if old_record.is_proper_pair and old_record.reference_id == old_record.next_reference_id \
                and abs(old_record.reference_start - old_record.next_reference_start) <= _MAX_INSERT_SIZE:
            record.is_proper_pair = True

    oa_tag = '{},{},{},{},{},{};'.format(old_record.reference_name, old_record.reference_start + 1,
        '-' if old_record.is_reverse else '+', old_record.cigarstring, old_record.mapping_quality,
        old_record.get_tag('NM') if old_record.has_tag('NM') else '')
    record.set_tag('OA', oa_tag)

    old_rg = old_record.get_tag('RG')
    record.set_tag('RG', rgs[old_rg])
    record.set_tag('st', status)
    return record


def process_read(record, header, rgs, tg_chrom, tg_start, tg_end, out):
    if record.reference_name == tg_chrom and tg_start - PADDING <= record.reference_start <= tg_end + PADDING:
        out.write(create_record(record, header, rgs,
            record.reference_name, record.reference_start, record.cigarstring, True, 0))

    if not record.has_tag('XA'):
        return
    for xa in filter(bool, record.get_tag('XA').split(';')):
        chrom, pos, cigar, _ = xa.split(',')
        if chrom != tg_chrom:
            continue

        pos = int(pos)
        strand = pos >= 0
        start = abs(pos)
        if tg_start - PADDING <= pos <= tg_end + PADDING:
            new_strand = not strand if record.is_reverse else strand
            out.write(create_record(record, header, rgs, chrom, pos - 1, cigar, new_strand, 1))


def process_aln(aln, header, rgs, dupl_regions, target_region, out):
    tg_chrom, tg_start, tg_end = target_region
    for dp_chrom, dp_start, dp_end in dupl_regions:
        for record in aln.fetch(dp_chrom, max(0, dp_start - PADDING), dp_end + PADDING):
            process_read(record, header, rgs, tg_chrom, tg_start, tg_end, out)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', nargs='+')
    parser.add_argument('-I', '--input-list', type=argparse.FileType())

    parser.add_argument('-b', '--bed', type=argparse.FileType(), required=True)
    parser.add_argument('-r', '--region', required=True)
    parser.add_argument('-f', '--fasta', required=True)
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()

    chrom, pos = args.region.split(':')
    start, end = map(int, pos.split('-'))
    target_region = (chrom, start - 1, end)
    dupl_regions = list(load_bed(args.bed))

    fasta = pysam.FastaFile(args.fasta)

    if args.input is not None:
        alns = [pysam.AlignmentFile(filename, reference_filename=args.fasta, require_index=True)
            for filename in args.input]
    else:
        alns = [pysam.AlignmentFile(filename.split()[0].strip(), reference_filename=args.fasta, require_index=True)
            for filename in args.input_list]

    header, rgs = create_header(alns, chrom)
    tmp_filename = args.output + '.tmp'
    out = pysam.AlignmentFile(tmp_filename, 'wb', text=header)
    header = out.header
    for aln in alns:
        process_aln(aln, header, rgs, dupl_regions, target_region, out)

    out.close()
    subprocess.run(['samtools', 'sort', '-o', args.output, tmp_filename], check=True)
    subprocess.run(['samtools', 'index', args.output], check=True)
    os.remove(tmp_filename)


if __name__ == '__main__':
    main()
