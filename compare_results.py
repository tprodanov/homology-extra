#!/usr/bin/env python3

import argparse
import sys
import numpy as np
from collections import defaultdict


class SummaryEntry:
    def __init__(self, line):
        self.line = line = line.strip().split('\t')
        self.chrom = line[0]
        self.start = int(line[1])
        self.end = int(line[2])
        self.sample = line[4]
        self.filter, self.copy_num, self.qual = line[6:9]
        self.paralog_filter, self.paralog_copy_num, self.paralog_qual = line[9:12]

    def covers(self, chrom, pos):
        return self.chrom == chrom and self.start <= pos <= self.end


def load_summary(inp, chrom, positions):
    summaries = defaultdict(lambda: [None] * len(positions))
    for line in inp:
        if line.startswith('#'):
            continue
        entry = SummaryEntry(line)
        for i, pos in enumerate(positions):
            if entry.covers(chrom, pos):
                summaries[entry.sample][i] = entry
    return summaries


def get_positions(gene):
    chrom = None
    pos = None
    if gene == 'FCGR3A':
        chrom = 'chr1'
        pos = (161_545_000,)
    elif gene == 'AMY1C':
        chrom = 'chr1'
        pos = (103_755_000,)
    elif gene == 'SMN1':
        chrom = 'chr5'
        pos = (70_940_000, 70_950_000)
    else:
        sys.stderr.write(f'Cannot find gene {gene}\n')
        exit(1)
    return chrom, pos


def load_populations(inp):
    res = defaultdict(lambda: '*')
    if inp is None:
        return res

    for line in inp:
        line = line.strip().split('\t')
        res[line[1]] = line[6]
    return res


def total_copy_num_distance(entry, b_copy_num):
    if b_copy_num is None:
        return np.nan
    try:
        b_copy_num = float(b_copy_num)
    except ValueError:
        return np.nan

    a_copy_num = entry.copy_num
    if a_copy_num == '*':
        return np.nan
    elif a_copy_num.startswith('>'):
        a_copy_num = float(a_copy_num[1:]) + 1
        return max(0, a_copy_num - b_copy_num)
    elif a_copy_num.startswith('<'):
        a_copy_num = float(a_copy_num[1:]) - 1
        return max(0, b_copy_num - a_copy_num)
    else:
        a_copy_num = float(a_copy_num)
        return abs(a_copy_num - b_copy_num)


def paralog_copy_num_distance(entry, b_paralog):
    if b_paralog is None:
        return np.nan
    a_paralog = entry.paralog_copy_num.split(',')
    assert len(a_paralog) == len(b_paralog)

    dist = 0
    present_both = 0
    for a_val, b_val in zip(a_paralog, b_paralog):
        if a_val == '?' or b_val is None:
            continue
        a_val = float(a_val)

        try:
            b_val = float(b_val)
        except ValueError:
            continue

        dist += abs(a_val - b_val)
        present_both += 1
    if present_both:
        return dist / present_both
    return np.nan


class ResEntry:
    def __init__(self, entry, method, b_copy_num, b_paralog=None):
        self.entry = entry
        self.method = method
        self.b_copy_num = b_copy_num
        self.b_paralog = b_paralog

    def write_to(self, out, populations):
        entry = self.entry
        b_copy_num_str = str(self.b_copy_num) if self.b_copy_num is not None else '*'
        b_paralog_str = ','.join(str(value) if value is not None else '?' for value in self.b_paralog) \
            if self.b_paralog is not None else '*'
        copy_num_dist = total_copy_num_distance(entry, self.b_copy_num)
        paralog_dist = paralog_copy_num_distance(entry, self.b_paralog)

        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(entry.sample, populations[entry.sample],
            entry.filter, entry.copy_num, entry.qual, entry.paralog_filter, entry.paralog_copy_num, entry.paralog_qual))
        out.write('{}\t{}\t{}\t{:.4g}\t{:.4g}\n'.format(
            self.method, b_copy_num_str, b_paralog_str, copy_num_dist, paralog_dist))


def compare_line_fcgr3a(line, entries):
    line = line.strip().split('\t')
    sample = line[0]
    if sample not in entries:
        return
    entry = entries[sample][0]

    # TaqMan
    taq_man_copy_num = line[2]
    taq_man_paralog = (line[3].count('A'), line[3].count('B'))
    yield ResEntry(entry, 'TaqMan', taq_man_copy_num, taq_man_paralog)

    # PRT-REDVR
    prt_copy_num = line[5]
    prt_paralog = (line[6].count('A'), line[6].count('B'))
    yield ResEntry(entry, 'PRT_REDVR', prt_copy_num, prt_paralog)

    # STR
    str_paralog = (None, line[8].count('B'))
    yield ResEntry(entry, 'STR', None, str_paralog)

    # SYBR Green
    sybr_copy_num = line[9]
    yield ResEntry(entry, 'SYBR_Green', sybr_copy_num)


def compare_line_amy1c_2(line, entries):
    line = line.strip().split('\t')
    sample = line[0]
    if sample not in entries:
        return
    entry = entries[sample][0]
    # TODO: What is this method?
    yield ResEntry(entry, '*', line[11])


def compare_line_amy1c_qpcr(line, entries):
    line = line.strip('\n').split('\t')
    sample = line[0]
    if sample not in entries:
        return
    entry = entries[sample][0]

    prt, qpcr, g1k, qpcr_14 = line[1:]
    yield ResEntry(entry, 'PRT', prt)
    yield ResEntry(entry, 'qPCR', qpcr)
    yield ResEntry(entry, 'g1k', g1k)
    yield ResEntry(entry, 'qPCR_14', qpcr_14)


def combine_smn_entries(entry_16, entry_78):
    try:
        paralog_78 = list(map(int, entry_78.paralog_copy_num.split(',')))
        total_16 = int(entry_16.copy_num)
    except ValueError:
        return
    if entry_16.paralog_copy_num != '?,?':
        return

    entry_16.paralog_filter += ';Inferred'
    entry_16.paralog_copy_num = '{},{}'.format(paralog_78[0], total_16 - paralog_78[0])


def compare_line_smn1_caller(line, entries):
    line = line.strip().split('\t')
    sample = line[0]
    if sample not in entries:
        return

    # Two subregions: 16: exons 1-6 and 78: exons 7-8. delta: without exons 7-8.
    entry_16, entry_78 = entries[sample]
    combine_smn_entries(entry_16, entry_78)
    smn1_cn, smn2_full_cn, smn2_delta_cn, total_16, total_78 = line[3:8]
    try:
        smn1_cn = int(smn1_cn)
        smn2_full_cn = int(smn2_full_cn)
        smn2_delta_cn = int(smn2_delta_cn)
        paralog_16 = (smn1_cn, smn2_full_cn + smn2_delta_cn)
        paralog_78 = (smn1_cn, smn2_full_cn)
    except ValueError:
        paralog_16 = None
        paralog_78 = None
    yield ResEntry(entry_16, 'caller_16', total_16, paralog_16)
    yield ResEntry(entry_78, 'caller_78', total_78, paralog_78)


def compare_line_smn1_mlpa(line, entries):
    line = line.strip().split('\t')
    sample = line[2]
    if sample not in entries:
        return

    entry_16, entry_78 = entries[sample]
    combine_smn_entries(entry_16, entry_78)
    total_16, total_78 = line[3:5]
    yield ResEntry(entry_16, 'mlpa_16', total_16)
    yield ResEntry(entry_78, 'mlpa_78', total_78)


def compare_all_smn1_quick_mer2(file, entries):
    header = next(file).strip().split('\t')
    smn2 = next(file).strip().split('\t')
    smn1 = next(file).strip().split('\t')

    chrom, (pos1, pos2) = get_positions('SMN1')
    assert chrom == smn1[0] and int(smn1[1]) <= pos1 < int(smn1[2])

    for i in range(7, len(header)):
        sample = header[i]
        if sample not in entries:
            continue
        entry_16, entry_78 = entries[sample]
        combine_smn_entries(entry_16, entry_78)

        cn1 = float(smn1[i])
        cn2 = float(smn2[i])
        yield ResEntry(entry_16, 'qm2', cn1 + cn2, (cn1, cn2))


def select_function(gene, method):
    if gene == 'FCGR3A':
        assert method is None
        return compare_line_fcgr3a
    elif gene == 'AMY1C' and method == '2':
        return compare_line_amy1c_2
    elif gene == 'AMY1C' and method == 'qPCR':
        return compare_line_amy1c_qpcr
    elif gene == 'SMN1' and method == 'caller':
        return compare_line_smn1_caller
    elif gene == 'SMN1' and method == 'MLPA':
        return compare_line_smn1_mlpa
    elif gene == 'SMN1' and method == 'qm2':
        return compare_all_smn1_quick_mer2
    else:
        sys.stderr.write(f'Cannot find gene {gene} and method {method}\n')
        exit(1)


def compare(b_in, entries, populations, gene, method, out):
    out.write('# {}\n'.format(' '.join(sys.argv)))
    out.write('sample\tpopulation\tcopy_num_filter\tcopy_num\tcopy_num_qual\t'
        'paralog_filter\tparalog_copy_num\tparalog_qual\t'
        'method\tb_copy_num\tb_paralog\ttotal_dist\tparalog_dist\n')

    if gene == 'SMN1' and method == 'qm2':
        b_in = (b_in,)

    compare = select_function(gene, method)
    for line in b_in:
        for res in compare(line, entries):
            res.write_to(out, populations)


def main():
    parser = argparse.ArgumentParser(
        description='Compare results with another method.',
        formatter_class=argparse.RawTextHelpFormatter, add_help=False,
        usage='%(prog)s -a <file> -b <file> [-p <file>] -g <gene> [-m <method>] -o <file>')
    req_args = parser.add_argument_group('Required arguments')
    req_args.add_argument('-a', type=argparse.FileType(), metavar='<file>', required=True,
        help='Homology tools summary.')
    req_args.add_argument('-b', type=argparse.FileType(), metavar='<file>', required=True,
        help='Results of another method.')
    req_args.add_argument('-p', '--populations', type=argparse.FileType(), metavar='<file>', required=False,
        help='Optional: sample populations. 2nd column is sample, 7th column is population.')
    req_args.add_argument('-g', '--gene', metavar='<string>', required=True,
        help='Gene name.')
    req_args.add_argument('-m', '--method', metavar='<string>', required=False,
        help='Method name. Sometimes required, depends on the gene.')
    req_args.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<file>', required=True,
        help='Output csv file.')

    oth_args = parser.add_argument_group('Other arguments')
    oth_args.add_argument('-h', '--help', action='help', help='Show this message and exit.')
    args = parser.parse_args()

    populations = load_populations(args.populations)
    chrom, positions = get_positions(args.gene)
    entries = load_summary(args.a, chrom, positions)
    compare(args.b, entries, populations, args.gene, args.method, args.output)


if __name__ == '__main__':
    main()
