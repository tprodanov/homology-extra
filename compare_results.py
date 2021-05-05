#!/usr/bin/env python3

import argparse
import sys
import numpy as np
from collections import defaultdict
import csv
import operator


class SummaryEntry:
    def __init__(self, row):
        self.row = row
        self.chrom = row['chrom']
        self.start = int(row['start'])
        self.end = int(row['end'])
        self.sample = row['sample']

        self.filter = row['copy_num_filter']
        self.copy_num = row['copy_num']
        self.qual = row['copy_num_qual']

        self.paralog_filter = row['paralog_filter']
        self.paralog_copy_num = row['paralog_copy_num']
        self.paralog_qual = row['paralog_qual']

    def covers(self, chrom, pos):
        return self.chrom == chrom and self.start <= pos <= self.end


def load_summary(inp, chrom, positions):
    summaries = defaultdict(dict)
    fieldnames = None
    for line in inp:
        if not line.startswith('##'):
            assert line.startswith('#')
            fieldnames = line.lstrip('#').strip().split('\t')
            break

    assert fieldnames is not None
    reader = csv.DictReader(inp, fieldnames, delimiter='\t')
    for row in reader:
        entry = SummaryEntry(row)
        for key, pos in positions.items():
            if entry.covers(chrom, pos):
                summaries[entry.sample][key] = entry
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
        pos = { 'SERF1A_middle': 70_902_000,
                'SERF1A_end': 70_918_500,
                'SMN1_middle': 70_940_000,
                'SMN1_end': 70_950_000 }
    elif gene == 'NPY4R':
        chrom = 'chr10'
        pos = (46_462_500,)
    elif gene == 'RHCE':
        chrom = 'chr1'
        pos = (25_402_700,)
    elif gene == 'PMS2':
        chrom = 'chr7'
        pos = { 'exon15': 5_973_500,
                'exon14': 5_977_650 }
    elif gene == 'C4A':
        chrom = 'chr6'
        pos = { 'start': 31_982_100, # Exon 1
                'middle': 31_988_000, # Intron 9
                'end': 31_995_000 } # Exon 22
    elif gene == 'FAM185A':
        chrom = 'chr7'
        pos = { 'before': 102_790_000,
                'del': 102_795_000,
                'after': 102_800_000 }
    elif gene == 'NBPF4':
        chrom = 'chr1'
        pos = { 'start': 108_229_000, # Exon 13
                'middle': 108_240_000, # Exon 5
                'end': 108_244_000 } # Exon 1
    elif gene == 'EIF3C':
        chrom = 'chr16'
        pos = { 'start': 28_690_000,
                'middle': 28_705_000,
                'end': 28_725_000 }
    elif gene == 'SRGAP2':
        chrom = 'chr1'
        pos = { 'a': 206_200_000,
                'b': 206_250_000,
                'c': 206_300_000,
                'd': 206_350_000,
                'e': 206_400_000
        }
    elif gene == 'STRC':
        chrom = 'chr15'
        pos = { 'gene_start': 43_600_000,
                '02-01': 43_580_000,
                '02-02': 43_610_000,
                '02-03': 43_629_000 }
    elif gene == 'NEB':
        chrom = 'chr2'
        pos = { 'middle': 151_584_000 }
    elif gene == 'GTF2I':
        chrom = 'chr7'
        pos = { '0': 74_750_000 }
    elif gene == 'APOBEC3A':
        chrom = 'chr22'
        pos = { '0': 38_961_000 }
    elif gene == 'ABCC6':
        chrom = 'chr16'
        pos = {
            '02-01': 16_210_000,
            '03-01': 16_220_000 }
    elif gene == 'OTOA':
        chrom = 'chr16'
        pos = { '0': 21_750_000 }
    elif gene == 'HYDIN':
        chrom = 'chr16'
        pos = { 'exon18': 71_060_550 }
    elif gene == 'CEL':
        chrom = 'chr9'
        pos = { 'exon9': 133_069_200 }
    elif gene == 'NCF1':
        chrom = 'chr7'
        pos = { 'exon7': 74_783_600 }
    else:
        sys.stderr.write(f'Cannot find gene {gene}\n')
        exit(1)
    if not isinstance(pos, dict):
        pos = dict(enumerate(pos))
    return chrom, pos


def load_populations(inp):
    res = defaultdict(lambda: '*\t*')
    if inp is None:
        return res

    for line in inp:
        line = line.strip().split('\t')
        res[line[0]] = '{}\t{}'.format(line[3], line[5])
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

        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(entry.sample, populations[entry.sample],
            entry.filter, entry.copy_num, entry.qual, entry.paralog_filter, entry.paralog_copy_num, entry.paralog_qual))
        out.write(str(self.method))

        if self.b_copy_num is not None or self.b_paralog is not None:
            copy_num_dist = total_copy_num_distance(entry, self.b_copy_num)
            paralog_dist = paralog_copy_num_distance(entry, self.b_paralog)
            out.write('\t{}\t{}\t{:.4g}\t{:.4g}'.format(
                b_copy_num_str, b_paralog_str, copy_num_dist, paralog_dist))
        out.write('\n')


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
    sample_entries = entries[sample]
    entry_16 = sample_entries.get('SMN1_middle')
    entry_78 = sample_entries.get('SMN1_end')
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

    sample_entries = entries[sample]
    entry_16 = sample_entries.get('SMN1_middle')
    entry_78 = sample_entries.get('SMN1_end')
    combine_smn_entries(entry_16, entry_78)
    total_16, total_78 = line[3:5]
    yield ResEntry(entry_16, 'mlpa_16', total_16)
    yield ResEntry(entry_78, 'mlpa_78', total_78)


def compare_smn1_qm2(file, entries):
    header = next(file).strip().split('\t')
    smn2 = next(file).strip().split('\t')
    smn1 = next(file).strip().split('\t')

    chrom, positions = get_positions('SMN1')
    assert chrom == smn1[0] and int(smn1[1]) <= positions['SMN1_middle'] < int(smn1[2])

    for i in range(7, len(header)):
        sample = header[i]
        if sample not in entries:
            continue
        sample_entries = entries[sample]
        entry_16 = sample_entries.get('SMN1_middle')
        entry_78 = sample_entries.get('SMN1_end')
        combine_smn_entries(entry_16, entry_78)

        cn1 = float(smn1[i])
        cn2 = float(smn2[i])
        yield ResEntry(entry_16, 'qm2', cn1 + cn2, (cn1, cn2))


def compare_line_npy4r(line, entries):
    line = line.strip('\n').split('\t')
    sample = line[0]
    if sample not in entries:
        return
    entry = entries[sample][0]

    freec, cnvnator, ddpcr = line[2:]
    yield ResEntry(entry, 'FREEC', freec)
    yield ResEntry(entry, 'CNVnator', cnvnator)
    yield ResEntry(entry, 'ddPCR', ddpcr)


def compare_line_rhd(line, entries):
    line = line.strip('\n').split('\t')
    sample = line[0].strip('*')
    if sample not in entries:
        return
    entry = entries[sample][0]

    rhd_1, rhce_1, rh_1, rhd_2, rhce_2 = line[2:]
    rh_1 = rh_1.split()[0]
    yield ResEntry(entry, 'WGS', rh_1, (rhce_1, rhd_1))
    yield ResEntry(entry, 'MIP', int(rhd_2) + int(rhce_2), (rhce_2, rhd_2))


def compare_line_pms2(line, entries):
    line = line.strip('\n').split('\t')
    sample = line[0]
    if sample not in entries:
        return
    entry = entries[sample][0]

    has_deletion = False
    if len(line) > 1 and line[1].strip():
        assert line[1] == 'Exon 13-14 deletion'
        has_deletion = True
    yield ResEntry(entry, 'exon_13_14', 4 - has_deletion)


def _reorder_srgap2(value):
    if ',' not in value:
        return value
    values = value.split(',')
    assert len(values) == 4
    #                SRGAP2A    SRGAP2B    SRGAP2C    SRGAP2D
    return ','.join((values[0], values[3], values[1], values[2]))


def compare_line_srgap2(line, entries):
    if line.startswith('#') or line.startswith('\t') or line.startswith(' '):
        return
    line = line.strip('\n').split('\t')
    sample = line[0].rstrip('*').strip()
    if sample not in entries:
        return

    entry = entries[sample]['e']
    entry.paralog_copy_num = _reorder_srgap2(entry.paralog_copy_num)
    entry.paralog_qual = _reorder_srgap2(entry.paralog_qual)

    yield ResEntry(entry, 'WGS', line[6], tuple(line[2:6]))
    mip_based = line[7:11]
    yield ResEntry(entry, 'MIP', sum(map(int, mip_based)), tuple(mip_based))

    if len(line) > 11 and len(line[11].strip()) > 0:
        fish = line[11:15]
        fish[-1] = fish[-1].split(' (')[0]
        if '2 or 3' in fish[-1]:
            fish[-1] = 2.5
        yield ResEntry(entry, 'FISH', sum(map(float, fish)), tuple(fish))


def get_qm2_order(gene):
    if gene == 'SMN1':
        return { 'SMN1_middle': ('SMN1', 'SMN2'),
                 'SMN1_end': ('SMN1', 'SMN2') }
    elif gene == 'ABCC6':
        return { '03-01': ('ABCC6', 'ABCC6P2', 'ABCC6P1') }
    elif gene == 'GTF2I':
        return { '0': ('GTF2I', 'GTF2IP4', 'GTF2IP1') }
    elif gene == 'CEL':
        return { 'exon9': ('CEL', 'CELP') }
    elif gene == 'C4A':
        return { 'end': ('C4A', 'C4B') }
    elif gene == 'NCF1':
        return { 'exon7': ('NCF1', 'NCF1B', 'NCF1C') }
    else:
        sys.stderr.write(f'Cannot find QuicK-mer2 entries for gene {gene}\n')
        exit(1)


def compare_qm2(file, gene, entries):
    order = get_qm2_order(gene)
    lines = [line.strip().split('\t') for line in file]
    header = lines[0]

    for region, qm2_names in order.items():
        qm2_estimates = []
        for name in qm2_names:
            for line in lines[1:]:
                if name == line[3]:
                    qm2_estimates.append(line)
                    break
            else:
                sys.stderr.write(f'Cannot find QuicK-mer2 entry {name}\n')
                exit(1)

        for i in range(7, len(header)):
            sample = header[i]
            if sample not in entries:
                continue
            sample_entries = entries[sample]
            entry = sample_entries[region]

            cns = tuple(float(curr_estimates[i]) for curr_estimates in qm2_estimates)
            yield ResEntry(entry, region, '{:.5f}'.format(sum(cns)), cns)


def output_values(sample, entries):
    for key, entry in entries[sample].items():
        yield ResEntry(entry, key, None)


def select_function(gene, method):
    if method is None:
        return output_values

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
        return compare_smn1_qm2
    elif gene == 'NPY4R' and method == 'exper':
        return compare_line_npy4r
    elif gene == 'RHCE' and method == 'exper':
        return compare_line_rhd
    elif gene == 'PMS2' and method == 'exper':
        return compare_line_pms2
    elif gene == 'SRGAP2' and method == 'exper':
        return compare_line_srgap2
    elif gene == 'ABCC6' and method == 'qm2':
        return compare_abcc6_qm2
    elif gene == 'GTF2I' and method == 'qm2':
        return compare_gtf2i_qm2
    else:
        sys.stderr.write(f'Cannot find gene {gene} and method {method}\n')
        exit(1)


def write_header(chrom, positions, gene, method, out):
    out.write('# {}\n'.format(' '.join(sys.argv)))
    if method is None:
        for key, value in sorted(positions.items(), key=operator.itemgetter(1)):
            out.write('# {:15}   {}:{:,}\n'.format(key, chrom, value))
    out.write('sample\tpopulation\tsuperpopulation\tcopy_num_filter\tcopy_num\tcopy_num_qual\t'
        'paralog_filter\tparalog_copy_num\tparalog_qual\t')
    if method is None or method == 'qm2':
        out.write('region')
    else:
        out.write('method')
    if method is not None:
        out.write('\tb_copy_num\tb_paralog\ttotal_dist\tparalog_dist')
    out.write('\n')


def compare(b_in, entries, populations, gene, method, out):
    if method == 'qm2':
        for res in compare_qm2(b_in, gene, entries):
            res.write_to(out, populations)
        return

    if method is None:
        b_in = list(entries.keys())

    compare_fn = select_function(gene, method)
    for line in b_in:
        for res in compare_fn(line, entries):
            res.write_to(out, populations)


def main():
    parser = argparse.ArgumentParser(
        description='Compare results with another method.',
        formatter_class=argparse.RawTextHelpFormatter, add_help=False,
        usage='%(prog)s -a <file> -b <file> [-p <file>] -g <gene> [-m <method>] -o <file>')
    req_args = parser.add_argument_group('Required arguments')
    req_args.add_argument('-a', type=argparse.FileType(), metavar='<file>', required=True, nargs='+',
        help='Homology tools summary (or summaries).')
    req_args.add_argument('-b', type=argparse.FileType(), metavar='<file>',
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
    assert args.method is None or args.b is not None

    populations = load_populations(args.populations)
    chrom, positions = get_positions(args.gene)
    write_header(chrom, positions, args.gene, args.method, args.output)

    entries = {}
    for input in args.a:
        entries.update(load_summary(input, chrom, positions))
    compare(args.b, entries, populations, args.gene, args.method, args.output)


if __name__ == '__main__':
    main()
