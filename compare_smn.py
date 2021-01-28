#!/usr/bin/env python3

import sys
import re
import argparse
from enum import Enum


class Result(Enum):
    Match = 0
    Partial = 1
    Mismatch = 2
ResultPlural = ['Matches', 'Partial', 'Mismatches']


def compare_entries_mlpa(det_line1, det_line2, mlpa_line):
    try:
        mlpa1, mlpa2 = map(round, map(float, mlpa_line[3:5]))
        det_ploidy1 = int(det_line1[7])
        det_ploidy2 = int(det_line2[7])
    except ValueError:
        return Result.Mismatch

    if det_ploidy1 != mlpa1 or det_ploidy2 != mlpa2:
        return Result.Mismatch
    return Result.Match


def compare_entries_smn_caller(det_line1, det_line2, smnc_line):
    try:
        smnc_x, smnc_y, smnc_z = map(int, smnc_line[3:6])
        det_ploidy1 = int(det_line1[7])
        det_ploidy2 = tuple(map(int, det_line2[10].split(',')))
    except ValueError:
        return Result.Mismatch

    if det_ploidy1 != smnc_x + smnc_y + smnc_z or det_ploidy2 != (smnc_x, smnc_y):
        return Result.Mismatch

    try:
        det_par_ploidy1 = tuple(map(int, det_line1[10].split(',')))
        if det_par_ploidy1 != (smnc_x, smnc_y + smnc_z):
            return Result.Partial
        return Result.Match
    except ValueError:
        return Result.Match


def compare_entries(det_line1, det_line2, tool_line, tool):
    if tool == 'SMN_caller':
        return compare_entries_smn_caller(det_line1, det_line2, tool_line)
    elif tool == 'MLPA':
        return compare_entries_mlpa(det_line1, det_line2, tool_line)
    assert False


def get_sample(line, tool):
    if tool == 'SMN_caller':
        return line[0].split('.')[0]
    elif tool == 'MLPA':
        return line[2]


def main():
    parser = argparse.ArgumentParser(
        description='Add columns with Quick-Mer2 results.',
        formatter_class=argparse.RawTextHelpFormatter, add_help=False,
        usage='%(prog)s -s <file> -q <dir> -o <file>')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-s', '--summary', type=argparse.FileType(), metavar='<file>', required=True,
        help='Homology tools summary.')
    io_args.add_argument('-t', '--tool', choices=['MLPA', 'SMN_caller'], required=True,
        help='Other tool.')
    io_args.add_argument('-r', '--tool-res', type=argparse.FileType(), metavar='<dir>', required=True,
        help='File with another results.')
    io_args.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<file>', required=False, default='-',
        help='Output file.')

    oth_args = parser.add_argument_group('Other arguments')
    oth_args.add_argument('-h', '--help', action='help', help='Show this message and exit.')
    args = parser.parse_args()

    positions = [70_940_000, 70_950_000]
    ploidies = [ {}, {} ]
    for line in args.summary:
        if line.startswith('#'):
            continue

        line = line.strip().split('\t')
        reg_start = int(line[1])
        reg_end = int(line[2])
        reg_name = line[3]
        assert reg_name == 'SMN1'
        sample = line[4]

        for i, pos in enumerate(positions):
            if reg_start <= pos < reg_end:
                ploidies[i][sample] = line

    results = [[] for _ in range(len(Result))]

    next(args.tool_res)
    for line in args.tool_res:
        tool_line = line.strip().split('\t')
        sample = get_sample(tool_line, args.tool)

        det_line1 = ploidies[0].get(sample)
        det_line2 = ploidies[1].get(sample)
        if det_line1 is None and det_line2 is None:
            continue

        res = compare_entries(det_line1, det_line2, tool_line, args.tool)
        results[res.value].append((tool_line, det_line1, det_line2))

    out = args.output
    total = sum(map(len, results))
    out.write('Total entries: {:4}\n'.format(total))
    for i in range(len(Result)):
        out.write('{:12}   {:4} ({:4.1f}%)\n'.format(ResultPlural[i] + ':', len(results[i]),
            len(results[i]) / total * 100))

    for i in range(len(Result)):
        if i == Result.Match.value:
            continue
        out.write('\n')
        out.write('{}:\n'.format(ResultPlural[i].upper()))
        for tool_line, det_line1, det_line2 in results[i]:
            out.write('    {}\n'.format('\t'.join(tool_line)))
            out.write('    {}\n'.format('\t'.join(det_line1) if det_line1 else '*'))
            out.write('    {}\n'.format('\t'.join(det_line2) if det_line2 else '*'))
            out.write('------\n')

if __name__ == '__main__':
    main()
