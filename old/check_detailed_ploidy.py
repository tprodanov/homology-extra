#!/usr/bin/env python3

import argparse
import numpy as np
import glob
import os
import csv


class DetailedPloidy:
    def __init__(self, inp):
        fieldnames = None
        while fieldnames is None:
            line = next(inp)
            if line.startswith('##'):
                continue
            assert line.startswith('#')
            fieldnames = line.lstrip('#').strip().split('\t')
        self.fieldnames = { name: i for i, name in enumerate(fieldnames) }

        self.starts = []
        self.ends = []
        self.lines = []
        for line in inp:
            line = line.strip().split('\t')
            self.starts.append(int(line[1]))
            self.ends.append(int(line[2]))
            self.lines.append(line)
        self.starts = np.array(self.starts)
        self.ends = np.array(self.ends)

    def get_entries(self, sample, start, end):
        start_ix = self.ends.searchsorted(start, side='right')
        end_ix = self.starts.searchsorted(end, side='left')
        col = self.fieldnames[sample]
        return [self.lines[i][col] for i in range(start_ix, end_ix)]


def load_detailed_ploidies(dir):
    res = {}
    for filename in glob.glob(os.path.join(dir, '*/extra/detailed_ploidy.bed')):
        region = filename.split('/')[-3]
        with open(filename) as inp:
            res[region] = DetailedPloidy(inp)
    return res


def iterate_comparison(reader, length, qual, all_filt):
    for row in reader:
        if int(row['length']) < length:
            continue

        filt1 = row['ploidy_filter1']
        filt2 = row['ploidy_filter2']
        if filt1 == '*' or filt2 == '*':
            continue
        if not all_filt and (filt1 != 'PASS' or filt2 != 'PASS'):
            continue
        qual1 = float(row['ploidy_qual1'])
        qual2 = float(row['ploidy_qual2'])
        if qual1 < qual or qual2 < qual:
            continue
        ploidy1 = row['ploidy1']
        ploidy2 = row['ploidy2']
        if ploidy1 != ploidy2:
            yield row


def analyze_entries(entries, res):
    counts = [0] * 3
    for entry in entries:
        p1, p2 = entry.split()
        prob = 10 ** -float(p2.split(':')[1])
        counts[0] += 1
        suffix = ''
        if prob < 0.001:
            suffix = ' !!!'
            counts[2] += 1
        elif prob < 0.01:
            suffix = ' !'
            counts[1] += 1
        res += '        {}{}\n'.format(entry, suffix)
    return counts, res


def analyze_discrepancy(row, detailed_ploidies, discr_only, out):
    res = ''
    res += '{}\n'.format('\t'.join(row.values()))
    dataset1 = row['dataset1']
    dataset2 = row['dataset2']
    start = int(row['start'])
    end = int(row['end'])
    region = row['region']
    sample = row['sample']

    res += '=== {}   Length {:,} bp,   sample {}\n'.format(region, end - start, sample)
    res += '    {:11} {} ({})\n'.format(dataset1 + ':', row['ploidy1'], row['ploidy_qual1'])
    res += '    {:11} {} ({})\n'.format(dataset2 + ':', row['ploidy2'], row['ploidy_qual2'])
    res += '------------------\n'

    res += '    {}:\n'.format(dataset1)
    counts1, res = analyze_entries(detailed_ploidies[dataset1][region].get_entries(sample, start, end), res)
    res += '    {}:\n'.format(dataset2)
    counts2, res = analyze_entries(detailed_ploidies[dataset2][region].get_entries(sample, start, end), res)

    has_discr = counts1[1] + counts1[2] > 0 and counts2[1] + counts2[2] > 0
    if has_discr:
        res += '   Discrepancy!\n'
    res += '\n==================\n'
    if not discr_only or has_discr:
        out.write(res)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input directores.', required=True, nargs='+')
    parser.add_argument('-c', '--comparison', help='Comparison file.', type=argparse.FileType(), required=True)
    parser.add_argument('-q', '--quality', metavar='<float>', type=float, default=30,
        help='Ploidy quality threshold [%(default)s].')
    parser.add_argument('-l', '--min-len', metavar='<int>', type=int, default=0,
        help='Minimal region length [%(default)s].')
    parser.add_argument('-f', '--ignore-filters', action='store_true',
        help='Look at all records, not only with PASS filter.')
    parser.add_argument('-d', '--discrepancy', action='store_true',
        help='Write only discrepancies.')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default='-')
    args = parser.parse_args()

    detailed_ploidies = {}
    for dataset in args.input:
        detailed_ploidies[dataset] = load_detailed_ploidies(dataset)

    reader = csv.DictReader(args.comparison, delimiter='\t')
    for row in iterate_comparison(reader, args.min_len, args.quality, args.ignore_filters):
        analyze_discrepancy(row, detailed_ploidies, args.discrepancy, args.output)


if __name__ == '__main__':
    main()
