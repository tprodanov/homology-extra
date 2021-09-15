#!/usr/bin/env python3

import argparse
from collections import defaultdict
import glob
import os
import sys

REGION_COL = 3
SAMPLE_COL = 4
N_COLS = 13


class Entry:
    def __init__(self, line):
        self.cols = line.strip().split('\t')
        self.cols[1] = int(self.cols[1])
        self.cols[2] = int(self.cols[2])

    @property
    def start(self):
        return self.cols[1]

    @start.setter
    def start(self, value):
        self.cols[1] = value

    @property
    def end(self):
        return self.cols[2]

    @end.setter
    def end(self, value):
        self.cols[2] = value

    @property
    def sample(self):
        return self.cols[SAMPLE_COL]

    @property
    def region(self):
        return self.cols[REGION_COL]

    def __getitem__(self, i):
        return self.cols[i]

    def __lt__(self, oth):
        return (self.start, self.end) < (oth.start, oth.end)

    def __len__(self):
        return self.end - self.start

    def __bool__(self):
        return self.end > self.start

    def __str__(self):
        return '\t'.join(map(str, self.cols))

    @classmethod
    def write_two(cls, first, second, out):
        assert first is None or first
        assert second is None or second
        assert first or second

        PREFIX_COLS = list(range(3))
        INFO_COLS = list(range(6, 12))
        EMPTY_SYMBOL = '\t*'
        has_first = first is not None
        has_second = second is not None

        for i in PREFIX_COLS:
            out.write('{}\t'.format(first[i] if has_first else second[i]))
        out.write('{}\t'.format(first.sample if has_first else second.sample))
        out.write('{}\t'.format(first.region if has_first else '*'))
        out.write(second.region if has_second else '*')

        if has_first:
            for i in INFO_COLS:
                out.write('\t{}'.format(first[i]))
        else:
            out.write(EMPTY_SYMBOL * len(INFO_COLS))

        if has_second:
            for i in INFO_COLS:
                out.write('\t{}'.format(second[i]))
        else:
            out.write(EMPTY_SYMBOL * len(INFO_COLS))
        out.write('\n')

    def copy(self):
        cls = self.__class__
        new = cls.__new__(cls)
        new.cols = list(self.cols)
        return new

    def split_by_pos(self, pos):
        assert self.start < pos < self.end - 1
        first = self.copy()
        first.end = pos
        self.start = pos
        return first


def load(dirname, by_sample):
    summary_path = os.path.join(dirname, 'summary.bed')
    if not os.path.exists(summary_path):
        sys.stderr.write('WARN: File {} does not exist\n'.format(summary_path))
        return
    with open(summary_path) as inp:
        for line in inp:
            if line.startswith('#'):
                continue
            entry = Entry(line)
            by_sample[entry.sample].append(entry)


def load_multiple(dirnames):
    by_sample = defaultdict(list)
    for dirname in dirnames:
        load(dirname, by_sample)
    for value in by_sample.values():
        value.sort()
    return by_sample


def compare_sample_entries(entries_a, entries_b, out):
    j = 0
    size_b = len(entries_b)

    for entry_a in entries_a:
        while j < size_b and entries_b[j].end <= entry_a.start:
            entry_b = entries_b[j]
            Entry.write_two(None, entry_b, out)
            j += 1

        while True:
            if j == size_b or entry_a.end <= entries_b[j].start:
                Entry.write_two(entry_a, None, out)
                break

            entry_b = entries_b[j]
            if entry_b.start < entry_a.start:
                Entry.write_two(None, entry_b.split_by_pos(entry_a.start), out)

            if entry_a.start < entry_b.start:
                Entry.write_two(entry_a.split_by_pos(entry_b.start), None, out)

            assert entry_a.start == entry_b.start
            if entry_a.end == entry_b.end:
                Entry.write_two(entry_a, entry_b, out)
                j += 1
                break

            if entry_a.end < entry_b.end:
                Entry.write_two(entry_a, entry_b.split_by_pos(entry_a.end), out)
                break

            if entry_a.end > entry_b.end:
                Entry.write_two(entry_a.split_by_pos(entry_b.end), entry_b, out)
                j += 1


def analyze_groups(dirnames_a, dirnames_b, out):
    by_sample_a = load_multiple(dirnames_a)
    by_sample_b = load_multiple(dirnames_b)

    samples = sorted(set(by_sample_a) | set(by_sample_b))
    for sample in samples:
        compare_sample_entries(by_sample_a[sample], by_sample_b[sample], out)


def find_all_groups(path):
    groups = defaultdict(lambda: ([], []))
    for dirname in glob.glob(os.path.join(path, "*-*", "")):
        basename = os.path.basename(os.path.dirname(dirname))
        gene, suffix = basename.split('-')
        groups[gene][suffix == 'full'].append(dirname)
    return groups


def main():
    parser = argparse.ArgumentParser(
        description='Compare results with splitted and non-splitted runs.',
        formatter_class=argparse.RawTextHelpFormatter, add_help=False,
        usage='%(prog)s -i <dir> -o <csv>')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', metavar='<file>', required=True,
        help='Input directory with directories <gene_name>-<suffix>.')
    io_args.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<file>', required=False, default='-',
        help='Output csv file.')

    oth_args = parser.add_argument_group('Other arguments')
    oth_args.add_argument('-h', '--help', action='help', help='Show this message and exit.')
    args = parser.parse_args()

    groups = find_all_groups(args.input)
    for dirnames_a, dirnames_b in groups.values():
        analyze_groups(dirnames_a, dirnames_b, args.output)


if __name__ == '__main__':
    main()
