#!/usr/bin/env python3

import argparse
from collections import defaultdict
import glob
import os
import sys
import copy
import re
import csv

import summarize_results
from summarize_results import dget


class Entry:
    def __init__(self, row):
        self.chrom = row['chrom']
        self.start = int(row['start'])
        self.end = int(row['end'])
        self.locus = dget(row, 'region', 'locus')
        self.sample = row['sample']
        self.info = '{}\t{}\t{}\t{}\t{}\t{}'.format(
            dget(row, 'copy_num_filter', 'agCN_filter'),
            dget(row, 'copy_num', 'agCN'),
            dget(row, 'copy_num_qual', 'agCN_qual'),
            dget(row, 'paralog_filter', 'psCN_filter'),
            dget(row, 'paralog_copy_num', 'psCN'),
            dget(row, 'paralog_qual', 'psCN_qual'))
        self.region_group = re.search(r'group=([0-9\-]+)', row['info']).group(1)
        reg2 = dget(row, 'other_regions', 'homologous_regions')
        self.ref_cn = reg2.count(',') + bool(reg2 != '*') + 1

    def __lt__(self, oth):
        return (self.start, self.end) < (oth.start, oth.end)

    def __len__(self):
        return self.end - self.start

    def __bool__(self):
        return self.end > self.start

    @classmethod
    def write_two(cls, first, second, dataset1, dataset2, out):
        assert first is None or len(first) != 0
        assert second is None or len(second) != 0
        any_entry = first or second
        assert any_entry

        if first and second:
            assert first.chrom == second.chrom and first.start == second.start and first.end == second.end \
                and first.sample == second.sample

        out.write('{}\t{}\t'.format(dataset1, dataset2))
        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(any_entry.locus, 2 * any_entry.ref_cn,
            any_entry.chrom, any_entry.start, any_entry.end, len(any_entry), any_entry.sample))
        if first:
            out.write('\t{}\t{}'.format(first.info, first.region_group))
        else:
            out.write('\t*' * 7)

        if second:
            out.write('\t{}\t{}'.format(second.info, second.region_group))
        else:
            out.write('\t*' * 7)
        out.write('\n')

    def split_by_pos(self, pos):
        assert self.start < pos < self.end - 1
        first = copy.copy(self)
        first.end = pos
        self.start = pos
        return first


class RegionEntries:
    def __init__(self, dataset, locus, summary_path):
        self.dataset = dataset
        self.locus = locus

        self.empty = True
        if not os.path.exists(summary_path):
            sys.stderr.write('WARN: File {} does not exist\n'.format(summary_path))
            return

        self.entries = defaultdict(list)
        with open(summary_path) as inp:
            fieldnames = None
            for line in inp:
                if line.startswith('##'):
                    continue
                assert line.startswith('#')
                fieldnames = line.lstrip('#').strip().split('\t')
                break

            assert fieldnames is not None
            reader = csv.DictReader(inp, fieldnames=fieldnames, delimiter='\t')
            for row in reader:
                entry = Entry(row)
                self.entries[entry.sample].append(entry)
                self.empty = False

        for entries in self.entries.values():
            entries.sort()

    def compare_with(self, other, sample, out):
        j = 0
        entries_a = self.entries[sample]
        entries_b = other.entries[sample]
        size_b = len(entries_b)

        for entry_a in entries_a:
            while j < size_b and entries_b[j].end <= entry_a.start:
                entry_b = entries_b[j]
                Entry.write_two(None, entry_b, self.dataset, other.dataset, out)
                j += 1

            while True:
                if j == size_b or entry_a.end <= entries_b[j].start:
                    Entry.write_two(entry_a, None, self.dataset, other.dataset, out)
                    break

                entry_b = entries_b[j]
                if entry_b.start < entry_a.start:
                    Entry.write_two(None, entry_b.split_by_pos(entry_a.start), self.dataset, other.dataset, out)

                if entry_a.start < entry_b.start:
                    Entry.write_two(entry_a.split_by_pos(entry_b.start), None, self.dataset, other.dataset, out)

                assert entry_a.start == entry_b.start
                if entry_a.end == entry_b.end:
                    Entry.write_two(entry_a, entry_b, self.dataset, other.dataset, out)
                    j += 1
                    break

                if entry_a.end < entry_b.end:
                    Entry.write_two(entry_a, entry_b.split_by_pos(entry_a.end), self.dataset, other.dataset, out)
                    break

                if entry_a.end > entry_b.end:
                    Entry.write_two(entry_a.split_by_pos(entry_b.end), entry_b, self.dataset, other.dataset, out)
                    j += 1


def load_all(dirnames):
    datasets = []
    for dirname in dirnames:
        if '=' in dirname:
            dataset_name, dirname = dirname.split('=', 1)
        else:
            dataset_name = os.path.basename(dirname.rstrip('/'))

        filenames = summarize_results.get_filenames(dirname)
        dataset = {}
        for locus, filename in filenames:
            entries = RegionEntries(dataset_name, locus, filename)
            if entries.empty:
                continue
            dataset[locus] = entries
        datasets.append(dataset)
    return datasets


def compare_pair(dataset1, dataset2, out):
    regions = sorted(set(dataset1.keys()) & set(dataset2.keys()))
    for region in regions:
        regions1 = dataset1[region]
        regions2 = dataset2[region]
        samples = set(regions1.entries.keys()) & set(regions2.entries.keys())
        for sample in samples:
            regions1.compare_with(regions2, sample, out)


def compare_all(datasets, out):
    for i, dataset1 in enumerate(datasets):
        for dataset2 in datasets[i + 1 :]:
            compare_pair(dataset1, dataset2, out)


def main():
    parser = argparse.ArgumentParser(
        description='Merges results from two datasets with the same set of samples.',
        formatter_class=argparse.RawTextHelpFormatter, add_help=False,
        usage='%(prog)s -i <dir> -o <csv>')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', metavar='<file>', required=True, nargs='+',
        help='Input directories (in format [name=]directory).')
    io_args.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<file>', required=False, default='-',
        help='Output csv file.')

    oth_args = parser.add_argument_group('Other arguments')
    oth_args.add_argument('-h', '--help', action='help', help='Show this message and exit.')
    args = parser.parse_args()

    datasets = load_all(args.input)
    args.output.write('# {}\n'.format(' '.join(sys.argv)))
    args.output.write('dataset1\tdataset2\tlocus\tref_cn\tchrom\tstart\tend\tlength\tsample')
    for i in range(1, 3):
        args.output.write('\tcopy_num_filter{0}\tcopy_num{0}\tcopy_num_qual{0}'.format(i))
        args.output.write('\tparalog_filter{0}\tparalog_copy_num{0}\tparalog_qual{0}\tregion_group{0}'.format(i))
    args.output.write('\n')
    compare_all(datasets, args.output)


if __name__ == '__main__':
    main()
