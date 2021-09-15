#!/usr/bin/env python3

import argparse
import csv
import sys
import collections


class Concordance:
    def __init__(self):
        self.concordant = 0
        self.total = 0

    def add(self, row):
        length = int(row['end']) - int(row['start'])
        self.total += length
        if row['ploidy1'] == row['ploidy2']:
            self.concordant += length

    def __str__(self):
        return '{}\t{}'.format(self.concordant, self.total)

    def __repr__(self):
        return '(conc={}, total={})'.format(self.concordant, self.total)


def calculate_concordance(inp, qual):
    # Skip comment line and check if it is really a comment.
    assert next(inp).startswith('#')
    reader = csv.DictReader(inp, delimiter='\t')
    concordance = collections.defaultdict(Concordance)

    for row in reader:
        if row['ploidy1'] == '*' or row['ploidy2'] == '*' \
                or float(row['ploidy_qual1']) < qual or float(row['ploidy_qual2']) < qual:
            continue

        if row['region_group1'] != row['region_group2']:
            print('Entry has different region groups:\n    {}'.format('\t'.join(row)))
        key = (row['dataset1'], row['dataset2'], row['region'], row['sample'], row['region_group1'])
        concordance[key].add(row)
    return concordance


def write_concordance(concordance, out):
    out.write('# {}\n'.format(' '.join(sys.argv)))
    out.write('dataset1\tdataset2\tgene\tsample\tregion_group\tconcordant\ttotal\n')
    for key in sorted(concordance.keys()):
        out.write('\t'.join(key))
        out.write('\t{}\n'.format(concordance[key]))


def main():
    parser = argparse.ArgumentParser(
        description='Calculate concordance between different runs.',
        formatter_class=argparse.RawTextHelpFormatter, add_help=False,
        usage='%(prog)s -i <csv> -o <csv>')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', type=argparse.FileType(), metavar='<file>', required=True,
        help='Input CSV file.')
    io_args.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<file>', required=True,
        help='Output CSV file.')

    opt_args = parser.add_argument_group('Optional arguments')
    opt_args.add_argument('-q', '--qual', type=float, metavar='<float>', default=30,
        help='Threshold quality [default: %(default)s].')

    oth_args = parser.add_argument_group('Other arguments')
    oth_args.add_argument('-h', '--help', action='help', help='Show this message and exit.')
    args = parser.parse_args()

    concordance = calculate_concordance(args.input, args.qual)
    write_concordance(concordance, args.output)


if __name__ == '__main__':
    main()
