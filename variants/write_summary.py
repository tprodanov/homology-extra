#!/usr/bin/env python3

import sys
import gzip
import os
from collections import defaultdict
import numpy as np
import argparse


def float_or_nan(s):
    if s == 'None':
        return np.nan
    return float(s)


def write_summary(filename, thresholds, only_thresholds):
    if os.path.isfile(filename):
        dirname = os.path.dirname(filename)
    else:
        dirname = filename
        filename = os.path.join(dirname, 'weighted_roc.tsv.gz')

    print(f'====  {dirname}  ====')
    if filename.endswith('.gz'):
        inp = gzip.open(filename, 'rt')
    else:
        inp = open(filename)

    n_variants = None
    roc_matrix = []
    for line in inp:
        line = line.strip().split()
        if line[0].startswith('#'):
            if line[1] == 'baseline':
                print(f'Total baseline:      {line[-1]}')
                n_variants = int(line[-1])
            elif line[1] == 'call':
                print(f'Total call variants: {line[-1]}')
            elif line[1] == 'true_positives_baseline':
                print('Threshold  True-pos  False-pos  False-neg  Precision  Recall  F1-score')
                print('----------------------------------------------------------------------')
        else:
            roc_matrix.append(list(map(float_or_nan, line)))
    inp.close()
    roc_matrix = np.array(roc_matrix)
    roc_matrix = np.delete(roc_matrix, 3, axis=1)

    ixs = defaultdict(list)
    for thresh in thresholds:
        curr_ixs = np.where(roc_matrix[:, 0] >= thresh)[0]
        if len(curr_ixs) > 0:
            curr_i = curr_ixs[-1]
            ixs[curr_i].append('≥ {:.0f}'.format(thresh))

    if roc_matrix.shape[0] == 1 and np.isnan(roc_matrix[0, 0]):
        roc_matrix[0, 0] = 0

    true_pos_rate = roc_matrix[:, 1] / n_variants
    if not only_thresholds:
        for col, flag in [(4, 'best precision'), (5, 'best recall'), (6, 'best F1')]:
            curr_i = np.argmax(roc_matrix[:, col] + true_pos_rate * 0.0001)
            ixs[curr_i].append(flag)

    for i, flags in sorted(ixs.items()):
        s = '{:9.3f}  {:8.0f}  {:9.0f}  {:9.0f}  {:9.4f}  {:6.4f}  {:8.4f}'.format(*roc_matrix[i])
        if flags:
            s += '  ({})'.format(', '.join(flags))
        print(s)
    print()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs='+', metavar='FILE', help='Input directories')
    parser.add_argument('-t', '--thresholds', metavar='FLOAT,FLOAT...', default='0,10,20,50,100',
        help='Quality thresholds, through comma [default %(default)s].')
    parser.add_argument('-T', '--only-thresholds', action='store_true',
        help='Print only accuracy values for quality thresholds.')
    args = parser.parse_args()

    thresholds = list(map(float, args.thresholds.split(',')))
    for path in args.input:
        write_summary(path, thresholds, args.only_thresholds)


if __name__ == '__main__':
    main()
