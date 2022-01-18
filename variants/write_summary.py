#!/usr/bin/env python3

import sys
import gzip
import os
from collections import defaultdict
import numpy as np


def write_summary(filename):
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

    THRESHOLDS = [0, 10, 20, 50, 100]

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
            roc_matrix.append(list(map(float, line)))
    inp.close()
    roc_matrix = np.array(roc_matrix)
    roc_matrix = np.delete(roc_matrix, 3, axis=1)

    ixs = defaultdict(list)
    for thresh in THRESHOLDS:
        curr_ixs = np.where(roc_matrix[:, 0] >= thresh)[0]
        if len(curr_ixs) > 0:
            curr_i = curr_ixs[-1]
            ixs[curr_i].append('≥ {}'.format(thresh))

    true_pos_rate = roc_matrix[:, 1] / n_variants
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
    for arg in sys.argv[1:]:
        write_summary(arg)


if __name__ == '__main__':
    main()