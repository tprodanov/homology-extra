#!/usr/bin/env python3

import sys
import numpy as np
from Bio import SeqIO


def true_runs(a):
    assert a.dtype == np.bool
    padded = np.hstack(([False], a, [False]))
    diff = np.diff(padded)
    ranges = np.where(diff)[0].reshape(-1, 2)
    return ranges


def main():
    fasta_file = sys.argv[1]
    out = open(sys.argv[2], 'w') if len(sys.argv) > 2 else sys.stdout

    for chrom in SeqIO.parse(fasta_file, 'fasta'):
        name = chrom.id
        seq = str(chrom.seq).encode()
        seq = np.frombuffer(seq, dtype=np.uint8)
        mask = np.logical_or(seq == ord('N'), seq >= ord('a'))
        ranges = true_runs(mask)
        for start, end in ranges:
            out.write('{}\t{}\t{}\n'.format(name, start, end))


if __name__ == '__main__':
    main()

