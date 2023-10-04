#!/usr/bin/env python3

import sys
import argparse
import collections
from tqdm import tqdm
from intervaltree import Interval, IntervalTree


class Stats:
    def __init__(self):
        self.kept = 0
        self.non_dupl = 0
        self.high_cn = 0
        self.tangled = 0
        self.too_diverse = 0
        self.too_short = 0
        self.duplicates = 0
        self.too_short_after = 0

    def discarded(self):
        return self.non_dupl + self.high_cn + self.tangled + self.too_short + self.too_diverse \
            + self.duplicates + self.too_short_after

    def summarize(self):
        sys.stderr.write(f'Kept            {self.kept:11,} bp\n')
        sys.stderr.write(f'Discarded       {self.discarded():11,} bp:\n')
        sys.stderr.write(f'   non dupl.:   {self.non_dupl:11,} bp\n')
        sys.stderr.write(f'   high CN:     {self.high_cn:11,} bp\n')
        sys.stderr.write(f'   tangled:     {self.tangled:11,} bp\n')
        sys.stderr.write(f'   too diverse: {self.too_diverse:11,} bp\n')
        sys.stderr.write(f'   too short:   {self.too_short:11,} bp\n')
        sys.stderr.write(f'   duplicates:  {self.duplicates:11,} bp\n')
        sys.stderr.write(f'   too short:   {self.too_short_after:11,} bp (after duplicate removal)\n')


def process(line, trees, stats, out, args):
    chrom, start, end, filt, cn, info, homol = line.strip().split('\t')
    start = int(start)
    end = int(end)
    length = end - start
    if length < args.length:
        stats.too_short += length
        return
    elif filt != 'PASS':
        stats.tangled += length
        return
    elif cn.startswith('>='):
        stats.high_cn += length
    cn = int(cn)
    if cn < 3:
        stats.non_dupl += length
        return
    elif cn > args.max_cn:
        stats.high_cn += length
        return

    max_simil = max(map(float, info.rsplit('=', 1)[1].split(',')))
    assert max_simil <= 1
    if max_simil < args.simil:
        stats.too_diverse += length
        return

    for overlap in sorted(trees[chrom][start:end]):
        curr_end = overlap.begin
        assert curr_end <= end
        curr_len = max(0, curr_end - start)
        if curr_len < args.length:
            stats.too_short_after += curr_len
        else:
            out.write(f'{chrom}\t{start}\t{curr_end}\n')
            stats.kept += curr_len
        stats.duplicates += max(0, min(end, overlap.end) - max(start, overlap.begin))
        start = min(end, overlap.end)

    rem_len = max(0, end - start)
    if rem_len < args.length:
        stats.too_short_after += rem_len
    else:
        out.write(f'{chrom}\t{start}\t{end}\n')
        stats.kept += rem_len

    for region in homol.split(','):
        oth_chrom, oth_coord, _ = region.split(':')
        oth_start, oth_end = map(int, oth_coord.split('-'))
        trees[oth_chrom][oth_start - 1:oth_end] = ()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True,
        help='Input BED file.')
    parser.add_argument('-o', '--output', required=True,
        help='Output BED file.')
    parser.add_argument('-l', '--length', type=int, default=3000,
        help='Minimal duplication length [%(default)s].')
    parser.add_argument('-s', '--simil', type=float, default=0.97,
        help='Minimal sequence similarity [%(default)s].')
    parser.add_argument('-c', '--max-cn', type=int, default=8,
        help='Maximum copy number [%(default)s].')
    args = parser.parse_args()
    assert args.length > 0 # will not work without it.

    trees = collections.defaultdict(IntervalTree)
    stats = Stats()
    with open(args.input) as inp, open(args.output, 'w') as out:
        for line in (inp):
            if line.startswith('#'):
                continue
            process(line, trees, stats, out, args)
    stats.summarize()


if __name__ == '__main__':
    main()
