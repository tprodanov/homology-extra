#!/usr/bin/env python3

import collections


class Locus:
    def __init__(self, line):
        self.chrom, start, end, self.names = line.strip().split('\t')[:4]
        self.start, self.end = map(int, (start, end))

    def __str__(self):
        return f'{self.chrom}\t{self.start}\t{self.end}\t{self.names}'

    def __len__(self):
        return self.end - self.start

    def overlap(self, oth):
        assert self.chrom == oth.chrom
        return max(min(self.end, oth.end) - max(self.start, oth.start), 0)

    def frac_overlap(self, oth):
        return self.overlap(oth) / min(len(self), len(oth))


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Merge overlapping genes.',
        usage='%(prog)s -i input.bed [-f NUM] -o output.bed')
    parser.add_argument('-i', '--input', metavar='FILE', default='/dev/stdin',
        help='Input BED flie, must be sorted!')
    parser.add_argument('-f', '--fraction', metavar='NUM', type=float, default=0.5,
        help='Merge loci, contained in another by this fraction.')
    parser.add_argument('-o', '--output', metavar='DIR', default='/dev/stdout',
        help='Output BED flie.')
    args = parser.parse_args()

    out = open(args.output, 'w')
    current = collections.deque()
    for line in open(args.input):
        locus = Locus(line)
        size = len(locus)
        while len(current) > 0 and (current[0].chrom != locus.chrom or current[0].end <= locus.start):
            out.write(f'{current.popleft()}\n')

        for curr in current:
            if curr.frac_overlap(locus) >= args.fraction:
                curr.end = locus.end
                curr.names += f',{locus.names}'
                break
        else:
            current.append(locus)
    for locus in current:
        out.write(f'{locus}\n')


if __name__ == '__main__':
    main()
