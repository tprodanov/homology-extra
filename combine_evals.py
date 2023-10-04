#!/usr/bin/env python3

import sys
import os
import argparse
from multiprocessing import Pool
import functools
from tqdm import tqdm
import re
import gzip


def open_stream(filename, mode='r'):
    assert mode == 'r' or mode == 'w'
    if filename is None or filename == '-':
        return sys.stdin if mode == 'r' else sys.stdout
    elif filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    else:
        return open(filename, mode)


def _recursive_find_tuples(path, matches, all_tuples, curr_tuple=(), depth=0, shift=0):
    if depth == len(matches):
        all_tuples.append(curr_tuple)
        return 0
    m = matches[depth]
    dir = path[ : m.start() + shift]
    if not os.path.isdir(dir):
        return 1

    suffix = path[m.end() + shift : ]
    subdirs = []
    with os.scandir(dir) as it:
        for entry in it:
            if not entry.name.startswith('.') and entry.is_dir():
                subdirs.append(entry)

    ignored = 0
    for entry in subdirs:
        ignored += _recursive_find_tuples(entry.path + suffix, matches, all_tuples,
            curr_tuple + (entry.name,), depth + 1, len(entry.path) - m.end())
    return ignored


def load_tags(path):
    path = os.path.abspath(path)
    matches = list(re.finditer(r'\{([a-zA-Z0-9_]+)\}', path))
    tags = [m.group(1) for m in matches]
    all_tuples = []
    ignored = _recursive_find_tuples(path, matches, all_tuples)
    sys.stderr.write('Found {} tag combinations, discarded {}.\n'.format(len(all_tuples), ignored))
    return tags, all_tuples


def load_and_process(tup, tags, input_fmt, quality):
    d = dict(zip(tags, tup))
    input_dir = input_fmt.format(**d)

    total_regions = 0
    sum_len = 0
    with open(os.path.join(input_dir, '..', 'comparison.bed')) as f:
        for line in f:
            line = line.strip()
            if line:
                total_regions += 1
                line = line.split('\t')
                sum_len += int(line[2]) - int(line[1])

    s = ''
    for ty, filename in [('all', 'weighted_roc.tsv.gz'), ('snp', 'snp_roc.tsv.gz'), ('indels', 'non_snp_roc.tsv.gz')]:
        total_call = None
        total_base = None
        with gzip.open(os.path.join(input_dir, filename), 'rt') as f:
            last_line = None
            for line in f:
                if line.startswith('#total baseline variants:'):
                    total_base = int(line.split(':')[1])
                elif line.startswith('#total call variants:'):
                    total_call = int(line.split(':')[1])
                elif not line.startswith('#'):
                    line = list(map(float, line.strip().split('\t')))
                    if line[0] >= quality:
                        last_line = line
                    else:
                        break

        for tag in tup:
            s += f'{tag}\t'
        s += f'{ty}\t{sum_len}\t{total_regions}\t{total_base}\t{total_call}\t'
        if last_line:
            s += '\t'.join(map('{:.5g}'.format, last_line))
        else:
            s += '\t'.join(('NA',) * 8)
        s += '\n'
    return s


def main():
    parser = argparse.ArgumentParser(
        description='Summarize evaluation result.',
        usage='%(prog)s -i path -o out.csv [-@ threads]')
    parser.add_argument('-i', '--input', metavar='STR', required=True,
        help='Path to evaluation results. Tags within `{tag}` are automatically found.')
    parser.add_argument('-o', '--output', metavar='FILE',  required=True,
        help='Output CSV file.')
    parser.add_argument('-q', '--quality', metavar='FLOAT', type=float, default=20,
        help='Quality threshold [%(default)s].')
    parser.add_argument('-@', '--threads', metavar='INT', type=int, default=8,
        help='Number of threads [%(default)s].')
    args = parser.parse_args()

    tags, tag_tuples = load_tags(args.input)

    out = open_stream(args.output, 'w')
    out.write('# {}\n'.format(' '.join(sys.argv)))
    for tag in tags:
        out.write(tag + '\t')
    out.write('type\tsum_len\tn_regions\tbase_vars\tcall_vars\tqual\ttp_base\tfp\ttp_call\tfn\tprecision\trecall\tf1\n')

    n = len(tag_tuples)
    threads = min(n, args.threads)
    f = functools.partial(load_and_process, tags=tags, input_fmt=args.input, quality=args.quality)

    if threads > 1:
        with Pool(threads) as pool:
            for line in tqdm(pool.imap(f, tag_tuples), total=n):
                out.write(line)
    else:
        for line in tqdm(map(f, tag_tuples), total=n):
            out.write(line)


if __name__ == '__main__':
    main()
