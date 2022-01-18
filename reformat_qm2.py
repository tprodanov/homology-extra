#!/usr/bin/env python3

import pysam
import argparse
import os
import collections
import operator
import numpy as np
import sys
import multiprocessing


Center = collections.namedtuple('Center', 'chrom start end locus paralog')


def load_gene_centers(inp):
    centers = []
    for line in inp:
        if line.startswith('#'):
            continue
        chrom, start, end, locus, paralog = line.strip().split('\t')
        start = int(start)
        end = int(end)
        paralog = int(paralog) - 1
        centers.append(Center(chrom, start, end, locus, paralog))
    return centers


def estimate_quality(cn):
    dist_to_int = abs(cn - np.round(cn))
    assert dist_to_int <= 0.5
    return '{:.1f}'.format(100 - 200 * dist_to_int)


def process_file(filename, gene_centers, n_paralogs):
    sample = os.path.basename(filename).split('.')[0]
    copy_nums = { locus: np.full(count, np.nan) for locus, count in n_paralogs.items() }
    with pysam.TabixFile(filename, parser=pysam.asTuple()) as inp:
        for center in gene_centers:
            entries = list(inp.fetch(center.chrom, center.start, center.end))
            if len(entries) != 1:
                raise ValueError(
                    f'Error at sample {sample}, locus {center.locus}, paralog {center.paralog}: {len(entries)} entries.')
            entry = entries[0]
            copy_nums[center.locus][center.paralog] = float(entry[3])

    res_lines = ''
    for center in gene_centers:
        if center.paralog != 0:
            continue

        pscn = copy_nums[center.locus]
        out_line = f'{center.chrom}\t{center.start}\t{center.end}\t{center.locus}\t{sample}\t'

        agcn = np.sum(pscn)
        out_line += 'PASS\t{:.0f}\t{}\t'.format(np.round(agcn), estimate_quality(agcn))

        pscn_int = np.round(pscn).astype(dtype=np.int16)
        out_line += 'PASS\t{}\t{}\t'.format(','.join(map(str, pscn_int)), ','.join(map(estimate_quality, pscn)))

        reg2 = ','.join('!' * (len(pscn) - 1))
        out_line += f'*\t{reg2}\t*\t*\t*\t{center.locus}\t'
        out_line += '{} -> {:.3f}\n'.format(' '.join(map('{:.3f}'.format, pscn)), agcn)
        res_lines += out_line
    return res_lines


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-I', '--input-list', type=argparse.FileType(), required=True)
    parser.add_argument('-g', '--gene-centers', type=argparse.FileType(), required=True)
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), required=True)
    parser.add_argument('-@', '--threads', type=int, required=False, default=4)
    args = parser.parse_args()

    centers = load_gene_centers(args.gene_centers)
    n_paralogs = collections.Counter(map(operator.attrgetter('locus'), centers))
    filenames = list(map(str.strip, args.input_list))
    n = len(filenames)

    pool = multiprocessing.Pool(args.threads)

    res = []
    def callback(lines):
        res.append(lines)
        if len(res) % 10 == 0:
            sys.stderr.write(f'Finished {len(res):>4} / {n} jobs.\n')
            sys.stderr.flush()

    def err_callback(exc):
        sys.stderr.write('Thread finished with an exception: {}'.format(exc))
        pool.terminate()

    for filename in filenames:
        pool.apply_async(process_file, args=(filename, centers, n_paralogs),
            callback=callback, error_callback=err_callback)
    pool.close()
    pool.join()
    print('Finished:', len(res))

    for s in res:
        args.output.write(s)


if __name__ == '__main__':
    main()
