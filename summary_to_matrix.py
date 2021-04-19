#!/usr/bin/env python3

import argparse
import numpy as np
import itertools
import collections
import operator
import re


class NonOverlappingSet:
    """
    Stores non-overlapping intervals and allows for fast retrieval of overlapping intervals.
    """
    def __init__(self, starts, ends):
        assert len(starts) == len(ends)
        self._starts = np.array(starts)
        self._ends = np.array(ends)
        assert np.all(self._ends[:-1] <= self._starts[1:])

    @classmethod
    def from_start_end_pairs(cls, pairs):
        self = cls.__new__(cls)
        self._starts = np.fromiter(map(operator.itemgetter(0), pairs), np.int32, len(pairs))
        self._ends = np.fromiter(map(operator.itemgetter(1), pairs), np.int32, len(pairs))
        assert np.all(self._ends[:-1] <= self._starts[1:])
        return self

    def select(self, start, end):
        start_ix = self._ends.searchsorted(start, side='right')
        end_ix = self._starts.searchsorted(end, side='left')
        return start_ix, end_ix


def load_summaries(inputs):
    # Dictionary (sample, chrom): [(start, end, rest_of_line)].
    results = collections.defaultdict(list)
    unique_events = collections.defaultdict(dict)
    for line in itertools.chain(*inputs):
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        chrom, start, end = line[:3]
        start = int(start)
        end = int(end)
        sample = line[4]
        results[(sample, chrom)].append((start, end, line[3:]))
        if (start, end) not in unique_events[chrom]:
            unique_events[chrom][(start, end)] = [line[3:]]
        else:
            unique_events[chrom][(start, end)].append(line[3:])

    results2 = {}
    for (sample, chrom), curr_results in results.items():
        curr_results.sort()
        searcher = NonOverlappingSet.from_start_end_pairs(curr_results)
        results2[(sample, chrom)] = (curr_results, searcher)
    return results2, unique_events


def sort_keys(results):
    samples = set(map(operator.itemgetter(0), results.keys()))
    chroms = set(map(operator.itemgetter(1), results.keys()))
    samples = sorted(samples)

    convert = lambda text: int(text) if text.isdigit() else text.lower()
    natural_key = lambda key: tuple(convert(c) for c in re.split('([0-9]+)', key))
    chroms = sorted(chroms, key=natural_key)
    return samples, chroms


def process_sample_entries(chrom, start, end, samples, results):
    CNF = 0
    CN = 1
    CNQ = 2
    PCNF = 3
    PCN = 4
    PCNQ = 5

    all_copy_nums = collections.Counter()
    all_paralog_cns = collections.Counter()
    sample_results = ''

    for sample in samples:
        sample_results += '\t'
        sample_entries, searcher = results.get((sample, chrom), (None, None))
        if sample_entries is None:
            sample_results += '*'
            continue

        start_ix, end_ix = searcher.select(start, end)
        if end_ix <= start_ix:
            sample_results += '*'
            continue

        info_sets = [set() for _ in range(6)]
        coverage = 0
        # Convert to int because islice does not work with numpy.int.
        for entry_start, entry_end, entry in itertools.islice(sample_entries, int(start_ix), int(end_ix)):
            curr_info = entry[3:9]
            for i, value in enumerate(curr_info):
                info_sets[i].add(value)
            curr_cov = min(end, entry_end) - max(start, entry_start)

            all_copy_nums[curr_info[CN]] += curr_cov
            all_paralog_cns[curr_info[PCN]] += curr_cov
            coverage += curr_cov

        coverage = 100 * coverage / (end - start)
        if len(info_sets[1]) != 1:
            sample_results += '! ! ! ! ! ! {:.1f}'.format(coverage)
            continue
        # Replace with constants.
        for i in (CN, CNF, CNQ, PCN, PCNF, PCNQ):
            if len(info_sets[i]) == 1:
                sample_results += '{} '.format(info_sets[i].pop())
            else:
                sample_results += '! '
        sample_results += '{:.1f}'.format(coverage)
    return all_copy_nums, all_paralog_cns, sample_results


def write_header(out, samples):
    out.write('## For each sample 7 values are stored: copy_num, copy_num_filter, copy_num_qual; '
        'paralog_copy_num, paralog_filter, paralog_qual; coverage.\n')
    out.write('## coverage - percentage of the region covered by the sample entries.\n')
    out.write('## Entries for sample can contain "!", that means that several entries '
        'cover the region and have different values.\n')
    out.write('#chrom\tstart\tend\tregion\tref\tcopy_num_freqs\tparalog_freqs\tinfo\tother_regions\t')
    out.write('\t'.join(samples))
    out.write('\n')


def create_matrix(results, unique_events, chrom, samples, out):
    unique_events = unique_events[chrom]
    for start, end in sorted(unique_events.keys()):
        templates = unique_events[(start, end)]
        template = templates[0]
        region_name = template[0]
        pos2 = template[2]
        if pos2 == '' or pos2 == '*':
            ref_copy_num = 2
        else:
            assert pos2.startswith('chr')
            ref_copy_num = 2 * (pos2.count(',') + 2)

        out.write('{}\t{}\t{}\t{}\t{}\t'.format(chrom, start, end, region_name, ref_copy_num))
        all_copy_nums, all_paralog_cns, sample_results = process_sample_entries(chrom, start, end, samples, results)

        copy_num_freqs = []
        for copy_num, freq in all_copy_nums.items():
            if copy_num.isdigit():
                # Store sorting key, copy number, and frequency.
                copy_num_freqs.append((int(copy_num), copy_num, freq))
            elif copy_num.startswith('<'):
                copy_num_freqs.append((-1, copy_num, freq))
            elif copy_num.startswith('>'):
                copy_num_freqs.append((1000, copy_num, freq))
            # else: ignore
        copy_num_freqs.sort()
        copy_num_freq_sum = sum(map(operator.itemgetter(2), copy_num_freqs))
        out.write(' '.join('{}={:.5g}'.format(copy_num, freq / copy_num_freq_sum)
            for _, copy_num, freq in copy_num_freqs))
        out.write('\t')

        paralog_freq_sum = sum(all_paralog_cns.values())
        out.write(' '.join('{}={:.5g}'.format(paralog, freq / paralog_freq_sum)
            for paralog, freq in all_paralog_cns.most_common()))

        info = 'len={:.1f}kb;samples={}:{}{}'.format((end - start) / 1000, len(templates),
            ','.join(entry[1] for entry in itertools.islice(templates, 0, 10)),
            ',...' if len(templates) > 10 else '')
        out.write('\t{}\t'.format(info))

        out.write(pos2)

        out.write(sample_results)
        out.write('\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar='<file>', nargs='+', type=argparse.FileType(),
        help='Input summary files.')
    parser.add_argument('-o', '--output', metavar='<file>', type=argparse.FileType('w'),
        help='Output summary matrix.')
    args = parser.parse_args()

    results, unique_events = load_summaries(args.input)
    samples, chroms = sort_keys(results)
    write_header(args.output, samples)
    for chrom in chroms:
        create_matrix(results, unique_events, chrom, samples, args.output)


if __name__ == '__main__':
    main()
