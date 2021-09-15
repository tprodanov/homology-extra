#!/usr/bin/env python3

import sys
import datetime
import collections


def parse_time(s):
    hours, minutes, seconds = map(float, s.split(':'))
    return datetime.timedelta(seconds=hours * 3600 + minutes * 60 + seconds)


def parse_entries(entries, res):
    prev_time = datetime.timedelta()
    gene = None

    for entry in entries:
        s_time, msg = entry.strip('\n').split(maxsplit=1)
        time = parse_time(s_time)
        if msg.startswith('Analyzing'):
            gene = msg.split()[-1][1:-1]
            continue
        job_time = time - prev_time
        prev_time = time
        res[gene][msg] = job_time


def summarize(res):
    total = datetime.timedelta()
    for gene, gene_times in res.items():
        gene_time = sum(gene_times.values(), datetime.timedelta())
        print(f'{gene}\t{gene_time}')
        total += gene_time
    print(f'Total\t{total}\t{total.total_seconds()} sec')


def main():
    res = collections.defaultdict(dict)
    for filename in sys.argv[1:]:
        sys.stderr.write(f'Analyzing {filename}\n')

        with open(filename) as inp:
            lines = []
            for line in inp:
                if line.startswith('====') and lines:
                    parse_entries(lines, res)
                    lines.clear()
                else:
                    lines.append(line)
    summarize(res)


if __name__ == '__main__':
    main()
