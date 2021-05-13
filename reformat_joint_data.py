#!/usr/bin/env python3

import sys
import gzip


def parse_line(line):
    line = line.strip().split('\t')
    info = {}
    if line[3] != '*':
        for entry in line[3].split(';'):
            key, value = entry.split('=')
            value = value.split(',')
            info[key] = value
    return line[:3], info


def reformat_line(line):
    (ty, reg1, reg2), info = parse_line(line)
    if ty == 'main':
        if len(info['reliable_threshold']) == 1:
            info['reliable_threshold'] = ['0.8', '0.95']
        if 'copy_num_jump' not in info:
            info['copy_num_jump'] = '6'
        if 'transition_prob' not in info:
            info['transition_prob'] = '-5'

    if ty == 'hmm_paths':
        info['initial'] = [val.lstrip('-') for val in info['initial']]

    if ty == 'window':
        if 'jumps' in info:
            n = len(info['jumps'])
            if n != 3:
                info['jumps'] = info['jumps'][n // 2 - 1 : n // 2 + 2]
            info['jumps'] = [val.lstrip('-') for val in info['jumps']]

    return '{}\t{}\t{}\t{}\n'.format(ty, reg1, reg2,
        ';'.join('{}={}'.format(key, ','.join(value)) for key, value in info.items()) if info else '*')


def main():
    in_filename = sys.argv[1]
    out_filename = sys.argv[2]
    with gzip.open(in_filename, 'rt') as inp, gzip.open(out_filename, 'wt') as out:
        for line in inp:
            out.write(reformat_line(line))


if __name__ == '__main__':
    main()
