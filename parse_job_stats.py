#!/usr/bin/env python3

import sys
import operator


def transform_mem(memory):
    assert memory.endswith('kb')
    memory = int(memory[:-2])
    return '{:10,} kb = {:5.2f} Gb'.format(memory, memory / 1024 / 1024)


def parse_resources(text):
    PREFIX = 'resources_used.'
    d = {}
    for entry in text.strip().split():
        if not entry.startswith(PREFIX):
            continue
        key, value = entry[len(PREFIX) : ].split('=')
        d[key] = value

    keys = [('walltime', None), ('cput', 'cputime'), ('mem', None), ('vmem', None)]
    for key, new_name in keys:
        if key in d:
            value = d[key]
            if value.endswith('kb'):
                value = transform_mem(value)
            print('{:8} = {}'.format(new_name or key, value))

    keys = set(map(operator.itemgetter(0), keys))
    for key, value in d.items():
        if key not in keys:
            print('{:8} = {}'.format(key, value))


def transform_time(time):
    assert time.endswith('s')
    minutes, seconds = time[:-1].split('m')
    minutes = int(minutes)
    seconds = float(seconds)
    return '{:>3}:{:02}:{:02.0f} ({:>4.0f} minutes)'.format(minutes // 60, minutes % 60, seconds, minutes + seconds / 60)


def parse_time(text):
    lines = list(filter(bool, map(str.strip, text.split('\n'))))
    for line in lines:
        key, value = line.split()
        if 'm' in value and value.endswith('s'):
            value = transform_time(value)
        print('{:4} = {}'.format(key, value))


def main():
    text = sys.stdin.read()
    if 'resources' in text:
        parse_resources(text)
    else:
        parse_time(text)


if __name__ == '__main__':
    main()
