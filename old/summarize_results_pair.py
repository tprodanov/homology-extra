#!/usr/bin/env python3

import argparse
import sys
import csv
from glob import glob
import os
import operator
import numpy as np
from collections import defaultdict

import summarize_results


def analyze_locus(locus, rows, writer, args, sample_dicts):
    qual = args.qual
    dicts = summarize_results.nested_dict(3, list)

    for row in rows:
        ref_cn = int(row['ref_cn']) // 2
        sample = row['sample']
        length = int(row['length'])

        ag_cn_present = row['copy_num_filter1'] == 'PASS' and row['copy_num_filter2'] == 'PASS' \
            and float(row['copy_num_qual1']) >= qual and float(row['copy_num_qual2']) >= qual
        dicts[ref_cn]['agCN_present:only_mean:'][sample].append((length, ag_cn_present))
        sample_dicts[sample][ref_cn]['agCN_present:only_mean:diff_len:'][locus].append((length, ag_cn_present))
        if ag_cn_present:
            obs_ag_cn = 0.5 * (int(row['copy_num1'].lstrip('<>')) + int(row['copy_num2'].lstrip('<>')))
            dicts[ref_cn]['obs_agCN:diff_len:only_mean:'][sample].append((length, obs_ag_cn / 100))

            ag_cn_match = row['copy_num1'] == row['copy_num2']
            dicts[ref_cn]['agCN:diff_len:'][sample].append((length, ag_cn_match))
            sample_dicts[sample][ref_cn]['agCN:diff_len:'][locus].append((length, ag_cn_match))
        else:
            ag_cn_match = False

        par_qual1 = row['paralog_qual1']
        par_qual2 = row['paralog_qual2']
        par_cn1 = row['paralog_copy_num1']
        par_cn2 = row['paralog_copy_num2']
        ps_cn_present = ag_cn_present and ag_cn_match \
            and row['paralog_filter1'] == 'PASS' and row['paralog_filter2'] == 'PASS' \
            and all(float(x) >= qual if x != '*' else False for x in row['paralog_qual1'].split(',')) \
            and all(float(x) >= qual if x != '*' else False for x in row['paralog_qual2'].split(',')) \
            and '?' not in par_cn1 and '?' not in par_cn2

        dicts[ref_cn]['psCN_present:only_mean:'][sample].append((length, ps_cn_present))
        sample_dicts[sample][ref_cn]['psCN_present:only_mean:diff_len:'][locus].append((length, ps_cn_present))
        if ps_cn_present:
            par_cn1 = par_cn1.split(',')
            par_cn2 = par_cn2.split(',')
            assert len(par_cn1) == len(par_cn2) == ref_cn
            ps_cn_match = all(x == y for x, y in zip(par_cn1, par_cn2))
            dicts[ref_cn]['psCN:diff_len:'][sample].append((length, ps_cn_match))
            sample_dicts[sample][ref_cn]['psCN:diff_len:'][locus].append((length, ps_cn_match))

    summarize_results.write_results(locus, 'gene', dicts, writer, args)


def load(inp):
    dataset1 = None
    dataset2 = None

    fieldnames = None
    for line in inp:
        if line.startswith('#'):
            continue
        fieldnames = line.strip().split('\t')
        break
    reader = csv.DictReader(inp, fieldnames, delimiter='\t')

    res = defaultdict(list)
    for row in reader:
        if dataset1 is None:
            dataset1 = row['dataset1']
            dataset2 = row['dataset2']
        else:
            assert dataset1 == row['dataset1']
            assert dataset2 == row['dataset2']
        res[row['locus']].append(row)
    return res


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar='<file>', required=True, type=argparse.FileType(),
        help='Input csv file with merged results.')
    parser.add_argument('-o', '--output', metavar='<file>', required=True, type=argparse.FileType('w'),
        help='Output csv file.')
    parser.add_argument('-q', '--qual', type=float, metavar='<float>', default=20.0,
        help='Quality threshold [%(default)s].')
    parser.add_argument('-p', '--percentiles', type=int, metavar='<int>', nargs='*', default=(0, 10, 25, 50, 75, 100),
        help='Output percentiles [%(default)s].')
    parser.add_argument('-t', '--thresholds', type=float, metavar='<float>', nargs='*', default=(99,),
        help='Output number of samples that overcome a threshold.')
    args = parser.parse_args()

    keys = ('obs_agCN:only_mean:', 'agCN_present:only_mean:', 'agCN', 'psCN_present:only_mean:', 'psCN')
    writer = summarize_results.create_writer(args.output, keys, args)
    merged_inp = load(args.input)

    sample_dicts = summarize_results.nested_dict(4, list)
    for locus, rows in merged_inp.items():
        analyze_locus(locus, rows, writer, args, sample_dicts)

    for sample in sorted(sample_dicts.keys()):
        summarize_results.write_results(sample, 'sample', sample_dicts[sample], writer, args)


if __name__ == '__main__':
    main()
