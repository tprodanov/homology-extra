#!/usr/bin/env python3

import argparse
import sys
import csv
from glob import glob
import os
import operator
import numpy as np
from collections import defaultdict


def get_filenames(in_dir, verbose=False):
    filenames = glob('{}/*/**/summary.bed'.format(in_dir), recursive=True) \
        + glob('{}/*/**/res.samples.bed'.format(in_dir), recursive=True)
    by_locus = defaultdict(list)
    for filename in filenames:
        rel_filename = os.path.relpath(filename, in_dir)
        filename_split = rel_filename.split('/')
        locus = filename_split[0]
        key = (len(filename_split), filename_split[1])
        by_locus[locus].append((key, filename))

    res = []
    for locus, filenames in by_locus.items():
        filenames.sort()
        res.append((locus, filenames[-1][1]))
        if 'v100' not in res[-1][-1]:
            print(res[-1][-1])
            assert False
        if len(filenames) > 1 and verbose:
            sys.stderr.write('For locus {} select filename {} (out of {})\n'.format(locus, res[-1][1],
                ', '.join(map(operator.itemgetter(1), filenames))))
    res.sort()
    return res


def sum_dicts(dicts):
    cns = sorted(dicts.keys())
    new_dict = nested_dict(2, list)
    for cn in cns:
        for key, values in dicts[cn].items():
            for entry, subval in values.items():
                new_dict[key][entry].extend(subval)
    dicts['all'] = new_dict


def write_results(name, ty, dicts, writer, args, all_cn=True):
    """
    dicts: { ref_cn: { key: { entry: [ (length, flag) ] } } }
    """
    cns = sorted(dicts.keys())
    if all_cn and len(cns) > 1:
        sum_dicts(dicts)
        cns.append('all')

    for cn in cns:
        curr_dicts = dicts[cn]
        entries = set()
        for subdict in curr_dicts.values():
            entries.update(subdict.keys())
        entries = sorted(entries)
        n_entries = len(entries)

        outrow = {}
        outrow['name'] = name
        outrow['type'] = ty
        outrow['entries'] = str(n_entries)
        outrow['refCN'] = cn if isinstance(cn, str) else str(cn * 2)

        length = 0
        for d in curr_dicts.values():
            try:
                length = max(length, max(sum(map(operator.itemgetter(0), values)) for values in d.values()))
            except ValueError:
                pass
        outrow['length'] = length

        for key, values in curr_dicts.items():
            diff_length = ':diff_len:' in key
            only_mean = ':only_mean:' in key
            key = key.split(':', 1)[0]

            cov = np.zeros(n_entries)
            for i, entry in enumerate(entries):
                total_len = sum(map(operator.itemgetter(0), values[entry])) if diff_length else length
                if total_len > 0:
                    cov[i] = sum(sublen * flag / total_len * 100.0 for (sublen, flag) in values[entry])
                else:
                    cov[i] = np.nan

            cov = cov[~np.isnan(cov)]
            if len(cov) == 0:
                continue
            if only_mean:
                outrow[key] = '{:.2f}'.format(np.mean(cov))
                continue

            outrow['{}_mean'.format(key)] = '{:.2f}'.format(np.mean(cov))
            for threshold in args.thresholds:
                outrow['{}_over{:.5g}'.format(key, threshold)] = np.sum(cov >= threshold)
            for i, val in enumerate(np.quantile(cov, np.array(args.percentiles) / 100)):
                outrow['{}_perc{}'.format(key, args.percentiles[i])] = '{:.2f}'.format(val)
        writer.writerow(outrow)


def dget(row, *keys):
    for key in keys:
        value = row.get(key)
        if value is not None:
            return value
    raise RuntimeError('There are no keys {} in the row {}'.format(', '.join(keys), row))


def analyze_locus(locus, filename, writer, args, sample_dicts):
    with open(filename) as inp:
        fieldnames = None
        for line in inp:
            if line.startswith('##'):
                continue
            fieldnames = line.lstrip('#').strip().split('\t')
            break

        assert fieldnames is not None
        reader = csv.DictReader(inp, fieldnames, delimiter='\t')

        dicts = nested_dict(3, list)
        for row in reader:
            reg2 = dget(row, 'other_regions', 'homologous_regions')
            ref_cn = reg2.count(',') + bool(reg2 != '*') + 1

            sample = row['sample']
            length = int(row['end']) - int(row['start'])
            ag_cn_filt = dget(row, 'copy_num_filter', 'agCN_filter')
            ag_cn_qual = float(dget(row, 'copy_num_qual', 'agCN_qual'))
            ag_cn_high_qual = ag_cn_filt == 'PASS' and ag_cn_qual >= args.qual
            dicts[ref_cn]['agCN'][sample].append((length, ag_cn_high_qual))
            if ag_cn_high_qual:
                obs_ag_cn = int(dget(row, 'copy_num', 'agCN').lstrip('<>'))
                dicts[ref_cn]['obs_agCN:diff_len:only_mean:'][sample].append((length, obs_ag_cn / 100))
            sample_dicts[sample][ref_cn]['agCN:diff_len:'][locus].append((length, ag_cn_high_qual))

            par_filter = dget(row, 'paralog_filter', 'psCN_filter')
            par_cn = dget(row, 'paralog_copy_num', 'psCN')
            par_qual = dget(row, 'paralog_qual', 'psCN_qual')
            par_qual = np.array([float(x) if x != '*' else 0.0 for x in par_qual.split(',')])
            ps_cn_high_qual = ag_cn_high_qual \
                and par_filter == 'PASS' and par_cn != '*' and '?' not in par_cn and np.all(par_qual >= args.qual)
            dicts[ref_cn]['psCN'][sample].append((length, ps_cn_high_qual))
            sample_dicts[sample][ref_cn]['psCN:diff_len:'][locus].append((length, ps_cn_high_qual))
        write_results(locus, 'gene', dicts, writer, args)


def create_writer(out, keys, args):
    out.write('# {}\n'.format(' '.join(sys.argv)))
    out.write('# Qual = {:.1f}\n'.format(args.qual))
    fieldnames = 'name type refCN length entries'.split()
    for key in keys:
        only_mean = ':only_mean:' in key
        key = key.split(':', 1)[0]
        if only_mean:
            fieldnames.append(key)
            continue

        fieldnames.append('{}_mean'.format(key))
        for p in args.percentiles:
            fieldnames.append('{}_perc{}'.format(key, p))
        for t in args.thresholds:
            fieldnames.append('{}_over{:.5g}'.format(key, t))
    writer = csv.DictWriter(out, fieldnames, delimiter='\t', dialect='unix', quoting=csv.QUOTE_MINIMAL, restval='NA')
    writer.writeheader()
    return writer


def nested_dict(depth, last_fn):
    if depth <= 1:
        return defaultdict(last_fn)

    def inner():
        return nested_dict(depth - 1, last_fn)
    return defaultdict(inner)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar='<dir>', required=True,
        help='Input directory.')
    parser.add_argument('-o', '--output', metavar='<file>', required=True, type=argparse.FileType('w'),
        help='Output csv file.')
    parser.add_argument('-q', '--qual', type=float, metavar='<float>', default=20.0,
        help='Quality threshold [%(default)s].')
    parser.add_argument('-p', '--percentiles', type=int, metavar='<int>', nargs='*', default=(0, 10, 25, 50, 75, 100),
        help='Output percentiles [%(default)s].')
    parser.add_argument('-t', '--thresholds', type=float, metavar='<float>', nargs='*', default=(99,),
        help='Output number of samples that overcome a threshold.')
    args = parser.parse_args()

    writer = create_writer(args.output, ('obs_agCN:only_mean:', 'agCN', 'psCN'), args)
    filenames = get_filenames(args.input)
    # filenames = [(locus, filename) for locus, filename in filenames if locus == 'ABCC6']

    sample_dicts = nested_dict(4, list)
    for locus, filename in filenames:
        analyze_locus(locus, filename, writer, args, sample_dicts)

    for sample in sorted(sample_dicts.keys()):
        write_results(sample, 'sample', sample_dicts[sample], writer, args)


if __name__ == '__main__':
    main()
