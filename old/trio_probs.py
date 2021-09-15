#!/usr/bin/env python3

import argparse
import csv
import numpy as np
import scipy.optimize
from scipy.special import logsumexp
import collections
import sys


def load_input(input):
    res = collections.defaultdict(dict)
    for inp in input:
        reader = csv.DictReader(inp, fieldnames=('chr start end locus sample reg2 agCN_filter agCN agCN_qual '
            'psCN_filter psCN psCN_qual gene exon').split(), delimiter='\t')
        for row in reader:
            gene = '{}:{}'.format(row['gene'], row['exon'])
            res[gene][row['sample']] = row
    return res


def load_trios(input):
    reader = csv.DictReader(input, delimiter='\t')
    res = []
    for row in reader:
        father = row['Paternal ID']
        mother = row['Maternal ID']
        if father == '0' or mother == '0':
            continue
        sample = row['Individual ID']
        pop = row['Population']
        res.append((sample, father, mother, pop))
    return res


def get_obs(row, paralog, qual):
    if row['agCN_filter'] != 'PASS' or float(row['agCN_qual']) < qual or row['psCN_filter'] != 'PASS':
        return None
    try:
        if float(row['psCN_qual'].split(',')[paralog]) < qual:
            return None
        return int(row['psCN'].split(',')[paralog])
    except (ValueError, IndexError):
        return None


def get_fractions(x):
    if np.sum(x) > 1:
        return None
    fractions = np.hstack((x, (0.,)))
    fractions[-1] = 1 - np.sum(fractions)
    fractions /= np.sum(fractions)
    return fractions


def get_exp_counts(fractions, n_values):
    n = len(fractions)
    exp_counts = np.zeros(n * 2 - 1)
    for i in range(n):
        mult = n_values * fractions[i]
        for j in range(n):
            exp_counts[i + j] += mult * fractions[j]
    return exp_counts


def arrstr(arr):
    return '[{}]'.format('  '.join('{:10.1f}'.format(x) for x in arr))


def estimate_pop_frequencies(values, log):
    max_val = max(values)
    counts = np.bincount(values, minlength=max_val * 2 + 1).astype(np.float)
    n_values = len(values)

    def inner(x):
        if np.min(x) < 0 or np.max(x) > 1:
            return 1e8

        fractions = get_fractions(x)
        if fractions is None:
            return 1e8

        exp_counts = get_exp_counts(fractions, n_values)
        ixs = np.where((exp_counts > 0) | (counts > 0))[0]
        # log.write('\n\nx = {}    fractions = {}\n'.format(x, fractions))
        # log.write('    counts = {}\n'.format(arrstr(counts)))
        # log.write('exp_counts = {}\n'.format(arrstr(exp_counts)))
        # log.write('       ixs = {}\n'.format(arrstr(ixs)))
        # log.write('diff   = {}\n'.format(arrstr(counts[ixs] - exp_counts[ixs])))
        # log.write('diff^2 = {}\n'.format(arrstr(np.power(counts[ixs] - exp_counts[ixs], 2.0))))
        # log.write('denom  = {}\n'.format(arrstr(exp_counts[ixs])))
        # log.write('frac   = {}\n'.format(arrstr(np.power(counts[ixs] - exp_counts[ixs], 2.0) / exp_counts[ixs])))
        chi2 = np.sum(np.power(counts[ixs] - exp_counts[ixs], 2.0) / exp_counts[ixs])
        # log.write('chi2 = {:.6f}\n'.format(chi2))
        return chi2

    x0 = np.full(max_val, 0.1 / max_val)
    x0[1] = 0.9
    sol = scipy.optimize.minimize(inner, x0, method='Nelder-Mead')
    fractions = get_fractions(sol.x)
    log.write('    Chi2 value {:.2f} is achieved with fractions {}\n'.format(sol.fun, fractions))
    if sol.fun > 1000:
        log.write('    WARN: Chi2 is very high!\n')
    log.write('    Observed counts: {}\n'.format(arrstr(counts)))
    assert fractions is not None
    exp_counts = get_exp_counts(fractions, n_values)
    log.write('    Expected counts: {}\n'.format(arrstr(exp_counts)))
    return fractions


QUARTER = np.log(0.25)
LOG10 = np.log(10)


def calculate_prob(cn_i, cn_f, cn_m, fractions, verbose=False):
    if verbose:
        print(f'Father: {cn_f},  Mother: {cn_m},  Indiv: {cn_i}')
    f_probs = np.full(cn_f + 1, np.nan)
    for f1 in range(cn_f + 1):
        f2 = cn_f - f1
        f_probs[f1] = fractions[f1] + fractions[f2]
    f_probs -= logsumexp(f_probs)
    if verbose:
        print(f'    father probs = {np.exp(f_probs)}')

    m_probs = np.full(cn_m + 1, np.nan)
    for m1 in range(cn_m + 1):
        m2 = cn_m - m1
        m_probs[m1] = fractions[m1] + fractions[m2]
    m_probs -= logsumexp(m_probs)
    if verbose:
        print(f'    mother probs = {np.exp(m_probs)}')

    res = [-np.inf]
    res_probs = np.zeros(max(cn_i, cn_f + cn_m) + 1)
    for f1, f_prob in enumerate(f_probs):
        f2 = cn_f - f1
        for m1, m_prob in enumerate(m_probs):
            m2 = cn_m - m1
            prob = QUARTER + f_prob + m_prob

            for x in (f1, f2):
                for y in (m1, m2):
                    if verbose:
                        print(f'    {x + y} ({x} + {y}) = {np.exp(prob)}')
                    res_probs[x + y] += np.exp(prob)
                    if x + y == cn_i:
                        res.append(prob)
            if verbose:
                print()
    if verbose:
        print(f'    -> {res_probs}')
    return logsumexp(res) / LOG10


def analyze_trios(gene, trios, paralog, observations, fractions, args):
    qual = args.quality
    output = args.output
    gene, exon = gene.split(':')
    fractions = np.log(fractions)

    all_probs = []
    for indiv, father, mother, population in trios:
        indiv_row = observations.get(indiv)
        father_row = observations.get(father)
        mother_row = observations.get(mother)
        if indiv_row is None or father_row is None or mother_row is None:
            continue

        indiv_cn = get_obs(indiv_row, paralog, qual)
        father_cn = get_obs(father_row, paralog, qual)
        mother_cn = get_obs(mother_row, paralog, qual)
        if indiv_cn is None or father_cn is None or mother_cn is None:
            continue

        prob = calculate_prob(indiv_cn, father_cn, mother_cn, fractions)
        output.write(f'{gene}\t{exon}\t{paralog + 1}\t{indiv}\t{father}\t{mother}\t{population}\t'
            f'{indiv_cn}\t{father_cn}\t{mother_cn}\t{prob:.3g}\n')
        all_probs.append(prob)
    return np.array(all_probs)


def analyze_gene(gene, paralog, trios, observations, args):
    args.log.write(f'==============\nGene {gene} paralog {paralog + 1}\n')
    qual = args.quality

    pscn_values = []
    for row in observations.values():
        pscn = get_obs(row, paralog, qual)
        if pscn is not None:
            pscn_values.append(pscn)
    args.log.write(f'    Observations present in {len(pscn_values)} / {len(observations)} samples\n')
    if len(pscn_values) < 100:
        return

    pscn_values = np.array(pscn_values)
    fractions = estimate_pop_frequencies(pscn_values, args.log)
    all_probs = analyze_trios(gene, trios, paralog, observations, fractions, args)
    gene, exon = gene.split(':')

    n_copies = len(next(iter(observations.values()))['psCN'].split(','))
    args.summary.write(f'{gene}\t{exon}\t{paralog + 1}\t{n_copies * 2}\t'
        f'{len(pscn_values)}\t{np.mean(pscn_values):.3g}\t{np.std(pscn_values):.3g}\t'
        f'{len(all_probs)}\t{sum(all_probs >= -3)}\t{sum(all_probs >= -2)}\t{sum(all_probs >= -1)}\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', nargs='+', type=argparse.FileType(), required=True,
        help='Input filtered res.samples.bed.')
    parser.add_argument('-t', '--trios', type=argparse.FileType(), required=True,
        help='Input file with trio information.')
    parser.add_argument('-l', '--log', type=argparse.FileType('w'), required=True,
        help='Output log file.')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), required=True,
        help='Output csv file.')
    parser.add_argument('-s', '--summary', type=argparse.FileType('w'), required=True,
        help='Output summary file.')
    parser.add_argument('-q', '--quality', type=float, default=20,
        help='Minimal quality.')
    args = parser.parse_args()

    np.set_printoptions(precision=6, linewidth=sys.maxsize, suppress=True, threshold=sys.maxsize)
    obs = load_input(args.input)
    trios = load_trios(args.trios)
    args.output.write('# {}\n'.format(' '.join(sys.argv)))
    args.output.write('gene\texon\tparalog\tindiv\tfather\tmother\tpopulation\t'
        'indiv_psCN\tfather_psCN\tmother_psCN\tprob\n')
    args.summary.write('# {}\n'.format(' '.join(sys.argv)))
    args.summary.write('# ge_x = log10(probability) >= -x\n')
    args.summary.write('gene\texon\tparalog\tref_cn\tn_obs\tcn_mean\tcn_std\tn_trios\tge_3\tge_2\tge_1\n')

    for gene, observations in obs.items():
        if not observations:
            continue
        n_copies = len(next(iter(observations.values()))['psCN'].split(','))
        assert n_copies > 1
        for paralog in range(n_copies):
            analyze_gene(gene, paralog, trios, observations, args)


if __name__ == '__main__':
    main()
