#!/usr/bin/env python3

import argparse
import csv
import numpy as np
import scipy.optimize
import scipy.stats
from scipy.special import logsumexp
from collections import defaultdict, namedtuple
import sys
import os
import operator
import warnings


def load_input(input):
    res = defaultdict(dict)
    for inp in input:
        reader = csv.DictReader(inp, fieldnames=(
            'chr start end locus sample agCN_filter agCN agCN_qual '
            'psCN_filter psCN psCN_qual info reg2 _chr _start _end gene').split(), delimiter='\t')
        for row in reader:
            res[row['gene']][row['sample']] = row
    return res


def load_populations(input):
    reader = csv.DictReader(input, delimiter='\t')
    res = {}
    for row in reader:
        pop_code = row['Population code'].strip()
        spop_code = row['Superpopulation code'].strip()
        if pop_code and spop_code:
            res[pop_code] = spop_code
    return res


Trio = namedtuple('Trio', 'sample father mother pop spop')


def load_trios(input, populations, super_population):
    reader = csv.DictReader(input, delimiter='\t')
    trios = []
    sample_pops = {}

    for row in reader:
        sample = row['Individual ID']
        pop = row['Population']
        sample_pops[sample] = spop = populations.get(pop, '*')
        if super_population is not None and spop != super_population:
            continue

        father = row['Paternal ID']
        mother = row['Maternal ID']
        if father != '0' and mother != '0':
            trios.append(Trio(sample, father, mother, pop, spop))

    for trio in trios:
        if sample_pops[trio.father] != trio.spop or sample_pops[trio.mother] != trio.spop:
            print(f'WARN: trio {trio} has conflicting populations')
    return trios, sample_pops


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
    assert np.min(x) >= 0
    fractions = x / np.sum(x)
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


def define_x0s(n):
    for z in [0.5, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999]:
        w = (1 - z) / (n - 1)
        x0 = np.full(n, w)
        x0[1] = z
        yield x0


def estimate_pop_frequencies(values, log):
    max_val = max(values)
    counts = np.bincount(values, minlength=max_val * 2 + 1).astype(np.float)
    n_values = len(values)

    best_fractions = None
    best_chi2 = [np.inf]

    def inner(x):
        assert np.min(x) >= 0
        fractions = get_fractions(x)
        exp_counts = get_exp_counts(fractions, n_values)
        ixs = np.where((exp_counts > 0) | (counts > 0))[0]
        if not np.all(np.isfinite(exp_counts[ixs])):
            return 1e8
        chi2 = np.sum(np.power(counts[ixs] - exp_counts[ixs], 2.0) / exp_counts[ixs])

        nonlocal best_chi2
        nonlocal best_fractions
        if chi2 < best_chi2:
            best_chi2 = chi2
            best_fractions = fractions
        return chi2

    for x0 in define_x0s(max_val + 1):
        tmp_sol = scipy.optimize.minimize(inner, x0, method='L-BFGS-B', bounds=((1e-8, 2),) * (max_val + 1),
            options=dict(eps=1e-5, maxiter=30000, maxfun=30000))
        log.write('       x = {:30s} -> chi2 {:.2f}\n'.format(str(get_fractions(tmp_sol.x)), tmp_sol.fun))

    log.write('    Chi2 value {:.2f} is achieved with fractions {}\n'.format(best_chi2, best_fractions))
    if best_chi2 > 20:
        log.write('    WARN: Chi2 is too high!\n')
    log.write('    Observed counts: {}\n'.format(arrstr(counts)))
    assert best_fractions is not None
    exp_counts = get_exp_counts(best_fractions, n_values)
    log.write('    Expected counts: {}\n'.format(arrstr(exp_counts)))
    return best_chi2, best_fractions


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

    res_probs = [[] for _ in range(max(cn_i, cn_f + cn_m) + 1)]
    for f1, f_prob in enumerate(f_probs):
        f2 = cn_f - f1
        for m1, m_prob in enumerate(m_probs):
            m2 = cn_m - m1
            prob = QUARTER + f_prob + m_prob

            for x in (f1, f2):
                for y in (m1, m2):
                    if verbose:
                        print(f'    {x + y} ({x} + {y}) = {np.exp(prob)}')
                    res_probs[x + y].append(prob)
            if verbose:
                print()
    if verbose:
        print(f'    -> {res_probs}')
    # max_val = np.max(res_probs)

    res_probs = [logsumexp(curr) if curr else -np.inf for curr in res_probs]
    prob = res_probs[cn_i] / LOG10
    norm_prob = (res_probs[cn_i] - np.max(res_probs)) / LOG10 if np.isfinite(prob) else prob
    return prob, norm_prob


def analyze_trios(gene, trios, paralog, observations, fractions, out_files, min_qual):
    output = out_files['output']
    fractions = np.log(fractions)

    all_probs = []
    norm_probs = []
    cn_matrix = []

    for trio in trios:
        indiv_row = observations.get(trio.sample)
        father_row = observations.get(trio.father)
        mother_row = observations.get(trio.mother)
        if indiv_row is None or father_row is None or mother_row is None:
            continue

        indiv_cn = get_obs(indiv_row, paralog, min_qual)
        father_cn = get_obs(father_row, paralog, min_qual)
        mother_cn = get_obs(mother_row, paralog, min_qual)
        if indiv_cn is None or father_cn is None or mother_cn is None:
            continue

        cn_matrix.append((indiv_cn, father_cn, mother_cn))
        prob, norm_prob = calculate_prob(indiv_cn, father_cn, mother_cn, fractions)
        output.write(f'{gene}\t{paralog + 1}\t{trio.sample}\t{trio.father}\t{trio.mother}\t{trio.pop}\t{trio.spop}\t'
            f'{indiv_cn}\t{father_cn}\t{mother_cn}\t'
            f'{prob:.3g}\t{10 ** prob:.3g}\t{norm_prob:.3g}\t{10 ** norm_prob:.3g}\n')
        all_probs.append(prob)
        norm_probs.append(norm_prob)
    return np.array(all_probs), np.array(norm_probs), np.array(cn_matrix)


THRESHOLDS = [0.01, 0.05, 0.1, 0.2, 0.5]


def analyze_gene(gene, paralog, n_copies, trios, sample_pops, observations, out_files, min_qual):
    log = out_files['log']
    summary = out_files['summary']

    for superpop in sorted(set(map(operator.attrgetter('spop'), trios))):
        log.write(f'==============\nGene {gene} paralog {paralog + 1} superpopulation {superpop}\n')
        curr_obs = { sample: row for sample, row in observations.items()
            if sample_pops.get(sample, '*') == superpop }
        curr_trios = [trio for trio in trios if trio.spop == superpop]

        pscn_values = []
        for row in curr_obs.values():
            pscn = get_obs(row, paralog, min_qual)
            if pscn is not None:
                pscn_values.append(pscn)

        log.write(f'    Observations present in {len(pscn_values)} / {len(curr_obs)} samples\n')
        if len(pscn_values) < 20:
            summary.write(f'{gene}\t{paralog + 1}\t{superpop}\t{n_copies * 2}\t' +
                f'{len(pscn_values)}\t{np.mean(pscn_values):.3g}\t{np.std(pscn_values):.3g}'.replace('\tnan', '\tNA') +
                ('\tNA' * (len(THRESHOLDS) * 2 + 4)) + '\n')
            return

        pscn_values = np.array(pscn_values)
        chi2, fractions = estimate_pop_frequencies(pscn_values, log)
        chi2_prob = scipy.stats.chi2.sf(chi2, df=len(fractions) - 1)
        all_probs, _norm_probs, cn_matrix = analyze_trios(gene, curr_trios, paralog, curr_obs, fractions, out_files,
            min_qual)
        if len(all_probs) > 0 and cn_matrix.shape != (len(all_probs), 3):
            raise ValueError('Error in locus {}, paralog {}: N probs = {},  CN matrix shape = {}'
                .format(gene, paralog + 1, len(all_probs), cn_matrix.shape))

        summary.write(f'{gene}\t{paralog + 1}\t{superpop}\t{n_copies * 2}\t'
            f'{len(pscn_values)}\t{np.mean(pscn_values):.3g}\t{np.std(pscn_values):.3g}\t'
            f'{chi2:.3f}\t{chi2_prob:.3g}')
        for i in range(2):
            curr_probs = all_probs if i == 0 or len(all_probs) == 0 else all_probs[cn_matrix[:, 0] != 2]
            summary.write(f'\t{len(curr_probs)}')
            if len(curr_probs) == 0:
                summary.write('\tNA' * len(THRESHOLDS))
            else:
                for thresh in THRESHOLDS:
                    log_thresh = np.log10(thresh)
                    summary.write('\t{}'.format(np.sum(curr_probs >= log_thresh)))
        summary.write('\n')
        summary.flush()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', nargs='+', type=argparse.FileType(), required=True,
        help='Input filtered res.samples.bed.')
    parser.add_argument('-t', '--trios', type=argparse.FileType(), required=True,
        help='Input file with trio information.')
    parser.add_argument('-p', '--populations', type=argparse.FileType(), required=False,
        help='Optional: populations.')
    parser.add_argument('-s', '--super-population', required=False,
        help='Optional: Limit the analysis to a single super-population.')
    parser.add_argument('-o', '--output', required=True,
        help='Output directory.')
    parser.add_argument('-q', '--quality', type=float, default=20,
        help='Minimal quality.')
    args = parser.parse_args()

    warnings.filterwarnings('ignore')

    try:
        os.mkdir(args.output)
    except FileExistsError:
        pass
    out_files = {
        'output': open(os.path.join(args.output, 'out.csv'), 'w'),
        'summary': open(os.path.join(args.output, 'summary.csv'), 'w'),
        'log': open(os.path.join(args.output, 'log.txt'), 'w'),
    }

    np.set_printoptions(precision=8, linewidth=sys.maxsize, suppress=True, threshold=sys.maxsize)
    obs = load_input(args.input)

    if args.populations:
        populations = load_populations(args.populations)
    else:
        populations = {}
    trios, sample_pops = load_trios(args.trios, populations, args.super_population)

    out_files['output'].write('# {}\n'.format(' '.join(sys.argv)))
    out_files['output'].write('gene\tparalog\tindiv\tfather\tmother\tpop\tsuperpop\t'
        'indiv_psCN\tfather_psCN\tmother_psCN\tlog10_prob\tprob\tlog10_norm_prob\tnorm_prob\n')

    summary_out = out_files['summary']
    summary_out.write('# {}\n'.format(' '.join(sys.argv)))
    summary_out.write('# ge_x: probability >= x\n')
    summary_out.write('gene\tparalog\tsuperpop\tref_cn\tn_obs\tcn_mean\tcn_std\tchi2\tchi2_pval\tn_trios')
    for thresh in THRESHOLDS:
        summary_out.write('\tprobs_ge_{0:.2g}'.format(thresh))
    summary_out.write('\tnonref_n_trios')
    for thresh in THRESHOLDS:
        summary_out.write('\tnonref_probs_ge_{0:.2g}'.format(thresh))
    summary_out.write('\n')
    summary_out.flush()

    for gene in sorted(obs):
        observations = obs[gene]
        assert observations
        reg2 = next(iter(observations.values()))['reg2']
        assert reg2 != '*'
        n_copies = 1 + len(reg2.split(','))
        for paralog in range(n_copies):
            analyze_gene(gene, paralog, n_copies, trios, sample_pops, observations, out_files, args.quality)


if __name__ == '__main__':
    main()
