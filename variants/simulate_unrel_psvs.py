#!/usr/bin/env python3

import argparse
import pysam
import numpy as np
import random
import os
import sys
import subprocess


_rev_comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C','a':'T', 't':'A', 'c':'G', 'g':'C', 'N':'N', 'n':'N' }
def rev_comp(seq): # reverse complement of string
    return ''.join(_rev_comp.get(nt, 'X') for nt in seq[::-1])


def cond_rev_comp(seq, *, strand):
    if strand:
        return seq
    return rev_comp(seq)


def create_vcf_header(genome):
    vcf_header = pysam.VariantHeader()
    vcf_header.add_line('##command="%s"' % ' '.join(sys.argv))

    for name, length in zip(genome.references, genome.lengths):
        vcf_header.add_line('##contig=<ID={},length={}>'.format(name, length))

    vcf_header.add_line('##INFO=<ID=psv,Number=1,Type=Integer,Description="Initial PSV position.')
    vcf_header.add_line('##INFO=<ID=fval,Number=1,Type=Float,Description="f-value for this copy.')
    vcf_header.add_line('##INFO=<ID=rel,Number=1,Type=Character,Description="PSV reliability (r|s|u).')
    vcf_header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    vcf_header.add_sample('sim')
    return vcf_header


def process_psv(psv, rate, header):
    all_pos = (f'{psv.chrom}:{psv.pos}:+:0',) + psv.info['pos2']
    if 'fval' not in psv.info:
        return
    fvalues = psv.info['fval']
    if np.any(np.isnan(fvalues)):
        return
    status = psv.info['rel']

    for pos_i, pos in enumerate(all_pos):
        pos = pos.split(':')
        chrom = pos[0]
        start = int(pos[1]) - 1
        strand = pos[2] == '+'
        allele_ix = 1 if len(pos) == 3 else int(pos[3])

        allele_ixs = [allele_ix] + [i for i in range(len(psv.alleles)) if i != allele_ix]
        record = header.new_record()
        record.chrom = chrom
        record.start = start
        record.alleles = [cond_rev_comp(psv.alleles[i], strand=strand) for i in allele_ixs]
        record.info['psv'] = psv.pos

        record.info['fval'] = curr_fval = fvalues[pos_i]
        record.info['rel'] = status

        q = np.clip((1 - curr_fval) * rate, 0, 1)
        p = 1 - q
        ref_hom = p * p
        alt_hom = q * q
        heteroz = 1 - ref_hom - alt_hom
        r = random.random()
        alt_allele = random.randrange(1, len(psv.alleles))
        if r <= ref_hom:
            gt = (0, 0)
        elif r <= ref_hom + alt_hom:
            gt = (alt_allele, alt_allele)
        else:
            gt = (alt_allele, 0) if random.randrange(2) else (0, alt_allele)
        record.samples[0]['GT'] = gt
        yield record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True,
        help='Input Parascopy working directory.')
    parser.add_argument('-f', '--fasta-ref', required=True,
        help='Input FASTA reference file.')
    parser.add_argument('-r', '--rate', type=float, default=1,
        help='Unreliable PSV rate. PSV with f-value = x will be unreliable with probability (1-x)*rate.')
    parser.add_argument('-o', '--output', required=True,
        help='Output VCF file.')
    args = parser.parse_args()

    psvs = pysam.VariantFile(os.path.join(args.input, 'psvs.vcf.gz'))
    genome = pysam.FastaFile(args.fasta_ref)
    header = create_vcf_header(genome)

    records = []
    for psv in psvs:
        records.extend(process_psv(psv, args.rate, header))

    chrom_dict = { chrom: i for i, chrom in enumerate(genome.references) }
    records.sort(key=lambda record: (chrom_dict[record.chrom], record.start))
    with pysam.VariantFile(args.output, 'w', header=header) as out:
        for record in records:
            out.write(record)


if __name__ == '__main__':
    main()
