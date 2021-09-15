#!/usr/bin/env python3

import pysam
import sys


def rev_comp(seq):
    nts = dict(['AT', 'CG', 'GC', 'TA'])
    return ''.join(nts[nt] for nt in reversed(seq))


def main():
    genome = pysam.FastaFile(sys.argv[1])
    vcf = pysam.VariantFile(sys.argv[2])
    for record in vcf:
        seq1 = genome.fetch(record.chrom, record.start, record.start + len(record.ref))
        if seq1 != record.ref:
            print('Problem at the record:\n    %s' % str(record).strip())
            print('    Reference do not match (VCF = %s, genome = %s)' % (record.ref, seq1))
            continue

        for pos2 in record.info['pos2']:
            pos2_split = pos2.split(':')
            chrom = pos2_split[0]
            start = int(pos2_split[1]) - 1
            strand = pos2_split[2] == '+'
            if len(pos2_split) > 3:
                allele_str = record.alleles[int(pos2_split[3])]
            else:
                allele_str = record.alleles[1]

            end = start + len(allele_str)
            seq2 = genome.fetch(chrom, start, end)
            if not strand:
                seq2 = rev_comp(seq2)
            if seq2 != allele_str:
                print('Problem at the record:\n    %s' % str(record).strip())
                print('    Allele of the pos2 %s' % pos2)
                print('    Allele do not match (VCF = %s, genome = %s)' % (allele_str, seq2))
                continue


if __name__ == '__main__':
    main()
