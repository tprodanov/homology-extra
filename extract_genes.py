#!/usr/bin/env python3

import sys
import os
import argparse


def parse_info(info):
    return dict(entry.split('=') for entry in info.split(';'))


def get_gene_name(info):
    if 'Name' in info:
        return info['Name']
    return None


def load_genes(gff_inp, gene_names):
    for line in gff_inp:
        line = line.strip().split('\t')
        if line[2] != 'gene':
            continue
        info = parse_info(line[8])
        gene_name = get_gene_name(info)
        if gene_name in gene_names:
            yield (line[0], line[3], line[4], gene_name)


def main():
    parser = argparse.ArgumentParser(
        description='Extract genes from GFF file',
        formatter_class=argparse.RawTextHelpFormatter, add_help=False,
        usage='%(prog)s -i <gff> (-g <txt> | -G <gene> [<gene> ...]) -o <bed> [arguments]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', metavar='<file>', required=True,
        help='Input GFF annotation file.')
    genes_ma = io_args.add_mutually_exclusive_group(required=True)
    genes_ma.add_argument('-g', '--gene-file', metavar='<file>',
        help='Input file with gene names.')
    genes_ma.add_argument('-G', '--genes', metavar='<gene>', nargs='+',
        help='At least one gene name.')
    io_args.add_argument('-o', '--output', metavar='<file>',
        help='Output bed file with genes.')

    oth_args = parser.add_argument_group('Other arguments')
    oth_args.add_argument('-h', '--help', action='help', help='Show this message and exit.')
    args = parser.parse_args()

    if args.gene_file:
        with open(args.gene_file) as inp:
            gene_names = set(map(str.strip, inp))
    else:
        gene_names = set(args.genes)
    not_found = set(gene_names)

    with open(args.input) as inp, open(args.output, 'w') as outp:
        outp.write('# %s\n' % ' '.join(sys.argv))
        for gene in load_genes(inp, gene_names):
            not_found.remove(gene[3])
            outp.write('\t'.join(gene))
            outp.write('\n')

    if not_found:
        print('Could not find %d genes' % len(not_found))
    for name in not_found:
        print(name)


if __name__ == '__main__':
    main()
