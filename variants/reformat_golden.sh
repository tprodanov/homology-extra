#!/bin/bash

set -eu

dir="$1"
regions="$2/calling.bed"

a="$dir/golden.vcf.gz"
b="$dir/golden_tmp.vcf.gz"
c="$dir/calls/golden.vcf.gz"

zcat "$a" | sed 's/INFO$/INFO\tFORMAT\tsim/; s/WP=\(.*\)/.\tGT\t\1/;
    s/#CHROM/##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM/' | bgzip > "$b"
tabix -p vcf "$b"
bcftools view -T "$regions" "$b" | bgzip > "$c"
tabix -p vcf "$c"

rm "$b" "$b".tbi
