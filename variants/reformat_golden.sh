#!/bin/bash

set -eu

if [[ $# -ge 1 ]]; then
    dir="$1"
else
    dir="."
fi

a="$dir/golden.vcf.gz"
b="$dir/golden_phased.vcf.gz"

zcat "$a" | sed 's/INFO$/INFO\tFORMAT\tsim/; s,WP=\([0-9]\+\)/\([0-9]\+\),.\tGT\t\1|\2,;
    s/#CHROM/##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM/' | bgzip > "$b"
tabix -p vcf "$b"
