#!/bin/bash

set -eu

if [[ $# -ge 1 ]]; then
    golden="$1"
else
    golden="golden.vcf.gz"
fi

out_golden="$(echo ${golden} | sed s/\.vcf\.gz/_phased.vcf.gz/)"
echo "Reformating golden VCF file: ${golden} -> ${out_golden}"

zcat "$golden" | grep -v "ID=WP," | \
    sed 's/INFO$/INFO\tFORMAT\tsim/; s,WP=\([0-9]\+\)/\([0-9]\+\),.\tGT\t\1|\2,;
        s/#CHROM/##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM/' | bgzip > "$out_golden"
tabix -p vcf "$out_golden"
