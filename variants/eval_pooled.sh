#!/bin/bash

set -eu

wdir="$(dirname "$0")"
genome_sdf="$1"

if [[ $# -ge 2 ]]; then
    dir="$2"
else
    dir="."
fi

if [[ $# -ge 3 ]]; then
    subdir="${dir}/$3"
else
    subdir="${dir}/calls"
fi

var_dir="$(dirname $(ls ${dir}/parascopy/*/variants_pooled.vcf.gz | head -1))"

zcat ${var_dir}/variants_pooled.vcf.gz | python <( echo '
import sys
for line in map(str.strip, sys.stdin):
    if "=GQ," in line:
        print(line)
        print("##FORMAT=<ID=GQ2,Number=1,Type=Float,Description=\"Genotype Quality OR 0 if GTfilter\">")
    elif line.startswith("#"):
        print(line)
    else:
        line = line.split("\t")
        line[8] += ":GQ2"
        if "GTfilter" in line[8]:
            line[9] += ":0"
        else:
            i = line[8].split(":").index("GQ")
            line[9] += ":" + line[9].split(":")[i]
        print("\t".join(line))
' ) | bgzip > ${var_dir}/variants_pooled2.vcf.gz
tabix -p vcf ${var_dir}/variants_pooled2.vcf.gz

function write_summary() {
    zcat "$1" | python <( echo '
import sys
print("Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure")
print("----------------------------------------------------------------------------------------------------")
for line in map(str.strip, sys.stdin):
    line = line.split("\t")
    if line[1] ==
' )
}

readarray -t cns < <(zgrep -v '^#' ${dir}/parascopy/res.matrix.bed.gz | cut -f5 | sort -n | uniq)
subdir2=${subdir}/pooled
rm -rf ${subdir2}
mkdir -p ${subdir2}

PADDING=5

for cn in "${cns[@]}"; do
    echo -e "=== Reference CN ${cn} ==="
    zgrep -v '^#' ${dir}/parascopy/res.matrix.bed.gz | awk -v cn=${cn} '$5 == cn' | cut -f1-3 | \
        bedtools merge -i - | \
        bedtools intersect -a - -b ${subdir}/calling.bed > ${subdir2}/calling_${cn}.bed

    echo -e "\n=== Filtering golden VCF ==="
    rtg vcffilter --include-bed=${subdir2}/calling_${cn}.bed -i golden_pooled.vcf.gz -o ${subdir2}/golden_${cn}.vcf.gz
    zgrep -v '^#' ${subdir2}/golden_${cn}.vcf.gz | \
        awk -v cn=${cn} -v padd=${PADDING} 'BEGIN {OFS="\t"} {
            split($10, a, "|"); if (length(a) != cn) { print $1, $2 - padd - 1, $2 + length($3) + padd - 1 }
        }' > ${subdir2}/exclude_${cn}.bed

    echo -e "\n=== Filtering variants VCF ==="
    rtg vcffilter --include-bed=${subdir2}/calling_${cn}.bed --exclude-bed=${subdir2}/exclude_${cn}.bed \
        -i ${var_dir}/variants_pooled2.vcf.gz -o ${subdir2}/parascopy_${cn}.vcf.gz

    rm -rf ${subdir2}/eval-${cn}
    rtg vcfeval \
        -b ${subdir2}/golden_${cn}.vcf.gz \
        -c ${subdir2}/parascopy_${cn}.vcf.gz \
        -t ${genome_sdf} \
        -o ${subdir2}/eval-${cn} \
        -f GQ2 \
        --sample-ploidy=${cn} > /dev/null

    rm -rf ${subdir2}/eval-${cn}-squash
    rtg vcfeval \
        -b ${subdir2}/golden_${cn}.vcf.gz \
        -c ${subdir2}/parascopy_${cn}.vcf.gz \
        -t ${genome_sdf} \
        -o ${subdir2}/eval-${cn}-squash \
        -f GQ2 \
        --squash-ploidy \
        --sample-ploidy=${cn} > /dev/null
    echo
done

${wdir}/write_summary.py ${subdir2}/eval-*
