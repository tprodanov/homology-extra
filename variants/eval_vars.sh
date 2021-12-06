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

var_dir="$(dirname $(ls ${dir}/parascopy/*/variants.vcf.gz | head -1))"

bcftools view -T ${subdir}/calling.bed -i'GQ >= 0' ${var_dir}/variants.vcf.gz | bgzip > ${subdir}/parascopy.vcf.gz
tabix -p vcf ${subdir}/parascopy.vcf.gz

if ! [[ -f ${dir}/golden_phased.vcf.gz ]]; then
    ${wdir}/reformat_golden.sh ${dir}
fi
bcftools view -T ${subdir}/calling.bed golden_phased.vcf.gz | bgzip > ${subdir}/golden.vcf.gz
tabix -p vcf ${subdir}/golden.vcf.gz

cd ${subdir}
rm -rf eval-*
rtg vcfeval -b golden.vcf.gz -c freebayes.vcf.gz -t ~/Data/hg38/genome/genome.fa.sdf -o eval-fb > /dev/null
rtg vcfeval -b golden.vcf.gz -c gatk.vcf.gz -t ~/Data/hg38/genome/genome.fa.sdf -o eval-gatk > /dev/null
rtg vcfeval -b golden.vcf.gz -c parascopy.vcf.gz -t ~/Data/hg38/genome/genome.fa.sdf -o eval-paras > /dev/null

${wdir}/write_summary.py eval-*
