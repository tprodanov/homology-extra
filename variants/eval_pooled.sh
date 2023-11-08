#!/bin/bash

set -eu
wdir="$(dirname "$0")"

USAGE="$(cat <<-END
Evaluate pooled genotypes.
    -p <dir>,   --parascopy  <dir>
        Parascopy directory.
    -o <dir>,   --output <dir>
        Output directory.
    -f <file>,   --fasta-ref <file>
        Fasta reference. Must contain <file>.sdf index.
    --padding <int>
        Padding around the benchmark variants with unexpected copy number [default: 5]
END
)"

padding=5

while (( "$#" )); do
    case "$1" in
        -p|--parascopy)
            input="$2"
            shift 2
            ;;
        -o|--output)
            output="$2"
            shift 2
            ;;
        -f|--fasta-ref)
            genome="$2"
            shift 2
            ;;
        --padding)
            padding="$2"
            shift 2
            ;;
        -h|--help)
            echo "${USAGE}"
            exit 0
            ;;
        *)
            echo "Error: Unexpected argument $1" >&2
            exit 1
            ;;
    esac
done

input="$(dirname $(ls ${input}/**/variants_pooled.vcf.gz | head -1))"

readarray -t cns < <(zgrep -v '^#' ${input}/res.matrix.bed.gz | cut -f5 | sort -n | uniq)
mkdir -p ${output}

for cn in "${cns[@]}"; do
    echo -e "=== Reference CN ${cn} ==="
    zgrep -v '^#' ${input}/res.matrix.bed.gz | awk -v cn=${cn} '$5 == cn' | cut -f1-3 | \
        bedtools merge -i - > ${output}/calling_${cn}.bed

    echo -e "\n=== Filtering benchmark VCF ==="
    rm -f ${output}/benchmark_${cn}.vcf.gz{,.tbi}
    rtg vcffilter --include-bed=${output}/calling_${cn}.bed -i benchmark_pooled.vcf.gz \
        -o ${output}/benchmark_${cn}.vcf.gz
    zgrep -v '^#' ${output}/benchmark_${cn}.vcf.gz | \
        awk -v cn=${cn} -v padding=${padding} 'BEGIN {OFS="\t"} {
            split($10, a, "|"); if (length(a) != cn) { print $1, $2 - padding - 1, $2 + length($3) + padding - 1 }
        }' > ${output}/exclude_${cn}.bed

    echo -e "\n=== Filtering variants VCF ==="
    rm -f ${output}/parascopy_${cn}.vcf.gz{,.tbi}
    rtg vcffilter --include-bed=${output}/calling_${cn}.bed --exclude-bed=${output}/exclude_${cn}.bed \
        -i ${input}/variants_pooled.vcf.gz -o ${output}/parascopy_${cn}.vcf.gz

    rm -rf ${output}/eval-${cn}
    rtg vcfeval \
        -b ${output}/benchmark_${cn}.vcf.gz \
        -c ${output}/parascopy_${cn}.vcf.gz \
        -t ${genome}.sdf \
        -o ${output}/eval-${cn} \
        --sample-ploidy=${cn} > /dev/null

    rm -rf ${output}/eval-${cn}-squash
    rtg vcfeval \
        -b ${output}/benchmark_${cn}.vcf.gz \
        -c ${output}/parascopy_${cn}.vcf.gz \
        -t ${genome}.sdf \
        -o ${output}/eval-${cn}-squash \
        --squash-ploidy \
        --sample-ploidy=${cn} > /dev/null
    echo
done

${wdir}/write_summary.py ${output}/eval-*
