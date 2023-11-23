#!/bin/bash

set -eu
wdir="$(dirname "$0")"

USAGE="$(cat <<-END
Evaluate pooled genotypes.
    -b/--benchmark <vcf> <bed>
        Benchmark pooled variants and filtering BED file.
    -c/--calls     <vcf> <bed>
        Input pooled variants and corresponding BED file.
    -d/--discard   <bed>
        Optionally, discard regions from this file.
    -e/--examine   <bed>
        Input BED file with region copy numbers.
    -o/--output    <dir>
        Output directory.
    -f/--fasta-ref <file>
        Fasta reference. Must contain <file>.sdf index.
    --padding      <int>
        Padding around variants with unexpected copy number [default: 5]
    --cns          <int,>
        List of copy numbers through comma [default: 4,6,8].
    --memory       <str>
        Available memory [10G].
END
)"

padding=5
incns="4,6,8"
discard=/dev/null
mem=10G

while (( "$#" )); do
    case "$1" in
        -c|--calls)
            calls_vcf="$2"
            calls_bed="$3"
            shift 3
            ;;
        -b|--benchmark)
            bench_vcf="$2"
            bench_bed="$3"
            shift 3
            ;;
        -e|--examine)
            examine="$2"
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
        -d|--discard)
            discard="$2"
            shift 2
            ;;
        --padding)
            padding="$2"
            shift 2
            ;;
        --cns)
            incns="$2"
            shift 2
            ;;
        --memory)
            mem="$2"
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

if [[ ! -f ${bench_vcf-} ]] || [[ ! -f ${bench_bed-} ]]; then
    >&2 echo "Error: Benchmark VCF or BED file is not provided or does not exist!"
    exit 1
elif [[ ! -f ${calls_vcf-} ]] || [[ ! -f ${calls_bed-} ]]; then
    >&2 echo "Error: Variant calls VCF or BED file is not provided or does not exist!"
    exit 1
elif [[ ! -f ${genome-} ]]; then
    >&2 echo "Error: Fasta reference is not provided or does not exist!"
    exit 1
elif [[ ! -f ${examine-} ]]; then
    >&2 echo "Error: Examine BED file is not provided or does not exist!"
    exit 1
elif [[ -z ${output} ]]; then
    >&2 echo "Error: Output directory (-o, --output) is not provided!"
    exit 1
fi

readarray -t cns < <(echo "$incns" | tr ',' '\n')
mkdir -p ${output}

for cn in "${cns[@]}"; do
    subdir="${output}/${cn}"
    rm -rf "${subdir}"
    mkdir -p "${subdir}"

    echo -e "=== Reference CN ${cn} ==="
    zgrep -v '^#' ${examine} | awk -v cn=${cn} '$4 == "PASS" && $5 == cn' | cut -f1-3 | \
        bedtools sort -i - -g "${genome}.fai" | bedtools merge -i - > "${subdir}/approp_cn.bed"
    rtg vcffilter --include-bed="${subdir}/approp_cn.bed" -i "${bench_vcf}" -o "${subdir}/bench.vcf.gz"
    rtg vcffilter --include-bed="${subdir}/approp_cn.bed" -i "${calls_vcf}" -o "${subdir}/calls.vcf.gz"
    zgrep -v '^#' "${subdir}/bench.vcf.gz" | cut -f1,2,4,10 | cut -f1 -d: | sed 's,/,|,g' | \
        awk -v cn=${cn} -v padd=${padding} \
                'BEGIN {OFS="\t"} { if (split($4, a, "|") != cn) { print $1,$2-1,$2-1+length($3)+padd; } }' | \
        bedtools sort -i - -g "${genome}.fai" | bedtools merge -i - > "${subdir}/excl_bench.bed"

    bedtools intersect -a "${subdir}/approp_cn.bed" -b "${bench_bed}" | \
        bedtools intersect -a - -b "${calls_bed}" | \
        bedtools subtract -a - -b "${subdir}/excl_bench.bed" | \
        bedtools subtract -a - -b "${discard}" | \
        bedtools sort -i - -g "${genome}.fai" | bedtools merge -i - > "${subdir}/eval.bed"

    readarray -t reg_stats < <(awk '{ n += 1; l += $3 - $2 } END { printf("%d\n%d\n",n,l) }' "${subdir}/eval.bed")
    if [[ "${reg_stats[0]}" -eq 0 ]]; then
        echo -e "\n=== No regions to analyze, continuing ===\n"
        continue;
    fi
    echo -e "\n=== Examining ${reg_stats[0]} regions, in total ${reg_stats[1]} bp ==="

    if [[ -f "${subdir}/eval/done" ]]; then
        echo "Skipping evaluation"
    else
        rtg "RTG_MEM=${mem}" vcfeval \
            -e "${subdir}/eval.bed" \
            --decompose \
            -b "${subdir}/bench.vcf.gz" \
            -c "${subdir}/calls.vcf.gz" \
            -t "${genome}.sdf" \
            -o "${subdir}/eval" \
            --sample-ploidy=${cn} &> "${subdir}/eval.log"
    fi
    ${wdir}/write_summary.py "${subdir}/eval"

    echo -e "\n=== Evaluating squashed variants ==="
    if [[ -f "${subdir}/eval_squash/done" ]]; then
        echo "Skipping evaluation"
    else
        rtg "RTG_MEM=${mem}" vcfeval \
            -e "${subdir}/eval.bed" \
            --decompose \
            -b "${subdir}/bench.vcf.gz" \
            -c "${subdir}/calls.vcf.gz" \
            -t "${genome}.sdf" \
            -o "${subdir}/eval_squash" \
            --squash-ploidy \
            --sample-ploidy=${cn} &> "${subdir}/eval_squash.log"
    fi
    ${wdir}/write_summary.py "${subdir}/eval_squash"
    echo
done
