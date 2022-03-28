#!/bin/bash

set -eu

USAGE="$(cat <<-END
Align reads using BWA.
    -b <path>,   --bwa <path>
        BWA executable [default: bwa].
    -i <prefix>, --input <prefix>
        Input prefix. There must be files <prefix>{1,2}.fq[.gz]
    -o <prefix>, --output <prefix>
        Output prefix.
    -f <fasta>,  --fasta-ref <fasta>
        Fasta reference.
    -s <string>, --sample <string>
        Sample name [default: none].
    -t <int>,    --threads <int>.
        Number of threads [default: 4].
END
)"

bwa=bwa
threads=4
sample=""

while (( "$#" )); do
    case "$1" in
        -b|--bwa)
            bwa="$2"
            shift 2
            ;;
        -i|--input)
            in_prefix="$2"
            shift 2
            ;;
        -o|--output)
            out_prefix="$2"
            shift 2
            ;;
        -f|--fasta-ref)
            fasta="$2"
            shift 2
            ;;
        -s|--sample)
            sample="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
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

compress=false
if [[ -f ${in_prefix}1.fq.gz ]]; then
    echo "Decompressing reads"
    pigz -d -p ${threads} ${in_prefix}{1,2}.fq.gz
    compress=true
fi

echo "Aligning reads"
if [ -z ${sample} ]; then
    bwa mem ${fasta} -R "@RG\tID:${sample}\tSM:${sample}" \
        ${in_prefix}{1,2}.fq -t ${threads} > ${out_prefix}.unsort.sam
else
    bwa mem ${fasta} ${in_prefix}{1,2}.fq -t ${threads} > ${out_prefix}.unsort.sam
fi

echo "Sorting alignments"
samtools sort -@ ${threads} -o ${out_prefix}.bwa.bam ${out_prefix}.unsort.sam

echo "Indexing alignments"
samtools index -@ ${threads} ${out_prefix}.bwa.bam
rm ${out_prefix}.unsort.sam

echo "Compressing reads back"
if [[ ${compress} = true ]]; then
    pigz -p ${threads} ${in_prefix}{1,2}.fq
fi

echo "Success"
