#!/bin/bash

set -eu

threads=4
inp=""
fasta=""
out=""
sample=""
subsample=""

USAGE="$(cat <<-END
Align reads using BWA.
    -i <file>, --input <file>
        Input BAM/CRAM file.
    -o <file>, --output <file>
        Output CRAM file.
    -f <fasta>,  --fasta-ref <fasta>
        Fasta reference.
    -s <string>, --sample <string>
        Sample name.
    --subsample <seed>.<frac>
        Subsample input file, see samtools view -s.
    -@ <int>,    --threads <int>.
        Number of threads [default: $threads].
END
)"

while (( "$#" )); do
    case "$1" in
        -i|--input)
            inp="$2"
            shift 2
            ;;
        -o|--output)
            out="$2"
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
        --subsample)
            subsample="$2"
            shift 2
            ;;
        -@|--threads)
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

if [[ -z ${inp} ]] || [[ ! -f ${inp} ]]; then
    >&2 echo "Error: Input file is not provided or does not exist!"
    exit 1
elif [[ -z ${fasta} ]] || [[ ! -f ${fasta} ]]; then
    >&2 echo "Error: Reference file is not provided or does not exist!"
    exit 1
elif [[ -z ${out} ]]; then
    >&2 echo "Error: Output file is not provided!"
    exit 1
elif [[ -z ${sample} ]]; then
    >&2 echo "Error: Sample name is not provided!"
    exit 1
elif [[ ${threads} -le 0 ]]; then
    >&2 echo "Error: Cannot work with 0 threads!"
    exit 1
fi

if [[ ! -z ${subsample} ]]; then
    subsample="-s ${subsample}"
fi

header1="${out}.header1.sam"
header2="${out}.header2.sam"
samtools view -T ${fasta} -H ${inp} > ${header1}
samt_threads=$((threads / 2 - 1))
samt_threads=$((samt_threads > 0 ? samt_threads : 0))

if grep -q "@RG" ${header1}; then
    # There are read groups in the input file.
    sed 's/^\(@RG.*SM\):[^\t]\+/\1:'"${sample}/" ${header1} > ${header2}
    (cat ${header2}; samtools view -T ${fasta} -@ ${samt_threads} $subsample ${inp}) | \
        samtools view -CT ${fasta} -@ ${samt_threads} > ${out}
else
    cp ${header1} ${header2}
    echo -e "@RG\tID:${sample}\tSM:${sample}" >> ${header2}
    (
        cat ${header2};
        samtools view -T ${fasta} -@ ${samt_threads} $subsample ${inp} | sed 's/$/\tRG:Z:'"${sample}/"
    ) | samtools view -CT ${fasta} -@ ${samt_threads} > ${out}
fi

rm ${header1} ${header2}
samtools index -@ $((threads-1)) ${out}

