#!/bin/bash

set -eu

function sum-len() {
    awk '{s += $3 - $2} END {print s}' "$@"
}

function fmt-len() {
    awk '$0~/^[0-9]+$/ { if ($0 >= 1e6) print $0/1e6,"Mb"; else if ($0 >= 1e3) print $0/1e3,"kb"; else print $0,"bp" }' "$@"
}

wdir="$(dirname "$0")"

n_copies=0

USAGE="$(cat <<-END
Extract regions for variant calling.
    -p <dir>,   --parascopy  <dir>
        Parascopy directory.
    -o <file>,   --output <file>
        Output file.
    -f <file>,  --fasta-ref <file>
        Reference fasta file.
    -n <int>,   --n-copies <int>
        Take only duplications with <int> copies [default: all].
END
)"

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
            fasta_ref="$2"
            shift 2
            ;;
        -n|--n-copies)
            n_copies="$2"
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

mkdir -p $(dirname ${output})

echo "Writing regions for variant calling to     ...     ${output}"

zgrep -v '^#' ${input}/res.samples.bed.gz | \
    awk -v n_copies=${n_copies} 'BEGIN{ FS="\t"; OFS="\t" }
    {
        if ($13 != "*") {
            split($13, hom_regions, ",");
        }
        curr_n_copies = length(hom_regions) + 1;
        if (n_copies > 0 && curr_n_copies != n_copies) {
            next;
        }

        print $1,$2,$3;
        for (i in hom_regions) {
             region = hom_regions[i];
             print gensub(/([^:]+):([0-9]+)-([0-9]+):.*/, "\\1\t\\2\t\\3", "g", region);
        }
    }' | \
    bedtools sort -i - -g ${fasta_ref}.fai | \
    bedtools merge -i - -d 0 > ${output}

echo "    Sum length = $(sum-len ${output} | fmt-len)"
