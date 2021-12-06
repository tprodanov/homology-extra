#!/bin/bash

function sum-len() {
    awk '{s += $3 - $2} END {print s}' "$@"
}

function fmt-len() {
    awk '$0~/^[0-9]+$/ { if ($0 >= 1e6) print $0/1e6,"Mb"; else if ($0 >= 1e3) print $0/1e3,"kb"; else print $0,"bp" }' "$@"
}

set -eu

wdir="$1"
genome_dir="$2"
genome="$genome_dir/genome.fa"
gen_fai="$genome_dir/genome.fa.fai"
table="$genome_dir/homology/table.bed.gz"

if [ $# -ge 3 ]; then
    n_copies="$3"
    out_dir="${wdir}/calls_${n_copies}"
else
    n_copies=0
    out_dir="${wdir}/calls"
fi

mkdir -p "$out_dir"
out_file="$out_dir/calling.bed"

echo "Writing regions for variant calling to     ...     $out_file"

zgrep -v '^#' "$wdir/parascopy/res.samples.bed.gz" | \
    awk -v n_copies="$n_copies" 'BEGIN{ FS="\t"; OFS="\t" }
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
    bedtools sort -i - -g "$gen_fai" | \
    bedtools merge -i - -d 0 > "$out_file"

echo "    Sum length = $(sum-len "$out_file" | fmt-len)"

