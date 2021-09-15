#!/bin/bash

function sum-len() {
    awk '{s += $3 - $2} END {print s}' "$@"
}

function fmt-len() {
    awk '$0~/^[0-9]+$/ { if ($0 >= 1e6) print $0/1e6,"Mb"; else if ($0 >= 1e3) print $0/1e3,"kb"; else print $0,"bp" }' "$@"
}

set -eu

input="$1"
genome_dir="$2"
output="$3"

genome="$genome_dir/genome.fa"
gen_fai="$genome_dir/genome.fa.fai"
table="$genome_dir/homology/table.bed.gz"
distance=1000

mkdir -p "$output"
out_file="$output/simulation.bed"

echo "Writing regions for read simulation to     ...     $out_file"

tabix "$table" -R "$input/bed/all.bed" | \
    awk 'BEGIN{FS="\t"; OFS="\t"} { print $1,$2,$3; if ($8 != "tangled") { print $5,$6,$7 } }' | \
    bedtools sort -i - -g "$gen_fai" | \
    bedtools slop -i - -g "$gen_fai" -b "$distance" | \
    bedtools merge -i - -d "$distance" | \
    bedtools sort -i - -g "$gen_fai" > "$out_file"

echo "    Sum length = $(sum-len "$out_file" | fmt-len)"
