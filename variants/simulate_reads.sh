#!/bin/bash

set -eu

script="gen_reads.py"
out="$1"
genome="$2"

mkdir -p "$out"
if [[ -f "$out/input.vcf" ]]; then
    input_vcf="-v $out/input.vcf"
else
    input_vcf=""
fi

# python -u to skip output buffering.
command="python3 -u $(which $script) -r $genome \\
    -R 150 -c 30 -to 0 --bam --vcf --rng $RANDOM --pe 450 100 \\
    -o $out/ -tr $(dirname $out)/simulation.bed $input_vcf"
echo "$command" | tee "$out/simulation.log"
eval "$command" |& tee -a "$out/simulation.log"

rename "_" "" "$out/_"*

echo
"$(dirname "$0")"/bwa_wrap.sh "$out/" "$genome"
