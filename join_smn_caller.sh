#!/bin/bash

echo 'Joining genomes/*.tsv -> joined.tsv'
out_tsv="./joined.tsv"
head -n1 "$(ls genomes/*.tsv | head -n1)" > "$out_tsv"
tail -q -n+2 genomes/*.tsv >> "$out_tsv"

echo 'Joining genomes/*.json -> joined.json'
out_json="./joined.json"
echo "{" > "$out_json"
awk '{ print (NR == FNR ? "" : ","), substr($0, 2, length($0) - 2) }' genomes/*.json >> "$out_json"
echo "}" >> "$out_json"

