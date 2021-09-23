#!/bin/bash

set -eu

prefix="$1"
genome="$2"
threads=10

if [[ -f "$prefix"read1.fq.gz ]]; then
    echo "Decompressing reads"
    pigz -d -p 16 "$prefix"read*.fq.gz
fi

echo "Aligning reads"
read_group="@RG\tID:sim\tSM:sim"
bwa mem "$genome" -R "$read_group" "$prefix"read{1,2}.fq -t "$threads" > "$prefix"unsort.sam
echo "Sorting alignments"
samtools sort -@ "$threads" -o "$prefix"bwa.bam "$prefix"unsort.sam
echo "Indexing alignments"
samtools index -@ "$threads" "$prefix"bwa.bam
rm "$prefix"unsort.sam

echo "Compressing reads back"
pigz -p 16 "$prefix"read*.fq

echo "Success"
