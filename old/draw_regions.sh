#!/bin/bash

set -e

in="$1"
out="$2"

mkdir -p "$out"

cat "$in" | while IFS=$'\t' read -r chrom start end name
do
    url='http://genome.ucsc.edu/cgi-bin/hgRenderTracks?position='
    url+="${chrom}:$((start+1))-${end}"
    # Previously: refGene
    url+='&hideTracks=1&wgEncodeGencodeBasicV37=pack&genomicSuperDups=pack&pubs=pack'
    echo "$name"
    curl --silent "$url" > "${out}/${name}.png"
done
