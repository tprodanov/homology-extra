#!/bin/bash

set -e

in="./genes.bed"
log="./log.txt"
out1="./genes_y.bed"
out2="./genes_n.bed"
plots="./plots"
view="/home/timofey/Code/HomologyTools/view.py"
table="../genome/homology/table.bed.gz"
out_log="./manual.log"

function write_int() {
    echo "$1" | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'
}

if [ $# -eq 1 ]; then
    start_gene="$1"
elif [ -f "$out1" ] || [ -f "$out2" ]; then
    echo "Output files already exist"
    exit 1
else
    start_gene="NONE"
fi

i=0
cat "$in" | while IFS=$'\t' read -r chrom start end name
do
    i=$((i + 1))
    if [ "$start_gene" != "NONE" ] && [ "$start_gene" != "$name" ]; then
        continue
    fi
    start_gene="NONE"

    bed_entry=$(echo -en "${chrom}\t${start}\t${end}\t${name}")
    if [ "$chrom" = "chrX" ] || [ "$chrom" = "chrY" ]; then
        echo "$bed_entry" >> "$out2"
        continue
    fi

    echo "==================="
    region="${chrom}:$(write_int $((start+1)))-$(write_int $end)"
    echo "[$i] $name   $region   ($(write_int $((end - start))) bp)"
    echo

    viewnior "${plots}/${name}.png" 2> /dev/null &

    echo -e "------- Log -------"
    grep -P "\b${name}\b" "$log" || true
    echo

    echo -e "------ Table ------"
    "$view" "$table" -r "$region" --pretty -i "ALENGTH >= 2000" | cut -f1-14 | csvr --skip-header -w 20
    echo

    echo -e "------ Keep? ------"
    read response < /dev/tty
    while ! (echo "$response" | grep -q "^[myn]"); do
        echo "Unrecognized answer $response"
        read response < /dev/tty
    done

    echo "Responded $response"
    if [ "$response" = "y" ]; then
        echo "$bed_entry" >> "$out1"
    elif [ "$response" = "n" ]; then
        echo "$bed_entry" >> "$out2"
    else
        read -r flag new_start new_end <<< "$response"
        echo "Response flag = $flag  start = $new_start  end = $new_end"
        if [ "$new_start" = "." ]; then
            new_start="$start"
        else
            new_start="$(echo "$new_start" | sed 's/,//g')"
        fi
        if [ "$new_end" = "." ]; then
            new_end="$end"
        else
            new_end="$(echo "$new_end" | sed 's/,//g')"
        fi
        echo "Moving region ${region} -> ${chrom}:$(write_int $((new_start+1)))-$(write_int $new_end)"
        echo -e "${chrom}\t${new_start}\t${new_end}\t${name}" >> "$out1"
    fi
    kill %1
    echo
done | tee -a "$out_log"
