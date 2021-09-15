#!/bin/bash

# Stdin - list of .bed.gz files,
# Stdout - subset of the processed files,
# Command arguments - regions.

tabix_extra() {
    sample="$(basename "$1" | sed 's/\.qm2.*//')"
    tabix "$@" | sed "s/$/\t$sample/"
}

export -f tabix_extra

cat - | parallel -P 10 --eta tabix_extra {} "$@"
