#!/bin/bash

out="$1"

function try_mv {
    if [ -e "$1" ]; then
        mv "$1" "$2"
    fi
}

for d in $(ls */summary.bed); do
    d="$(dirname ${d})"
    o="$d/$out"
    if [ -d "$o" ]; then
        echo "Cannot complete, directory $o already exists"
        exit 1
    fi
    echo "Moving $d/* -> $o"
    mkdir "$o"
    try_mv "$d/depth" "$o"
    try_mv "$d/paralog_ploidy" "$o"
    try_mv "$d/extra" "$o"
    mv "$d/summary.bed" "$o"
done

