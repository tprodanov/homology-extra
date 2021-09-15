#!/bin/bash

set -e

if [ "$#" -ne 3 ]; then
    >&2 echo "Incorrect number of arguments"
    exit 1
fi

run_dir="$(dirname "${BASH_SOURCE[0]}")"
data_dir="$1"
plots="$2"
suffix="$3"

>&2 echo "Drawing copy number variations"
mkdir -p ${plots}/${suffix}/depth
ls -d ${data_dir}/./*/${suffix}/extra | xargs -P11 -i ${run_dir}/final_plots/cnv.r -i "{}" -o ${plots}/${suffix}/depth -w 20
>&2 echo

#exit

>&2 echo "Drawing paralog ploidy profiles"
mkdir -p ${plots}/${suffix}/paralog
ls -d ${data_dir}/./*/${suffix}/extra | xargs -P10 -i ${run_dir}/psv_sample_matrix.r "{}" ${plots}/${suffix}/paralog

