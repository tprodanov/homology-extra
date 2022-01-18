#!/bin/bash

set -eu

wdir="$(dirname "$0")"

USAGE="$(cat <<-END
Evaluate diploid genotypes.
    -i <dir>,   --input <dir>
        Input directory. Must contain golden_phased.vcf.gz and bwa.bam.
    -p <dir>,   --parascopy <dir>
        Parascopy directory.
    -g <file>,  --golden <file>
        Golden VCF file.
    -o <dir>,   --output <dir>
        Output directory. Must contain calling.bed file.
    -f <file>,   --fasta-ref <file>
        Fasta reference. Must contain <file>.sdf index.
END
)"

while (( "$#" )); do
    case "$1" in
        -i|--input)
            input="$2"
            shift 2
            ;;
        -p|--parascopy)
            par_dir="$2"
            shift 2
            ;;
        -o|--output)
            output="$2"
            shift 2
            ;;
        -f|--fasta-ref)
            genome="$2"
            shift 2
            ;;
        -h|--help)
            echo "${USAGE}"
            exit 0
            ;;
        *)
            echo "Error: Unexpected argument $1" >&2
            exit 1
            ;;
    esac
done

par_dir="$(dirname $(ls ${par_dir}/**/variants.vcf.gz | head -1))"

rm -f ${output}/{parascopy,golden}.vcf.gz{,.tbi}
rtg vcffilter --include-bed=${output}/calling.bed -i ${input}/golden_phased.vcf.gz -o ${output}/golden.vcf.gz
rtg vcffilter --include-bed=${output}/calling.bed -i ${par_dir}/variants.vcf.gz -o ${output}/parascopy.vcf.gz

~/Code/freebayes/build/freebayes -f ${genome} -t ${output}/calling.bed --genotype-qualities ${input}/bwa.bam | \
    bgzip > ${output}/freebayes.vcf.gz
tabix -p vcf ${output}/freebayes.vcf.gz
~/Code/gatk-4.2.2.0/gatk HaplotypeCaller -I ${input}/bwa.bam -R ${genome} \
    -L ${output}/calling.bed -O ${output}/gatk.vcf.gz

rm -rf ${output}/eval-{fb,gatk,paras}
rtg vcfeval -b ${output}/golden.vcf.gz -c ${output}/freebayes.vcf.gz \
    -t ${genome}.sdf -o ${output}/eval-fb > /dev/null
rtg vcfeval -b ${output}/golden.vcf.gz -c ${output}/gatk.vcf.gz \
    -t ${genome}.sdf -o ${output}/eval-gatk > /dev/null
rtg vcfeval -b ${output}/golden.vcf.gz -c ${output}/parascopy.vcf.gz \
    -t ${genome}.sdf -o ${output}/eval-paras > /dev/null

${wdir}/write_summary.py ${output}/eval-{fb,gatk,paras}
