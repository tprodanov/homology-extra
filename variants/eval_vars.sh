#!/bin/bash

set -e

wdir="$(dirname "$0")"

USAGE="$(cat <<-END
Evaluate diploid genotypes.
    -a <file>,   --alignment <file>
        Alignment file (BAM or CRAM).
    -g <file>,   --golden <file>
        Golden vcf.gz file.
    -p <dir>,   --parascopy <dir>
        Parascopy directory.
    -P <str>,   --parascopy-name <str>
        Parascopy name. By default, use basename of the Parascopy dir.
    -f <file>,  --fasta-ref <file>
        Fasta reference. Must contain <file>.sdf index.
    -t <path> <path>,   --tools <path> <path>
        Freebayes & GATK paths. Use "freebayes" and "gatk" if not set.
    -r,  --rerun
        Rerun everything.
    -o <dir>,   --output <dir>
        Output directory. Must contain calling.bed file.
END
)"

alignment=""
freebayes="freebayes"
gatk="gatk"
rerun=false
command="$0 $*"

while (( "$#" )); do
    case "$1" in
        -a|--alignment)
            alignment="$2"
            shift 2
            ;;
        -g|--golden)
            golden="$2"
            shift 2
            ;;
        -p|--parascopy)
            parascopy="$2"
            shift 2
            ;;
        -P|--parascopy-name)
            par_name="$2"
            shift 2
            ;;
        -f|--fasta-ref)
            genome="$2"
            shift 2
            ;;
        -o|--output)
            output="$2"
            shift 2
            ;;
        -t|--tools)
            freebayes="$2"
            gatk="$3"
            shift 3
            ;;
        -r|--rerun)
            rerun=true
            shift 1
            ;;
        -h|--help)
            echo "${USAGE}"
            exit 0
            ;;
        *)
            >&2 echo "Error: Unexpected argument $1"
            exit 1
            ;;
    esac
done

if [[ -z ${golden} ]]; then
    >&2 echo "Error: Golden VCF file (-g, --golden) is not provided!"
    exit 1
elif [[ -z ${parascopy} ]]; then
    >&2 echo "Error: Parascopy directory (-p, --parascopy) is not provided!"
    exit 1
elif [[ -z ${genome} ]]; then
    >&2 echo "Error: Fasta reference (-f, --fasta-ref) is not provided!"
    exit 1
elif [[ -z ${output} ]]; then
    >&2 echo "Error: Output directory (-o, --output) is not provided!"
    exit 1
fi

if [[ -z ${par_name} ]]; then
    par_name="$(basename ${parascopy})"
fi

set -u

mkdir -p ${output}
regions="${output}/calling.bed"
log="${output}/eval.log"
echo "Command: $command" > ${log}

if [[ ! -f ${regions} ]]; then
    >& 2 echo "** Generating variant calling regions"
    ${wdir}/get_call_regions.sh -i ${parascopy} -f ${genome} -o ${regions} &>> ${log}
fi

if [[ ${rerun} = true ]]; then
    rm -f ${output}/{golden,${par_name},freebayes,gatk}.vcf.gz{,.tbi}
fi

if [[ ! -f ${output}/golden.vcf.gz ]]; then
    >&2 echo "** Filtering Golden VCF file."
    rtg vcffilter --include-bed=${regions} -i ${golden} -o ${output}/golden.vcf.gz &>> ${log}
else
    >&2 echo "!! Skipping Golden VCF file."
fi

if [[ ! -f ${output}/${par_name}.vcf.gz ]]; then
    >&2 echo "** Filtering Parascopy VCF file."
    rtg vcffilter --include-bed=${regions} -i ${parascopy}/variants.vcf.gz -o ${output}/${par_name}.vcf.gz &>> ${log}
else
    >&2 echo "!! Skipping Parascopy VCF file."
fi

if [[ -n ${alignment} ]] && [[ ! -f ${output}/freebayes.vcf.gz ]]; then
    >&2 echo "** Running Freebayes."
    ${freebayes} -f ${genome} -t ${regions} --genotype-qualities ${alignment} 2>> ${log} > ${output}/freebayes.vcf
    bgzip ${output}/freebayes.vcf
    tabix -p vcf ${output}/freebayes.vcf.gz
else
    >&2 echo "!! Skipping Freebayes."
fi

if [[ -n ${alignment} ]] && [[ ! -f ${output}/gatk.vcf.gz ]]; then
    >&2 echo "** Running GATK."
    ${gatk} HaplotypeCaller -I ${alignment} -R ${genome} -L ${regions} -O ${output}/gatk.vcf.gz &>> ${log}
else
    >&2 echo "!! Skipping GATK."
fi

if [[ ${rerun} = true ]]; then
    rm -rf ${output}/eval-{${par_name},freebayes,gatk}
fi

if [[ -f ${output}/gatk.vcf.gz ]] && [[ ! -d ${output}/eval-freebayes ]]; then
    >&2 echo "** Running Freebayes vcfeval."
    rtg vcfeval -b ${output}/golden.vcf.gz -c ${output}/freebayes.vcf.gz \
        -t ${genome}.sdf -e ${regions} -o ${output}/eval-freebayes &>> ${log}
else
    >&2 echo "!! Skipping Freebayes vcfeval."
fi

if [[ -f ${output}/gatk.vcf.gz ]] && [[ ! -d ${output}/eval-gatk ]]; then
    >&2 echo "** Running GATK vcfeval."
    rtg vcfeval -b ${output}/golden.vcf.gz -c ${output}/gatk.vcf.gz \
        -t ${genome}.sdf -e ${regions} -o ${output}/eval-gatk &>> ${log}
else
    >&2 echo "!! Skipping GATK vcfeval."
fi

if [[ ! -d ${output}/eval-${par_name} ]]; then
    >&2 echo "** Running Parascopy vcfeval."
    rtg vcfeval -b ${output}/golden.vcf.gz -c ${output}/${par_name}.vcf.gz \
        -t ${genome}.sdf -e ${regions} -o ${output}/eval-${par_name} &>> ${log}
else
    >&2 echo "!! Skipping Parascopy vcfeval."
fi

printf "\n\n"
${wdir}/write_summary.py ${output}/eval-*
