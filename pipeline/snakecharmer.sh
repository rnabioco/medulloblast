#!/usr/bin/env bash
#BSUB -J 10xscRNA 
#BSUB -o logs/10x_%J.out
#BSUB -e logs/10x_%J.err
#BSUB -R "select[mem>4] rusage[mem=4] " 

set -o nounset -o pipefail -o errexit -x

mkdir -p logs

. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load cellranger/3.0.2

args=' 
  -q rna 
  -o {log}.out 
  -e {log}.err 
  -J {params.job_name} 
  -R "{params.memory} span[hosts=1] " 
  -n {threads} ' 


#snakemake \
#  --drmaa "$args" \
#  --snakefile Snakefile \
#  --jobs 2 \
#  --latency-wait 50 \
#  --rerun-incomplete  \
#  --configfile config_181114.yaml 
#
#snakemake \
#  --drmaa "$args" \
#  --snakefile Snakefile \
#  --jobs 2 \
#  --latency-wait 50 \
#  --rerun-incomplete  \
#  --configfile config_190124.yaml 
#
#snakemake \
#  --drmaa "$args" \
#  --snakefile Snakefile \
#  --jobs 2 \
#  --latency-wait 50 \
#  --rerun-incomplete  \
#  --configfile config_190315.yaml 
#
#snakemake \
#  --drmaa "$args" \
#  --snakefile Snakefile \
#  --jobs 2 \
#  --latency-wait 50 \
#  --rerun-incomplete  \
#  --configfile config_190401.yaml 
#
#snakemake \
#  --drmaa "$args" \
#  --snakefile Snakefile \
#  --jobs 2 \
#  --latency-wait 50 \
#  --rerun-incomplete  \
#  --configfile config_190403.yaml 

#snakemake \
#  --drmaa "$args" \
#  --snakefile Snakefile \
#  --jobs 2 \
#  --latency-wait 50 \
#  --rerun-incomplete  \
#  --configfile config_190429.yaml 
#
#
#snakemake \
#  --drmaa "$args" \
#  --snakefile Snakefile \
#  --jobs 2 \
#  --latency-wait 50 \
#  --rerun-incomplete  \
#  --configfile config_190415.yaml
#
#
#snakemake \
#  --drmaa "$args" \
#  --snakefile Snakefile \
#  --jobs 2 \
#  --latency-wait 50 \
#  --rerun-incomplete  \
#  --configfile config_190104.yaml
#
#snakemake \
#  --drmaa "$args" \
#  --snakefile Snakefile \
#  --jobs 3 \
#  --latency-wait 50 \
#  --rerun-incomplete  \
#  --configfile config_190503.yaml

snakemake \
  --drmaa "$args" \
  --snakefile Snakefile \
  --jobs 3 \
  --latency-wait 50 \
  --rerun-incomplete  \
  --configfile config_190510.yaml
