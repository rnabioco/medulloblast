#! /usr/bin/env bash
#BSUB -J prep 
#BSUB -n 6 
#BSUB -q rna
#BSUB -R "select[mem>100] rusage[mem=100] span[hosts=1]"
#BSUB -e logs/err_%J.txt
#BSUB -o logs/out_%J.txt

mkdir -p logs

. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load anaconda

unset PYTHONPATH # avoid enum error
export PYTHONNOUSERSITE=True

source activate cnmf_env

export OPENBLAS_NUM_THREADS=2

cnmf="python /beevol/home/riemondy/src/cNMF/cnmf.py"

$cnmf prepare \
    --output-dir ./full_data_cnmf \
    --name cnmf_3 \
    -c full_count_matrix.tsv \
    -k 6 7 8 9 \
    --n-iter 100 \
    --total-workers 12 \
    --seed 42 \
    --numgenes 3000

