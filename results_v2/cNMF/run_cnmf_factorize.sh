#! /usr/bin/env bash
#BSUB -J factorize[1-4]
#BSUB -n 16
#BSUB -q rna 
#BSUB -R "select[mem>20] rusage[mem=20] span[hosts=1]"
#BSUB -e logs/err_%J.txt
#BSUB -o logs/out_%J.txt

mkdir -p logs

. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load anaconda

unset PYTHONPATH # avoid enum error
export PYTHONNOUSERSITE=True

source activate cnmf_env

export OPENBLAS_NUM_THREADS=16

workers=($(seq 0 1 4))
idx=${workers[$(($LSB_JOBINDEX - 1))]}


cnmf="python /beevol/home/riemondy/src/cNMF/cnmf.py"

$cnmf \
    factorize \
    --output-dir ./full_data_cnmf \
    --name cnmf_3 \
    --worker-index $idx 

