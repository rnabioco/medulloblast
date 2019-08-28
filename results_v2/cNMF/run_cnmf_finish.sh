#! /usr/bin/env bash
#BSUB -J finsh
#BSUB -n 6
#BSUB -q rna
#BSUB -R "select[mem>50] rusage[mem=50] span[hosts=1]"
#BSUB -e logs/err_%J.txt
#BSUB -o logs/out_%J.txt


. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load anaconda

unset PYTHONPATH # avoid enum error
export PYTHONNOUSERSITE=True

source activate cnmf_env

export OPENBLAS_NUM_THREADS=6


cnmf="python $HOME/src/cNMF/cnmf.py"

#$cnmf combine --output-dir ./full_data_cnmf --name cnmf_3
#$cnmf k_selection_plot --output-dir ./full_data_cnmf --name cnmf_3
$cnmf consensus \
    --output-dir ./full_data_cnmf \
    --name cnmf_2 \
    --components 12 \
    --local-density-threshold 0.01 \
    --show-clustering

