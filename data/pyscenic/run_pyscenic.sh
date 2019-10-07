#! /usr/bin/env bash
#BSUB -n 44 
#BSUB -J pyscenic
#BSUB -e err_%J.txt
#BSUB -o out_%J.txt
#BSUB -R "select[mem>450] rusage[mem=450] span[hosts=1]"

source /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load anaconda
source activate py37

python pyscenic_analysis.py
#python pyscenic_analysis_shh.py
