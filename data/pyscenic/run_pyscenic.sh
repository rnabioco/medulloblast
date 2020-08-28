#! /usr/bin/env bash
#BSUB -n 24
#BSUB -J pyscenic
#BSUB -e logs/err_%J.txt
#BSUB -o logs/out_%J.txt
#BSUB -R "select[mem>450] rusage[mem=450] span[hosts=1]"
#BSUB -q rna

mkdir -p logs

source /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load anaconda
source activate py37

python pyscenic_analysis.py "shh_expr_matrix.tsv.gz" "output_shh"
python pyscenic_analysis.py "gp3_expr_matrix.tsv.gz" "output_gp3_only"
python pyscenic_analysis.py "gp4_expr_matrix.tsv.gz" "output_gp4_only"
python pyscenic_analysis.py "gp34_expr_matrix.tsv.gz" "output_gp34"
python pyscenic_analysis.py "immune_expr_matrix.tsv.gz" "output_immune"

