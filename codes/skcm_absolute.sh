#!/bin/bash
#SBATCH -N 2
#SBATCH -n 4
#SBATCH -p cn
source $HOME/miniconda3/bin/activate
echo "start absolute.R" `date`
Rscript $HOME/R/rscripts/skcm_absolute.R
echo "end absolute.R" `date`
