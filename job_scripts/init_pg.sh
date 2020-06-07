#!/bin/bash

#$ -l h_vmem=32G
#$ -pe smp 1
#$ -binding linear:1
#$ -j y
#$ -l h_rt=72:00:00
#$ -o logs/
#$ -N pegasus_tm

source /broad/software/scripts/useuse
reuse Python-3.6
source /ahg/regevdata/projects/scqc/myenv/bin/activate
cd 

python3 pegasus/mc_functions.py mc_tm 12