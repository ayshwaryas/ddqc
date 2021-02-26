#!/bin/bash

#$ -l h_vmem=50G
#$ -pe smp 1
#$ -binding linear:1
#$ -j y
#$ -l h_rt=72:00:00
#$ -o logs/
#$ -t 1-140

source /broad/software/scripts/useuse
reuse R-3.5
reuse .python-3.8.3
source /home/unix/malperov/myenv3.8/bin/activate

python3 mc.py "mc_hca" $SGE_TASK_ID