#!/bin/bash

#$ -l h_vmem=15G
#$ -pe smp 6
#$ -binding linear:6
#$ -j y
#$ -l h_rt=48:00:00
#$ -o /ahg/regevdata/projects/scqc/code/logs/

source /broad/software/scripts/useuse
use .python-3.8.3
use R-3.5

source /home/unix/malperov/myenv3.8/bin/activate

cd /ahg/regevdata/projects/scqc/code
python run_jc_old.py $1 $2

