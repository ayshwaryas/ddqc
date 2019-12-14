#!/bin/bash

#$ -l h_vmem=32G
#$ -pe smp 2
#$ -binding linear:2
#$ -j y
#$ -l h_rt=12:00:00
#$ -o logs/

source /broad/software/scripts/useuse
reuse R-3.5
reuse Anaconda3
source activate /broad/hptmp/malperov/myenv

R_PROJECT=$1
export R_PROJECT

Rscript scripts/mc.R
