#!/bin/bash

source /broad/software/scripts/useuse
reuse R-3.5
reuse Anaconda3
source activate /broad/hptmp/malperov/myenv
cd $1

Rscript $2 $3
