#!/bin/bash

source /broad/software/scripts/useuse
reuse R-3.5
source ~/myenv/bin/activate
cd $1

Rscript $2 $3