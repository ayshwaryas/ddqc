#!/bin/bash

source /broad/software/scripts/useuse
reuse R-3.5
reuse Anaconda3
source activate /ahg/regevdata/projects/scqc/myenv2

cd /broad/hptmp/malperov/data/kb_test/
cd $1

kb count -i ../Mus_musculus.GRCm38.cdna.all.idx -g ../transcripts_to_genes.txt -x 10xv2 -t 2 *.fastq.gz

Rscript /home/unix/malperov/method_comparison/kb/convert.R $1