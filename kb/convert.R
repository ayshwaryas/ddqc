directory <<- commandArgs(trailingOnly = TRUE)[1]
source("scripts/settings.R")
source("scripts/local_settings.R")
library(readr)

data.path <- paste0(data.dir, directory, "/counts_unfiltered/")
output.path <- paste0(data.dir, "results/", directory, ".rds")

mtx <- t(readMM(paste0(data.path, "matrix.mtx")))
genes=read_tsv(paste0(data.path, "genes.tsv"), col_names=F)
barcodes=read_tsv(paste0(data.path, "barcodes.tsv"), col_names=F)
rownames(mtx)=genes$X1
colnames(mtx)=barcodes$X1
tiss.kb <- CreateSeuratObject(mtx) 

saveRDS(tiss, output.path)

