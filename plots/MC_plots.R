source("plots/tmp_utils.R")
source("plots/plotting.R")


task.name <<- commandArgs(trailingOnly = TRUE)[1]
results.dir <<- commandArgs(trailingOnly = TRUE)[2]
message("Starting R script to generate results")

tiss <- read.csv(paste0(results.dir, "!cells.csv"))

markers <- read.csv(paste0(results.dir, "!markers.csv"))
obj.markers <- data.frame("gene"=markers$feature, "avg_logFC"=markers$log2FC, "p_val_adj"=markers$t_qval, 
                          "cluster"=factor(markers$cluster - 1))
obj.markers <- obj.markers %>% filter(avg_logFC > 0.25)
tmp <- assignCellTypes(tiss, obj.markers, getAnnotations(tiss), record.stats = TRUE) #assign cell types
#unpack returned object
clusters <<- tmp$clusters
obj.markers <<- tmp$markers

generatePlots(tiss, task.name, clusters$cell.type) #make plots

write.csv(clusters, file=paste0(results.dir, "/!clusters.csv"))