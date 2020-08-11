data.from.pg <- TRUE
remake.plots <- FALSE

source("../scripts/mc_functions.R")
source("../scripts/readers.R")
source("../scripts/settings.R")
source("../scripts/local_settings.R")
project <<- commandArgs(trailingOnly = TRUE)[1]
tissue <<- commandArgs(trailingOnly = TRUE)[2]
res <<- commandArgs(trailingOnly = TRUE)[3]
method <<- commandArgs(trailingOnly = TRUE)[4]
param <<- commandArgs(trailingOnly = TRUE)[5]
message("Starting R script to generate results")

task.directory <- paste0(res, "-", method, "-", param)
task.name <<- paste0(tissue, "-", task.directory)
results.dir <<- paste0(output.dir, project, "/", tissue, "/", task.directory, "/") #directory for saving all other output

tiss <- read.csv(paste0(results.dir, "!cells.csv"))
tiss <- tiss %>% rename(nFeature_RNA = n_genes, nCount_RNA = n_counts, percent.mt = percent_mito, percent.rb = percent_ribo, seurat_clusters = louvain_labels)
tiss$seurat_clusters <- factor(tiss$seurat_clusters - 1)

markers <- read.csv(paste0(results.dir, "!markers.csv"))
obj.markers <- data.frame("gene"=markers$feature, "avg_logFC"=markers$log_fold_change, "p_val_adj"=markers$t_qval, 
                          "cluster"=factor(markers$cluster - 1))
obj.markers <- obj.markers %>% filter(avg_logFC > 0.25)

if (!remake.plots) {
  tmp <- assignCellTypes(tiss, obj.markers, getAnnotations(tiss), record.stats = TRUE) #assign cell types
  #unpack returned object
  clusters <<- tmp$clusters
  obj.markers <<- tmp$markers
} else {
  clusters <<- read.csv(results.dir + "!clusters.csv")
}

generatePlots(tiss, task.name, clusters$cell.type, clusters$annotation) #make plots

if (!remake.plots) {
  saveResults(tiss, clusters, obj.markers, save.cells = FALSE, save.markers = FALSE, mc_specific=FALSE)
}
