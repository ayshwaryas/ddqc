data.from.pg <- TRUE

split.special <- function(x) {
  x <- as.character(x)
  return(paste0(strsplit(x, "-")[[1]][2], "_", strsplit(x, "-")[[1]][3]))
}

source("../scripts/mc_functions.R")
source("../scripts/fc_plots.R")
source("../scripts/readers.R")
source("../scripts/settings.R")
source("../scripts/local_settings.R")
project <<- commandArgs(trailingOnly = TRUE)[1]
task.id <<- as.integer(commandArgs(trailingOnly = TRUE)[2])
tissue <<- commandArgs(trailingOnly = TRUE)[3]

message("Starting R script to generate results")

tiss <- AutoReader(project, cells.filter, features.filter, tasks.per.tiss)

tasks.per.res <- tasks.per.tiss #how many different methods per one resolution
res <<- commandArgs(trailingOnly = TRUE)[4]
metric <<- "all" #switch(task.id %% tasks.per.res + 1, "counts", "genes", "mito", "ribo", "all")
results.dir <<- paste0(output.dir, project, "/", tissue, "/filtered_cells_plots/", metric, "/") #directory for saving all other output

dir.create(paste0(results.dir, "additional_plots/"), showWarnings=FALSE)

cells <- read.csv(paste0(results.dir, "!cells.csv"))
colnames(cells) <- c("X", "Channel", "color", "passed_qc", "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "seurat_clusters", "pca1", "pca2", "tsne1", "tsne2", "umap1",  "umap2") 
markers <- read.csv(paste0(results.dir, "!markers.csv"))

all.cells <- colnames(tiss$RNA)
retained <- cells$X
pass <- list()
pass[all.cells] <- FALSE
pass[intersect(all.cells, sapply(retained, split.special))] <- TRUE
tiss[["pass"]] <- factor(as.character(pass))
tiss <<- tiss[,pass == TRUE]

rownames(cells) <- sapply(cells$X, split.special)
tiss$seurat_clusters <- factor(cells[colnames(tiss), ]$seurat_clusters - 1)
Idents(tiss) <- tiss$seurat_clusters
tiss[["color"]] <- cells[colnames(tiss), ]$color
tiss[["tsne1"]] <- cells[colnames(tiss), ]$tsne1
tiss[["tsne2"]] <- cells[colnames(tiss), ]$tsne2
tiss[["umap1"]] <- cells[colnames(tiss), ]$umap1
tiss[["umap2"]] <- cells[colnames(tiss), ]$umap2

obj.markers <- data.frame("gene"=markers$feature, "avg_logFC"=markers$log_fold_change, "p_val_adj"=markers$t_qval, 
                          "cluster"=factor(markers$cluster - 1))
obj.markers <- obj.markers %>% filter(avg_logFC > 0.25)

tmp <- assignCellTypes(tiss, obj.markers, getAnnotations(tiss), record.stats = TRUE) #assign cell types
#unpack returned object
clusters <<- tmp$clusters
obj.markers <<- tmp$markers

generateFCPlots(tiss, clusters)
#generatePlots(tiss, "", clusters$cell.type, clusters$annotation, sig.plots = FALSE) #make plots

saveResults(tiss, clusters, obj.markers, save.cells = FALSE, save.markers = FALSE, mc_specific=FALSE)

