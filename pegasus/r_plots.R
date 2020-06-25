data.from.pg <- TRUE

split.special <- function(x) {
  x <- as.character(x)
  return(paste0(strsplit(x, "-")[[1]][2], "_", strsplit(x, "-")[[1]][3]))
}

source("../scripts/mc_functions.R")
source("../scripts/readers.R")
source("../scripts/settings.R")
source("../scripts/local_settings.R")
project <<- commandArgs(trailingOnly = TRUE)[1]
task.id <<- as.integer(commandArgs(trailingOnly = TRUE)[2])
tissue <<- commandArgs(trailingOnly = TRUE)[3]
res1 <<- commandArgs(trailingOnly = TRUE)[4]

message("Starting R script to generate results")

tiss <- AutoReader(project, cells.filter, features.filter, tasks.per.tiss)

parse.task.id()
res <<- res1
task.directory <- paste0(res, "-", method, "-", param)
task.name <<- paste0(tissue, "-", task.directory)
results.dir <<- paste0(output.dir, project, "/", tissue, "/", task.directory, "/") #directory for saving all other output

cells <- read.csv(paste0(results.dir, "!cells.csv"))
colnames(cells) <- c("X", "Channel", "passed_qc", "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", 
                     "seurat_clusters", "pca1", "pca2", "tsne1", "tsne2", "umap1",  "umap2") 
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

generatePlots(tiss, task.name, clusters$cell.type, clusters$annotation) #make plots
# sgenerateMarkerPlots(tiss, obj.markers %>% group_by(cluster) %>% top_n(n = 9, wt = avg_logFC))

saveResults(tiss, clusters, obj.markers, save.cells = FALSE, save.markers = FALSE, mc_specific=FALSE)

