data.from.pg <- TRUE

source("../scripts/mc_functions.R")
source("../scripts/settings.R")
source("../scripts/local_settings.R")
project <<- commandArgs(trailingOnly = TRUE)[1]
task.id <<- as.integer(commandArgs(trailingOnly = TRUE)[2])
tissue <<- commandArgs(trailingOnly = TRUE)[3]

message("Starting R script to generate results")

parse.task.id()
task.directory <- paste0(res, "-", method, "-", param)
task.name <<- paste0(tissue, "-", task.directory)
results.dir <<- paste0("~/Documents/primes_storage/output_pg/", project, "/", tissue, "/", task.directory, "/") #directory for saving all other output

cells <- read.csv(paste0(results.dir, "!cells.csv"))
colnames(cells) <- c("X", "Channel", "passed_qc", "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", 
                     "seurat_clusters", "pca1", "pca2", "tsne1", "tsne2", "umap1",  "umap2") 
markers <- read.csv(paste0(results.dir, "!markers.csv"))

cells$seurat_clusters <- factor(cells$seurat_clusters - 1)
cells["idents"] <- cells$seurat_clusters
cells["annotations"] <- "Unknown"

obj.markers <- data.frame("gene"=markers$feature, "avg_logFC"=markers$percentage_fold_change, "p_val_adj"=markers$t_pval, 
                          "cluster"=factor(markers$cluster - 1))
obj.markers <- obj.markers %>% filter(avg_logFC > 0.25)

tmp <- assignCellTypes(cells, obj.markers, getAnnotations(cells), record.stats = TRUE) #assign cell types
#unpack returned object
clusters <<- tmp$clusters
sm <- tmp$sm
obj.markers <<- tmp$markers

generatePlots(cells, task.name, clusters$cell.type, clusters$annotation, sig.plots = FALSE) #make plots
# sgenerateMarkerPlots(tiss, obj.markers %>% group_by(cluster) %>% top_n(n = 9, wt = avg_logFC))

saveResults(cells, clusters, obj.markers, save.cells = FALSE, save.markers = FALSE, mc_specific=FALSE, sm=sm)
