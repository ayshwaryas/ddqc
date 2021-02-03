data.from.pg <- TRUE

source("../scripts/mc_functions.R")
source("../scripts/readers.R")
source("../scripts/settings.R")
source("../scripts/local_settings.R")
project <<- "mc_tm"
tissue <<- "Spleen"
res <<- 1.4
method <<- "mad"
param <<- 2
message("Starting R script to generate results")

output.dir <- "~/Documents/primes_storage/output_pg/"
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

clusters <<- read.csv(paste0(results.dir, "!clusters.csv"))

generatePlots(tiss, task.name, clusters$cell.type, clusters$annotation) #make plots

