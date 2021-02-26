data.from.pg <- TRUE

source("../scripts/mc_functions.R")
source("../scripts/fc_plots.R")
source("../scripts/readers.R")
source("../scripts/settings.R")
source("../scripts/local_settings.R")
project <<- "mc_tm"
tissue <<- "Kidney" 
res <<- 1.4
metric <<- "no_outlier"
message("Starting R script to generate results")

output.dir <- "~/Documents/primes_storage/output_pg/"
results.dir <<- paste0(output.dir, project, "/", tissue, "/filtered_cells_plots/", metric, "/") #directory for saving all other output
dir.create(paste0(results.dir, "additional_plots/"), showWarnings=FALSE)

tiss <- read.csv(paste0(results.dir, "!cells.csv"))
tiss <- tiss %>% rename(nFeature_RNA = n_genes, nCount_RNA = n_counts, percent.mt = percent_mito, percent.rb = percent_ribo, seurat_clusters = louvain_labels)
tiss$seurat_clusters <- factor(tiss$seurat_clusters - 1)

clusters <<- read.csv(paste0(results.dir, "!clusters.csv"))

generateFCPlots(tiss, clusters)
