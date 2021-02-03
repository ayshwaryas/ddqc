data.from.pg <- TRUE


source("../scripts/mc_functions.R")
source("../scripts/fc_plots.R")
source("../scripts/settings.R")
source("../scripts/local_settings.R")
project <<- "mc_tm"
tissue <<- "Lung"
res <<- 1.4
metric <<- "no_outlier"

results.dir <<-  paste0("/ahg/regevdata/projects/scqc/figure2_data/", project, "_", tissue, "/")
dir.create(paste0(results.dir, "additional_plots/"), showWarnings=FALSE)

tiss <- read.csv(paste0(results.dir, "!cells.csv"))
tiss <- tiss %>% rename(nFeature_RNA = n_genes, nCount_RNA = n_counts, percent.mt = percent_mito, percent.rb = percent_ribo, seurat_clusters = louvain_labels)
tiss$seurat_clusters <- factor(tiss$seurat_clusters - 1)

clusters <<- read.csv(paste0(results.dir, "!clusters.csv"))

generateFCPlots(tiss, clusters)

