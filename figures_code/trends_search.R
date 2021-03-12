results.dir <- "/ahg/regevdata/projects/scqc/output_pg/"
fout_ng_lower <- "trends_n_genes_lower.csv"
fout_ng_upper <- "trends_n_genes_upper.csv"
fout_mito <- "trends_mito_upper.csv"
fout_ribo <- "trends_ribo_upper.csv"

n_genes_lower_bound <- 200
n_genes_upper_bound <- 2000
mito_upper_bound <- 20
ribo_upper_bound <- 30

HEADER <- "project,tissue,cluster,annotation,cell_type,n_cells,genes_mean,genes_median,mito_mean,mito_median,ribo_mean,ribo_median,markers"
write(HEADER, fout_ng_lower)
write(HEADER, fout_ng_upper)
write(HEADER, fout_mito)
write(HEADER, fout_ribo)

for (project in dir(results.dir, recursive=FALSE)) {
  for (tissue in dir(paste0(results.dir, project, "/"), recursive=FALSE)) {
    source.dir <- paste0(results.dir, project, "/", tissue, "/1.4-mad-2/")
    
    if (!file.exists(paste0(source.dir, "!cells.csv")) || !file.exists(paste0(source.dir, "!clusters.csv"))) {
      print(paste0(project, "-", tissue))
      next
    }
    
    cells <- read.csv(paste0(source.dir, "!cells.csv"))
    rownames(cells) <- cells$barcodekey
    clusters <- read.csv(paste0(source.dir, "!clusters.csv"))
    rownames(clusters) <- clusters$cluster
    
    for (cl in rownames(clusters)) {
      cluster.cells <- subset(cells, louvain_labels == cl)
      
      if (length(cluster.cells$barcodekey) >= 0) {
        mrk <- strsplit(as.character(clusters[cl,]$markers), split=";")[[1]]
        info <- paste(project, tissue, cl, as.character(clusters[cl,]$annotation), as.character(clusters[cl,]$cell_type), length(cluster.cells$barcodekey), as.character(clusters[cl,]$genes_mean), as.character(clusters[cl,]$genes_median), as.character(clusters[cl,]$mito_mean), as.character(clusters[cl,]$mito_median), as.character(clusters[cl,]$ribo_mean), as.character(clusters[cl,]$ribo_median), paste(mrk[1:min(25, length(mrk))], collapse=";", sep=""),  sep=",")
        if (clusters[cl,]$genes_median <= n_genes_lower_bound) {
          write(info, file=fout_ng_lower, append = TRUE)
        }
        if (clusters[cl,]$genes_median >= n_genes_upper_bound) {
          write(info, file=fout_ng_upper, append = TRUE)
        }
        if (clusters[cl,]$mito_median >= mito_upper_bound) {
          write(info, file=fout_mito, append = TRUE)
        }
        if (clusters[cl,]$ribo_median >= ribo_upper_bound) {
          write(info, file=fout_ribo, append = TRUE)
        }
      }
    }
  }
}
