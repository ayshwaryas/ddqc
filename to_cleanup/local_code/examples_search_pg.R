#results.dir <- "~/Documents/primes_storage/output_pg/"
#results.dir <- "/Volumes/scqc/output_pg/"
results.dir <- "~/Downloads/"
write("project,tissue,cluster,annotation,cell.type,n.cells,percent.unique,qc.criteria,genes.median,mito.median,ribo.median,markers", "result.csv")

for (project in dir(results.dir,  pattern="mc", recursive=FALSE)) {
  for (tissue in dir(paste0(results.dir, project, "/"), recursive=FALSE)) {
    if (tissue == "stats_summary.csv") {
      next
    }
    source.dir <- paste0(results.dir, project, "/", tissue, "/filtered_cells_plots/no_outlier/")
    
    if (!file.exists(paste0(source.dir, "!cells.csv")) || !file.exists(paste0(source.dir, "!clusters.csv"))) {
      print(tissue)
      next
    }
    
    cells <- read.csv(paste0(source.dir, "!cells.csv"))
    cells$louvain_labels <- cells$louvain_labels - 1
    
    rownames(cells) <- cells$barcodekey
    clusters <- read.csv(paste0(source.dir, "!clusters.csv"))
    rownames(clusters) <- clusters$cluster
    
    for (cl in rownames(clusters)) {
      cluster.cells <- subset(cells, louvain_labels == cl)
      ribo.median <- round(median(cluster.cells$percent_ribo), 3)
      mad.only.cells <- rownames(subset(cluster.cells, color %in% c("MAD2 only")))
      # print(length(mad.only.cells) / length(cluster.cells$barcodekey))     
      if (length(cluster.cells$barcodekey) >= 30 && length(mad.only.cells) / length(cluster.cells$barcodekey) >= 0.9) {
        mrk <- strsplit(as.character(clusters[cl,]$markers), split=";")[[1]]
        qc.criteria.line <- ""
        if (clusters[cl,]$genes.median <= 200) {
          qc.criteria.line <- "n_genes"
        }
        if (clusters[cl,]$mito.median >= 10) {
          qc.criteria.line <- "percent_mito"
        }
        if (clusters[cl,]$genes.median <= 200 && clusters[cl,]$mito.median >= 10) {
          qc.criteria.line <- "n_genes and percent_mito"
        }
        write(paste(project, tissue, cl, as.character(clusters[cl,]$annotation), as.character(clusters[cl,]$cell.type), length(cluster.cells$barcodekey), round(length(mad.only.cells) / length(cluster.cells$barcodekey), 4), qc.criteria.line, as.character(clusters[cl,]$genes.median), as.character(clusters[cl,]$mito.median), ribo.median, paste(mrk[1:min(25, length(mrk))], collapse=";", sep=""),  sep=","), file = "result.csv", append = TRUE)
      }
    }
  }
}
