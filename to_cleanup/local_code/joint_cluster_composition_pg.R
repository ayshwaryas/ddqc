results.dir <- "~/Documents/primes_storage/output_pg/"
write("project,tissue,cluster,annotation,cell.type,n.cells,cells.all,percent.all,cells.mad,percent.mad,cells.cutoff,percent.cutoff,genes.median,mito.median,ribo.median,markers", "~/Downloads/joint_cluster_composition.csv")

for (project in dir(results.dir,  pattern="mc_human", recursive=FALSE)) {
  if (project == "mc_other_10X" || project == "mc_other_10X" || project == "mc_PanglaoDB") {
    next
  }
  for (tissue in dir(paste0(results.dir, project, "/"), recursive=FALSE)) {
    if (tissue == "stats_summary.csv") {
      next
    }
    if (tissue != "adipose" && tissue != "liver") {
      next
    }
    source.dir <- paste0(results.dir, project, "/", tissue, "/filtered_cells_plots/no_outlier/")
    
    if (!file.exists(paste0(source.dir, "!cells.csv"))) {
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
      cutoff.only.cells <- rownames(subset(cluster.cells, color %in% c("C10 only")))
      all.only.cells <- rownames(subset(cluster.cells, color %in% c("All")))
     
      mrk <- strsplit(as.character(clusters[cl,]$markers), split=";")[[1]]
      write(paste(project, tissue, cl, as.character(clusters[cl,]$annotation), as.character(clusters[cl,]$cell.type), 
                  length(cluster.cells$barcodekey),
                  length(all.only.cells),
                  round(length(all.only.cells) / length(cluster.cells$barcodekey), 4) * 100,
                  length(mad.only.cells),
                  round(length(mad.only.cells) / length(cluster.cells$barcodekey), 4) * 100, 
                  length(cutoff.only.cells),
                  round(length(cutoff.only.cells) / length(cluster.cells$barcodekey), 4) * 100, 
                  as.character(clusters[cl,]$genes.median), 
                  as.character(clusters[cl,]$mito.median), ribo.median, 
                  paste(mrk[1:min(25, length(mrk))], collapse=";", sep=""),  
                  sep=","), file = "~/Downloads/joint_cluster_composition.csv", append = TRUE)
      
    }
  }
}
