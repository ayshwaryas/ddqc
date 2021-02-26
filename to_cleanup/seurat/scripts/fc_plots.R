tasks.per.tiss <<- 1 #How many different res/methods per one tissue

readFilterCsvMethod <- function(method, metric) {
  filtered.cells <- NULL
  if (metric == "all") {
    for (metric1 in c("counts", "genes", "mito", "ribo")) {
      filtered.cells <- c(filtered.cells, tryCatch({
        as.character(read.csv(paste0(source.dir, res, "-", method, "/!filtered_", metric1, ".csv"))[["cell"]])}, 
        error = function(e) {warning(paste(method, metric1, "filtered cells not found"))}))
    }
  }
  else {
    filtered.cells <- c(filtered.cells, tryCatch({
      as.character(read.csv(paste0(source.dir, res, "-", method, "/!filtered_", metric, ".csv"))[["cell"]])}, 
                                                 error = function(e) {warning(paste(method, metric, "filtered cells not found"))}))
  }
  return(filtered.cells)
}


#function which reads data about filtered cells and categorizes them
readFilterCsv <- function(obj, metric) {
  message("Reading Filtered Cells")
  if (!exists("data.from.pg")) {
    all.cells = colnames(obj$RNA)
  }
  else {
    all.cells <- obj$cell
  }
  
  cutoff10 <- setdiff(all.cells, readFilterCsvMethod("cutoff-10", metric))
  outlier <- setdiff(all.cells, readFilterCsvMethod("outlier-0", metric))
  mad <- setdiff(all.cells, readFilterCsvMethod("mad-2", metric))
  
  #categorize cells
  color <- list()
  color[all.cells] <- "Did not pass"
  color[mad] <- "MAD2 only"
  color[cutoff10] <- "C10 only"
  color[outlier] <- "Outlier only"
  color[intersect(mad, cutoff10)] <- "MAD2 and C10"
  color[intersect(mad, outlier)] <- "MAD2 and Outlier"
  color[intersect(outlier, cutoff10)] <- "Outlier and C10"
  color[intersect(mad, intersect(cutoff10, outlier))] <- "All"
  obj[["color"]] <- factor(as.character(color))
  
  obj <- subset(obj ,color != "Did not pass") #filter out cells that did not pass any qc method
  return(obj)
}


#function that makes filtered cells plots
generateFCPlots <- function(obj, clusters) {
  message("Making FC Plots")
  lbls <- NULL #create labels in the following format: cluster #, Panglao Cell Type \n annotated Cell Type
  for (i in 1:length(clusters$cell.type)) {
    lbls <- c(lbls, paste0((i - 1), " ", clusters$cell.type[i]))#, "\n", clusters$annotation[i]))
  }
  names(lbls) <- 0:(length(lbls) - 1)  #rename labels with cluster #
  
  name <- paste0(res, "-", metric, "-", tissue)
  
  generatePlotsByMetric(obj, name, lbls, "nCount_RNA", "Number of UMIS", "count") #nUMI plots
  generatePlotsByMetric(obj, name, lbls, "nFeature_RNA", "Number of Genes", "genes") #nGenes plots
  generatePlotsByMetric(obj, name, lbls, "percent.mt", "percent.mt", "mito") #%mito plots
  generatePlotsByMetric(obj, name, lbls, "percent.rb", "percent.rb", "ribo") #%ribo plots
  
  #extract UMAP coordinates for plotting
  if (!exists("data.from.pg")) {
    cells <- colnames(obj)
    data <- as.data.frame(Embeddings(obj$umap)[cells, c(1, 2)])
  } else {
    data <- data.frame(UMAP_1 = obj$umap1, UMAP_2 = obj$umap2)
  }
  
  plot.cols <- c("Cutoff only" = "#4A3933", "MAD2 only" = "#16697A","All" = "#FFA62B" , "Neither" = "#AAAAAA") #to keep consistent plot colors
  
  
  data <- data.frame(UMAP_1 = data$UMAP_1, UMAP_2 = data$UMAP_2, cluster = obj$seurat_clusters, color=obj$color, annotation=obj$annotations)
  
  t1 <- theme(axis.text.x = element_text(size=15), axis.title.x = element_text(size=15), axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), plot.title = element_text(size = 20, face = "bold"), legend.title = element_text(size = 15), legend.text = element_text(size = 10))
  t2 <- guides(colour = guide_legend(override.aes = list(size=2)))
  t3 <- scale_color_discrete(labels = lbls)
  t4 <- scale_fill_discrete(labels = lbls)
  
  plot1 <- ggplot(data, aes(x=UMAP_1, y=UMAP_2, color=color, fill=cluster)) + geom_point(size = 1) + t1 + t2 + t3 + t4 + scale_color_manual(values = plot.cols) + ggtitle(name) #umap colored by filtering category
  plot2 <- ggplot(data, aes(x=UMAP_1, y=UMAP_2, color=cluster)) + geom_point(size = 1) + t1 + t2 + t3 + t4 + ggtitle(name) #umap colored by cluster
  for (cl in levels(obj$seurat_clusters)) { #add labels to plots
    cluster.umap <- subset(data, cluster == cl)
    plot1 <- plot1 + annotate("text", x = mean(cluster.umap$UMAP_1), y = mean(cluster.umap$UMAP_2), label = cl, size = 7, fontface=2)
    plot2 <- plot2 + annotate("text", x = mean(cluster.umap$UMAP_1), y = mean(cluster.umap$UMAP_2), label = cl, size = 7, fontface=2)
  }
  
  #create barplot which shows category distribution
  table.color <- NULL
  table.cluster <- NULL
  table.freq <- NULL
  for (cl in unique(data$cluster)) {
    data.cluster <- subset(data, cluster == cl)
    table.tmp <- as.data.frame(table(data.frame(color=factor(data.cluster$color), cluster = as.character(data.cluster$cluster))))
    table.tmp$Freq <- table.tmp$Freq / sum(table.tmp$Freq)
    table.color <- c(table.color, as.character(table.tmp$color))
    table.cluster <- c(table.cluster, as.integer(as.character(table.tmp$cluster)))
    table.freq <- c(table.freq, as.character(table.tmp$Freq))
  }
  data1 <- data.frame(color=table.color, cluster = as.factor(table.cluster), freq=as.double(as.character(table.freq)) * 100)
  
  plot3 <- ggplot(data1, aes(x=cluster, y=freq, fill=color)) + geom_bar(stat="identity") + 
    scale_fill_manual(values = plot.cols) + scale_x_discrete(labels=lbls) + theme(axis.text.x = element_text(angle = 45, size=15, hjust=1, face="bold"), legend.title = element_text(size = 15), legend.text = element_text(size = 10), axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), axis.title.x=element_blank(), plot.title = element_text(size = 20, face = "bold")) + ggtitle(name)
  
  
  #write plots
  n.clusters <- length(unique(obj$seurat_clusters))
  ggsave1(filename = paste0(results.dir, res, "-!filterplot.pdf"), plot=plot1, n.clusters = length(unique(obj$seurat_clusters)), type = "u")
  ggsave1(filename = paste0(results.dir, res, "-!clusterplot.pdf"), plot=plot2, n.clusters = length(unique(obj$seurat_clusters)), type = "u")
  ggsave1(filename = paste0(results.dir, res, "-!barplot.pdf"), plot=plot3, n.clusters=n.clusters)
}


#main function
FCPlotsMain <- function() {
  tasks.per.res <- tasks.per.tiss #how many different methods per one resolution
  res <<- 1 #0.5 * (1 + (task.id %% tasks.per.tiss) %/% tasks.per.res) #clustering resolution
  metric <<- "all" #switch(task.id %% tasks.per.res + 1, "counts", "genes", "mito", "ribo", "all")
  message(paste("Starting task.id:", task.id, "- tissue:", tissue, "res:", res, "metric:", metric, "project:", project))
  
  source.dir <<- paste0(source.dir.prefix, project, "/", tissue, "/") #directory where csv with filtered cells are located
  results.dir <<- paste0(output.dir, project, "/", tissue, "/filtered_cells_plots/", metric, "/") #directory where output will go
  dir.create(paste0(output.dir, project, "/", tissue, "/filtered_cells_plots/"), showWarnings=FALSE)
  dir.create(results.dir, showWarnings=FALSE)
  dir.create(paste0(results.dir, "additional_plots/"), showWarnings=FALSE)
  
  tiss <<- readFilterCsv(tiss, metric) #add cell categories
  
  tmp <- clusterize(tiss, res, compute.reductions=TRUE, compute.markers = TRUE) #cluster cells
  #unpack returned object
  tiss <<- tmp$obj
  obj.markers <<- tmp$markers
  tmp <- assignCellTypes(tiss, obj.markers, getAnnotations(tiss)) #assign cell types
  #unpack returned object
  clusters <<- tmp$clusters
  obj.markers <<- tmp$markers
  
  generateFCPlots(tiss, clusters)
  #generatePlots(tiss, "", clusters$cell.type, clusters$annotation, sig.plots = FALSE) #make plots
  
  saveResults(tiss, clusters, obj.markers, mc_specific=FALSE)
  message(paste("Finished task.id:", task.id, "- tissue:", tissue, "res:", res, "metric:", metric, "project:", project))
}

