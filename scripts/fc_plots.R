tasks.per.tiss <<- 1 #How many different res/methods per one tissue

readFilterCsvMethod <- function(method) {
  filtered.cells <- NULL
  for (metric in c("counts")) {#c("counts", "genes", "mito", "ribo")) {
    filtered.cells <- c(filtered.cells, tryCatch({
      as.character(read.csv(paste0(source.dir, res, "-", method, "/!filtered_", metric, ".csv"))[["cell"]])}, 
                                                 error = function(e) {warning(paste(method, metric, "filtered cells not found"))}))
  }
  return(filtered.cells)
}


#function which reads data about filtered cells and categorizes them
readFilterCsv <- function(obj) {
  message("Reading Filtered Cells")
  all.cells = colnames(obj$RNA)
  
  cutoff10 <- setdiff(all.cells, readFilterCsvMethod("cutoff-10"))
  zscore2 <- setdiff(all.cells, readFilterCsvMethod("z_score-2"))
  mad <- setdiff(all.cells, readFilterCsvMethod("mad-2"))
  
  #categorize cells
  color <- list()
  color[all.cells] <- "Did not pass"
  color[mad] <- "MAD2 only"
  color[cutoff10] <- "C10 only"
  color[zscore2] <- "ZSC2 only"
  color[intersect(mad, cutoff10)] <- "MAD2 and C10"
  color[intersect(mad, zscore2)] <- "MAD2 and ZSC2"
  color[intersect(zscore2, cutoff10)] <- "ZSC2 and C10"
  color[intersect(mad, intersect(cutoff10, zscore2))] <- "All"
  obj[["color"]] <- factor(as.character(color))
  
  obj <- obj[,color != "Did not pass"] #filter out cells that did not pass any qc method
  return(obj)
}


#function that makes filtered cells plots
generateFCPlots <- function(obj, clusters) {
  message("Making FC Plots")
  lbls <- NULL #create labels in the following format: cluster #, Panglao Cell Type \n annotated Cell Type
  for (i in 1:length(clusters$cell.type)) {
    lbls <- c(lbls, paste0((i - 1), " ", clusters$cell.type[i], "\n", clusters$annotation[i]))
  }
  names(lbls) <- 0:(length(lbls) - 1)  #rename labels with cluster #
  
  #extract UMAP coordinates for plotting
  cells <- colnames(obj)
  data <- as.data.frame(Embeddings(obj$umap)[cells, c(1, 2)])
  data <- data.frame(UMAP_1 = data$UMAP_1, UMAP_2 = data$UMAP_2, cluster = obj$seurat_clusters, color=obj$color, annotation=obj$annotations)
  
  t1 <- theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  t2 <- guides(colour = guide_legend(override.aes = list(size=2)))
  plot.cols <- c("C10 only" = "#D55E00", "MAD2 and C10" = "#E69F00", "ZSC2 and C10" = "#CC79A7", "MAD2 only" = "#56B4E9",
                 "ZSC2 only" = "#0072B2", "MAD2 and ZSC2" = "#009E73", "All" = "#999999" , "Did not pass" = "#FFFFFF") #to keep consistent plot colors
  
  plot1 <- ggplot(data, aes(x=UMAP_1, y=UMAP_2, color=color)) + geom_point(size = 0.5) + t1 + t2 + scale_fill_manual(values = plot.cols) + scale_color_manual(values = plot.cols) #umap colored by filtering category
  plot2 <- ggplot(data, aes(x=UMAP_1, y=UMAP_2, color=cluster)) + geom_point(size = 0.5) + t1 + t2 #umap colored by cluster
  for (cl in levels(obj$seurat_clusters)) { #add labels to plots
    cluster.umap <- subset(data, cluster == cl)
    plot1 <- plot1 + annotate("text", x = mean(cluster.umap$UMAP_1), y = mean(cluster.umap$UMAP_2), label = lbls[(as.numeric(cl) + 1)], size = 3, fontface=2)
    plot2 <- plot2 + annotate("text", x = mean(cluster.umap$UMAP_1), y = mean(cluster.umap$UMAP_2), label = lbls[(as.numeric(cl) + 1)], size = 3, fontface=2)
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
    scale_fill_manual(values = plot.cols) + scale_x_discrete(labels=lbls) + theme(axis.text.x = element_text(angle = 45, size=10, hjust=1, face="bold"),  axis.title.x=element_blank())
  
  #write plots
  n.clusters <- length(unique(obj$seurat_clusters))
  ggsave1(filename = paste0(results.dir, res, "-", mito.cutoff, "-filterplot.pdf"), plot=plot1)
  ggsave1(filename = paste0(results.dir, res, "-", mito.cutoff, "-clusterplot.pdf"), plot=plot2)
  ggsave1(filename = paste0(results.dir, res, "-", mito.cutoff, "-barplot.pdf"), plot=plot3, n.clusters=n.clusters)
}


#main function
FCPlotsMain <- function() {
  tasks.per.res <- tasks.per.tiss #how many different methods per one resolution
  res <<- 1 #0.5 * (1 + (task.id %% tasks.per.tiss) %/% tasks.per.res) #clustering resolution
  mito.cutoff <<- 80 #switch(task.id %% tasks.per.res + 1, 100, 80) #cutoff for percent mito
  tiss <<- subset(tiss, percent.mt <= mito.cutoff) #subset mito genes
  message(paste("Starting task.id:", task.id, "- tissue:", tissue, "res:", res, "mito.cutoff", mito.cutoff, "project:", project))
  
  source.dir <<- paste0(source.dir, project, "/", tissue, "/") #directory where csv with filtered cells are located
  results.dir <<- paste0(output.dir, project, "/", tissue, "/filtered_cells_plots/") #directory where output will go
  dir.create(results.dir, showWarnings=FALSE)
  
  tiss <<- readFilterCsv(tiss) #add cell categories
  
  tmp <- clusterize(tiss, res, compute.reductions=TRUE, compute.markers = TRUE) #cluster cells
  #unpack returned object
  tiss <<- tmp$obj
  obj.markers <<- tmp$markers
  clusters <<- assignCellTypes(tiss, obj.markers, getAnnotations(tiss))
  
  generateFCPlots(tiss, clusters) #make plots
  
  saveResults(tiss, clusters, obj.markers)
  message(paste("Finished task.id:", task.id, "- tissue:", tissue, "res:", res, "mito.cutoff", mito.cutoff, "project:", project))
}

