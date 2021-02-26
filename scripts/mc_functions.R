tasks.per.tiss <<- 4 #How many different res/methods per one tissue


#plots
ggsave1 <- function(filename, plot, n.clusters=30, type="h") {
  if (type == "h") {
    height = 10
    width = 14 / 30 * max(n.clusters, 30)
  }
  if (type == "v") {
    height = 10 / 30 * max(n.clusters, 30)
    width = 14
  }
  if (type == "u") {
    height = 10
    width = 10 + 2 * ceiling(n.clusters / 13)
  }
  no_bkg <- theme(axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank()) 
  ggsave(filename = filename, plot = plot + no_bkg, width = width, height = height) #saves plot with custom dimensions 
}

GetDimPlotPoints <- function(obj, reduction, metric.name) { #extracts UMAP/TSNE points for DimPlot
  cells <- colnames(obj)
  if (reduction == "umap") {
    if (!exists("data.from.pg")) {
      data <- as.data.frame(Embeddings(obj$umap)[cells, c(1, 2)])
    }
    else {
      data <- data.frame(UMAP_1 = obj$umap1, UMAP_2 = obj$umap2)
    }
    data <- data.frame(axis1 = data$UMAP_1, axis2 = data$UMAP_2, color = obj[[metric.name]], cluster = obj$seurat_clusters)
  }
  if (reduction == "tsne") {
    if (!exists("data.from.pg")) {
      data <- data <- as.data.frame(Embeddings(obj$tsne)[cells, c(1, 2)])
    }
    else {
      data <- data.frame(tSNE_1 = obj$tsne1, tSNE_2 = obj$tsne2)
    }
    data <- data.frame(axis1 = data$tSNE_1, axis2 = data$tSNE_2, color = obj[[metric.name]], cluster = obj$seurat_clusters)
  }
  return(data)
}

DimPlotContinuous <- function(obj, metric.name, lbls, name, reduction, log2=FALSE) { #DimPlot with continious colors by metric
  name <- paste0(name, "_", reduction)
  data <- GetDimPlotPoints(obj, reduction, metric.name)
  t <- theme(axis.text.x = element_text(size=15), axis.title.x = element_text(size=15), axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), plot.title = element_text(size = 20, face = "bold"), legend.title = element_text(size = 15), legend.text = element_text(size = 10)) 
  cols <- scale_colour_gradientn(colours=rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")))
  if (log2) {
    plot <- ggplot(data, aes(x=axis1, y=axis2, color=log2(color), fill = cluster)) + geom_point(size = 1) + labs(color=paste0("log2(", metric.name, ")")) + t + cols + scale_fill_discrete(labels = lbls) + ggtitle(name) + labs(x = "UMAP1", y = "UMAP2")
  } else {
    plot <- ggplot(data, aes(x=axis1, y=axis2, color=color, fill = cluster)) + geom_point(size = 1) + labs(color=metric.name) + t + cols + scale_fill_discrete(labels = lbls) + ggtitle(name) + labs(x = "UMAP1", y = "UMAP2")
  }
  for (cl in levels(obj$seurat_clusters)) { #add cluster labels
    cluster.red <- subset(data, cluster == cl)
    plot <- plot + annotate("text", x = mean(cluster.red$axis1), y = mean(cluster.red$axis2), label = cl, size = 7, fontface=2)
  } 
  return(plot)
}

DimPlotCluster <- function(obj, lbls, name, reduction) { #DimPlot colored by cluster
  name <- paste0(name, "_", reduction)
  data <- GetDimPlotPoints(obj, reduction, "seurat_clusters")
  t <- theme(axis.text.x = element_text(size=15), axis.title.x = element_text(size=15), axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), plot.title = element_text(size = 20, face = "bold"))
  plot <- ggplot(data, aes(x=axis1, y=axis2, color="Red")) + geom_point(size = 1) + ggtitle(name) + t + labs(x = "UMAP1", y = "UMAP2")
  for (cl in levels(obj$seurat_clusters)) { #add cluster labels
    cluster.red <- subset(data, cluster == cl)
    plot <- plot + annotate("text", x = mean(cluster.red$axis1), y = mean(cluster.red$axis2), cl, size = 7, fontface=2)
  }
  return(plot)
}

generatePlotsByMetric <- function(obj, name, lbls, metric.name.seurat, metric.name, name.suffix, add.dir="", save.log2=TRUE) { #different plots for QC metrics
  #themes and axis labels for plots
  t <- scale_x_discrete(labels=lbls)
  t1 <- scale_y_discrete(labels=lbls)
  t2 <- theme(axis.text.x = element_text(angle = 45, size=15, hjust=1, face="bold"), axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), legend.position="none", axis.title.x=element_blank(), plot.title = element_text(size = 20, face = "bold"))
  t3 <- ggtitle(name)
  t4 <- theme(legend.position="none")
  t5 <- facet_wrap(. ~ clusters, ncol=5, labeller = as_labeller(lbls))
  t6 <- stat_summary(fun.y=mean, geom="point", shape=23, fill="blue", size=3)
  t6h <- stat_summary(fun.x=mean, geom="point", shape=23, fill="blue", size=3)
  t7 <- theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15, face="bold"),  axis.title.x = element_text(size=15), legend.position="none", axis.title.y=element_blank(), plot.title = element_text(size = 20, face = "bold"))
  l1 <- labs(y=metric.name)
  l2 <- labs(y=paste0("log2(", metric.name, ")"))
  l3 <- labs(x=metric.name)
  l4 <- labs(x=paste0("log2(", metric.name, ")"))
  l5 <- labs(y=paste0("Average ", metric.name))
  l6 <- labs(y=paste0("Average log2(", metric.name, ")"))
  if (metric.name.seurat == "percent.mt") {
    a1 <- scale_y_continuous(breaks=seq(0, 80, 5))
    a2 <- scale_x_continuous(breaks=seq(0, 80, 5))
    hl <- geom_hline(yintercept=10, color="red", size=0.5)
    vl <- geom_vline(xintercept=10, color="red", size=0.5)
  }
  else {
    if (metric.name.seurat == "nFeature_RNA") {
      hl <- geom_hline(yintercept=log2(200), color="red", size=0.5)
      vl <- geom_vline(xintercept=log2(200), color="red", size=0.5)
    } else { 
      hl <- NULL
      vl <- NULL
    }
    a1 <- NULL
    a2 <- NULL
  }
  
  n.clusters <- length(unique(obj$seurat_clusters))
  plot.cols <- scales::hue_pal(c = 100, l = 65, h.start = 0)(n.clusters)
  names(plot.cols) <- 0:(n.clusters - 1)
  c1 <- scale_fill_manual(values = plot.cols) 
  c2 <- scale_color_manual(values = plot.cols)
  
  data <- data.frame(metric=obj[[metric.name.seurat]], clusters=obj$seurat_clusters) 
  colnames(data) <- c("metric", "clusters") #rename data columns
  
  #all plots except combined density and tsne/umap are saved in to vesions: no trasformation and log2, if save.log2 is true
  
  name.prefix <- paste0(results.dir, add.dir)
  
  if (metric.name.seurat == "nCount_RNA" || metric.name.seurat == "nFeature_RNA") {
    data$clusters = with(data, reorder(clusters, -log2(metric), mean)) #order data by cluster mean
    
    #boxplot by cluster
    ggsave1(filename=paste0(name.prefix, "box_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=clusters, y=log2(metric))) + geom_boxplot(aes()) + t + t2 + t3 + t6 + c1 + c2 + l2 + hl, n.clusters = n.clusters) 
    #joyplot by cluster
    ggsave1(filename=paste0(name.prefix, "density2_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=log2(metric), y=clusters)) + geom_density_ridges(aes()) + t1 + t3 + t6h + t7 + c1 + c2 + l4 + vl, n.clusters = n.clusters, type = "v") 
    #violin plot by cluster
    ggsave1(filename=paste0(name.prefix, "violin_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=clusters, y=log2(metric))) + geom_violin(aes()) + t + t2 + t3 + t6 + c1 + c2 + l2 + hl, n.clusters = n.clusters)
    #umap
    ggsave1(filename=paste0(name.prefix, "umap_", name.suffix, ".pdf"), plot=DimPlotContinuous(obj, metric.name.seurat, lbls, name, "umap", log2 = TRUE),  n.clusters = n.clusters, type = "u")
  } else {
    data$clusters = with(data, reorder(clusters, -metric, median)) #order data by cluster mean
    
    #boxplot by cluster
    ggsave1(filename=paste0(name.prefix, "box_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=clusters, y=metric)) + geom_boxplot(aes()) + t + t2 + t3 + t6 + c1 + c2 + l1 + a1 + hl, n.clusters = n.clusters) 
    #joyplot by cluster
    ggsave1(filename=paste0(name.prefix, "density2_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=metric, y=clusters)) + geom_density_ridges(aes()) + t1 + t3 + t6h + t7 + c1 + c2 + l3 + a2 + vl, n.clusters = n.clusters, type = "v") 
    #violin plot by cluster
    ggsave1(filename=paste0(name.prefix, "violin_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=clusters, y=metric)) + geom_violin(aes()) + t + t2 + t3 + t6 + c1 + c2 + l1 + a1 + hl, n.clusters = n.clusters)
    #umap
    ggsave1(filename=paste0(name.prefix, "umap_", name.suffix, ".pdf"), plot=DimPlotContinuous(obj, metric.name.seurat, lbls, name, "umap"), n.clusters = n.clusters, type = "u")
  }
}

generatePlots <- function(obj, name, cell.types, annotations) { #main plots function
  message("Making Plots")
  lbls <- NULL #create labels in the following format: cluster #, Panglao Cell Type \n annotated Cell Type
  for (i in 1:length(cell.types)) {
    lbls <- c(lbls, paste0((i - 1), " ", cell.types[i]))#, "\n", annotations[i]))
  }
  names(lbls) <- 0:(length(lbls) - 1) #rename labels with cluster #
  
  generatePlotsByMetric(obj, name, lbls, "nCount_RNA", "Number of UMIS", "count") #nUMI plots
  generatePlotsByMetric(obj, name, lbls, "nFeature_RNA", "Number of Genes", "genes") #nGenes plots
  generatePlotsByMetric(obj, name, lbls, "percent.mt", "percent.mt", "mito") #%mito plots
  generatePlotsByMetric(obj, name, lbls, "percent.rb", "percent.rb", "ribo") #%ribo plots
  
  #cluster colored dimplots
  ggsave1(filename=paste0(results.dir, "/umap_clusters.pdf"), plot=DimPlotCluster(obj, lbls, name, "umap"), n.clusters = length(unique(obj$seurat_clusters)), type = "u")
}

generateMarkerPlots <- function(obj, top.markers) { #generate plots of top marker genes. Does not work properly
  dir.create(paste0(results.dir, "/marker_plots/"), showWarnings = FALSE)
  n.clusters <- length(unique(obj$seurat_clusters))
  for (cl in unique(top.markers$cluster)) {
    features <- subset(top.markers, cluster == cl)$gene
    ggsave1(filename=paste0(results.dir, "/marker_plots/", cl, "-violin.pdf"), plot=VlnPlot(obj, features = features, pt.size = 0.25), n.clusters=n.clusters)
    ggsave1(filename=paste0(results.dir, "/marker_plots/", cl, "-scatter.pdf"), plot=FeaturePlot(obj, features = features, pt.size = 0.25))
    ggsave1(filename=paste0(results.dir, "/marker_plots/", cl, "-dot.pdf"), plot=DotPlot(obj, features = features))  
    message(paste("Cluster", cl, "marker plots finished"))
  }
}


#clustering & finding markers 
clusterize <- function(obj, res, compute.reductions=TRUE, compute.markers=TRUE, n.pieces=50, tsne.perplexity=30) { #function that will perform standart clustering procedures. If compute.markers == TRUE will compute DE genes. If compute.reductions = TRUE will calculate TSNE and UMAP
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(x = obj)
  obj <- ScaleData(obj, features = all.genes)
  obj <- RunPCA(obj, npcs=n.pieces)
  
  obj <- FindNeighbors(obj, dims = 1:n.pieces)
  obj <- FindClusters(obj, resolution = res)
  
  if (compute.reductions) {
    obj <- RunUMAP(obj, dims = 1:n.pieces)
    obj <- RunTSNE(obj, dims = 1:n.pieces, check_duplicates = FALSE, perplexity = tsne.perplexity)
  }
  
  if (compute.markers) {
    obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    Misc(obj, slot = "markers") <- obj.markers #store markers in seurat object
    return(list("obj "= obj, "markers" =  obj.markers))
  }
  else {
    return(obj)
  }
}


#annotations & cell types assingment
formatMarkers <- function(lst, mx=50) {
  st <- ""
  for (marker in unlist(lst)[1:min(length(unlist(lst)), mx)]) {
    st <- paste0(st, "; ", marker)
  }
  return(substring(st, 3))
}

getAnnotations <- function(obj) { #calculate most common annotations in each cluster 
  message("Assigning Annotations to Clusters")
  cluster.labels <- NULL
  cluster.labels2 <- NULL
  percents1 <- NULL
  percents2 <- NULL
  for (cl in unique(obj$seurat_clusters)) {
    cluster <- subset(obj, seurat_clusters == cl)
    #calculate first and second most frequent annotation and percentage of cells that have them
    first.frequent <- sort(table(cluster$annotations),decreasing=TRUE)[1]
    cluster.labels <- c(cluster.labels, names(first.frequent)[1])
    second.frequent <- sort(table(cluster$annotations),decreasing=TRUE)[2]
    if (is.na(second.frequent)){ #if no second annotation found
      cluster.labels2 <- c(cluster.labels2, "")
      percents2 <- c(percents2, 0)
    }
    else {
      cluster.labels2 <- c(cluster.labels2, names(second.frequent)[1])
      percents2 <- c(percents2, second.frequent / length(cluster$seurat_clusters))
    }
    percents1 <- c(percents1, first.frequent / length(cluster$seurat_clusters))
  }
  return(list("a1" = cluster.labels, "a2" = cluster.labels2, "p1" = percents1, "p2" = percents2))
}

assignCellTypes <- function(obj, markers, annotations, record.stats=TRUE, min.pval=0.05) { #function that assigns cell types based on marker genes using gene to cell.type dictionary
  message("Assigning Cell Types")
  genes <- read.csv(paste0(data.dir, "markers.tsv"), sep="\t") #read cell type markers
  genes <- data.frame(genes)
  genes <- select(genes, official.gene.symbol, cell.type, organ)
  
  markers <- markers %>% group_by(cluster) %>% filter(avg_logFC > 0, p_val_adj < min.pval) #leave only markers that have positive fold change and pval < 0.05
  markers <- markers[order(markers$cluster, -markers$avg_logFC),]
  
  cell.types <- NULL
  cell.types2 <- NULL
  cluster.ids <- NULL
  scores <- NULL
  scores2 <- NULL
  score.genes <- NULL
  
  if (record.stats) { #additional cluster stats
    cells <- NULL
    genes.mean <- NULL
    genes.median <- NULL
    mito.mean <- NULL
    mito.median <- NULL
    
    formatted.markers <- NULL #list of markers separated with ; for csv recording
  }
  
  for (cl in levels(markers$cluster)) {
    cluster.markers <- subset(markers, cluster == cl) #subset cells and markers for this cluster
    obj.cluster <- subset(obj, seurat_clusters == cl)
    cell.type <- list()
    sg <- list() #list of significant gebes
    for (i in 1:nrow(cluster.markers)) { #cycle though all cluster markers
      marker <- cluster.markers[i,]
      gene <- toupper(as.character(marker$gene))
      gene.info <- subset(genes, official.gene.symbol == gene)
      if (nrow(gene.info) > 0) {
        for (j in 1:nrow(gene.info)) { #cycle through all cell types that marker indicates for
          info <- gene.info[j,]
          type <- as.character(info$cell.type) 
          if (is.null(cell.type[[type]])) { #if this cell type was not encountored before create new list element
            cell.type[[type]] = 0
            sg[[type]] =  NULL
          }
          cell.type[[type]] = cell.type[[type]] + marker$avg_logFC #add the logFC to the score
          sg[[type]] =  c(sg[[type]], gene) #add the gene to the list of significant genes for this celltype
        }
      }
    }
    
    cluster.ids <- c(cluster.ids, cl) 
    
    if (length(cell.type) == 0) { #if no cell types were found
      cell.types <- c(cell.types, "Unknown")
      cell.types2 <- c(cell.types2, "")
      scores <- c(scores, 0)
      scores2 <- c(scores2, 0)
      score.genes <- c(score.genes, "")
    } else {
      cell.type <- sort(unlist(cell.type, use.names = TRUE), decreasing = TRUE)
      score.genes <- c(score.genes, formatMarkers(sg[[names(cell.type[1])]]))
      if (length(sg[[names(cell.type[1])]]) < 3) { #if most common celltype has less than 3 significant genes switch it to unkwown
        cell.types <- c(cell.types, "Unknown")
        if (length(cell.type) >= 1) {
          cell.types2 <- c(cell.types2, names(cell.type[1]))
        } else {
          cell.types2 <- c(cell.types2, "")
        }
        scores <- c(scores, 0)
        scores2 <- c(scores2, 0)
      } else {
        if (length(cell.type) > 1) { #if more than one cell type was found
          cell.types <- c(cell.types, names(cell.type[1]))
          cell.types2 <- c(cell.types2, names(cell.type[2]))
          scores <- c(scores, cell.type[[1]]) #absoute score
          scores2 <- c(scores2, (cell.type[[1]] - cell.type[[2]]) / cell.type[[1]]) #score accounting for the difference btw first and second cell type
        }
        if (length(cell.type) == 1) { #if only one cell type was found
          cell.types <- c(cell.types, names(cell.type[1]))
          cell.types2 <- c(cell.types2, "")
          scores <- c(scores, cell.type[[1]])
          scores2 <- c(scores2, 1)
        }
      }
    }
    
    if (record.stats) { #recorbd additional statistics about cluster
      cells <- c(cells, length(obj.cluster$seurat_clusters))
      genes.mean <- c(genes.mean, mean(obj.cluster$nFeature_RNA))
      genes.median <- c(genes.median, median(obj.cluster$nFeature_RNA))
      mito.mean <- c(mito.mean, mean(obj.cluster$percent.mt))
      mito.median <- c(mito.median, median(obj.cluster$percent.mt))
      
      formatted.markers <- c(formatted.markers, formatMarkers(cluster.markers$gene)) #format markers for csv recording by separating them with ;
    }
  }
  if (record.stats) {
    clusters <- tibble("cluster" = cluster.ids, "annotation" = annotations[["a1"]], "annotation2" = annotations[["a2"]], "%annotation1" = round(annotations[["p1"]], 3), 
                  "%annotation2" = round(annotations[["p2"]], 3), "cell.type" = cell.types, "cell.type2" = cell.types2, "sum.avg_logFC" = round(scores, 3), 
                  "delta.score" = round(scores2, 3), "cells" = cells, "genes.mean" = round(genes.mean, 3), "genes.median" = round(genes.median, 3), 
                  "mito.mean" = round(mito.mean, 3), "mito.median" = round(mito.median, 3), "markers" = formatted.markers, "score.genes" = score.genes)
    return(list("clusters" = clusters, "markers" = markers))
  }
  else { #return only cell types and annotations
    res <- tibble("cluster" = cluster.ids, "annotation" = annotations[["a1"]], "cell.type" = cell.types)
    return(res)
  }
}


saveResults <- function(obj, clusters, obj.markers, save.cells=TRUE, save.markers=TRUE, mc_specific=FALSE) {
  message("Writing Results")

  write.csv(clusters, file=paste0(results.dir, "/!clusters.csv"))
  if (save.markers) {
    write.csv(obj.markers, file=paste0(results.dir, "/!markers.csv"))
  }
  
  all.cells <- colnames(obj)
  
  if (!exists("data.from.pg")) {
    umap <- as.data.frame(Embeddings(obj$umap)[all.cells, c(1, 2)])
    tsne <- as.data.frame(Embeddings(obj$tsne)[all.cells, c(1, 2)])
  }
  else {
    umap <- data.frame(UMAP_1 = obj$umap1, UMAP_2 = obj$umap2)
    tsne <- data.frame(tSNE_1 = obj$tsne1, tSNE_2 = obj$tsne2)
  } 
  
  if (save.cells) {
    cells <- data.frame("cell"=all.cells, "cluster"=obj$seurat_clusters, "cell.type"=clusters$cell.type[obj$seurat_clusters], 
                        "cell.type2"=clusters$cell.type2[obj$seurat_clusters], 
                        "annotation" = clusters$annotation[obj$seurat_clusters], "annotation2" = clusters$annotation2[obj$seurat_clusters],
                        "nCount_RNA" = obj$nCount_RNA, "nFeature_RNA" = obj$nFeature_RNA, "percent.mt" = obj$percent.mt, "percent.rb" = obj$percent.rb, 
                        "UMAP_1"=umap$UMAP_1, "UMAP_2"=umap$UMAP_2,
                       "TSNE_1"=tsne$tSNE_1, "TSNE_2"=tsne$tSNE_2)
    write.csv(cells, file=paste0(results.dir, "/!cells.csv"))
  }
  
  if (mc_specific) {
    if (res == 1 && method == "none") { 
      message("Recording Additional Stats")
      #record all cells with QC stats into csv for summary plots
      t <- tibble("tissue" = tissue, "cell" = colnames(obj), "cluster" = obj$seurat_clusters, "nCount_RNA" = obj$nCount_RNA, "nFeature_RNA" = obj$nFeature_RNA, "percent.mt" = obj$percent.mt, "percent.rb" = obj$percent.rb)
      write.table(t, paste0(results.dir, "../../stats_summary.csv"), sep=",", append=TRUE, col.names=FALSE)
    }
    
    if (check.save()) {
      saveRDS(obj, file=paste0(robjs.dir, task.name, ".rds"))
    }
  }
}


check.save <- function() {
  if (save.res.1 && res == 1 && method == "none") {
    return(TRUE)
  }
  return(tryCatch({
      files.to.save <- read.table("files_to_save.txt")
      for (file in files.to.save$V1) {
        if (tolower(paste(project, tissue, res, method, param, sep="-")) == tolower(file)) {
          return(TRUE)
        }
      }
      return(FALSE)
  }, error=function(e) FALSE))
}


parse.task.id <- function() {
  tasks.per.res <- tasks.per.tiss #how many different methods per one resolution
  res <<- 1 #0.5 * (1 + ((task.id %% tasks.per.tiss) %/% tasks.per.res)) #clustering resolution
  method <<- switch(task.id %% tasks.per.res + 1, "none", "cutoff", "outlier", "mad") #filtering method
  param <<- switch(task.id %% tasks.per.res + 1, 0, 10, 0, 2) #filtering parameter
}

#main function
MCMain <- function() {
  parse.task.id()
  info.msg <- paste0("task.id:", task.id, " - tissue:", tissue, " res:", res, " mehtod:", method, " param:", param, 
                     " project:", project, " do.counts:", do.counts, " do.genes:", do.genes, " do.mito:", do.mito,
                     " do.ribo:", do.ribo)
  message(paste0("Starting ", info.msg))
  
  task.directory <- paste0(res, "-", method, "-", param)
  task.name <<- paste0(tissue, "-", task.directory)
  robjs.dir <<- paste0(output.dir, "robjs/", project, "/") #directory for saving R objects
  results.dir <<- paste0(output.dir, project, "/", tissue, "/", task.directory, "/") #directory for saving all other output
  
  message("Creating Output Direcrories")
  dir.create(paste0(output.dir, project, "/"), showWarnings=FALSE)
  dir.create(paste0(output.dir, project, "/", tissue), showWarnings=FALSE)
  dir.create(paste0(results.dir), showWarnings=FALSE)
  dir.create(paste0(results.dir, "/signatures"), showWarnings=FALSE)
  dir.create(paste0(output.dir, "robjs/"), showWarnings=FALSE)
  dir.create(robjs.dir, showWarnings=FALSE)
  
  tiss <<- filterCells(tiss, method, param, do.counts, do.genes, do.mito, do.ribo)
  
  tmp <- clusterize(tiss, res, compute.reductions = TRUE, compute.markers = TRUE) #cluster cells
  #unpack returned object
  tiss <<- tmp$obj
  obj.markers <<- tmp$markers
  
  tmp <- assignCellTypes(tiss, obj.markers, getAnnotations(tiss)) #assign cell types
  #unpack returned object
  clusters <<- tmp$clusters
  obj.markers <<- tmp$markers
  
  generatePlots(tiss, task.name, clusters$cell.type, clusters$annotation) #make plots
  # sgenerateMarkerPlots(tiss, obj.markers %>% group_by(cluster) %>% top_n(n = 9, wt = avg_logFC))
  
  saveResults(tiss, clusters, obj.markers, mc_specific=TRUE)
  
  message(paste0("Finished task.id:", info.msg))
}
