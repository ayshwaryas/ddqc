tasks.per.tiss <<- 5 #How many different res/methods per one tissue


#plots
ggsave1 <- function(filename, plot, n.clusters=30) {
  n.clusters <- max(10, n.clusters)
  no_bkg <- theme(axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank()) 
  ggsave(filename = filename, plot = plot + no_bkg, width = 14 / 30 * n.clusters, height = 10) #saves plot with custom dimensions 
}

GetDimPlotPoints <- function(obj, reduction, metric.name) { #extracts UMAP/TSNE points for DimPlot
  cells <- colnames(obj)
  if (reduction == "umap") {
    data <- as.data.frame(Embeddings(obj$umap)[cells, c(1, 2)])
    data <- data.frame(axis1 = data$UMAP_1, axis2 = data$UMAP_2, color = obj[[metric.name]], cluster = obj$seurat_clusters)
  }
  if (reduction == "tsne") {
    data <- as.data.frame(Embeddings(obj$tsne)[cells, c(1, 2)])
    data <- data.frame(axis1 = data$tSNE_1, axis2 = data$tSNE_2, color = obj[[metric.name]], cluster = obj$seurat_clusters)
  }
  return(data)
}

DimPlotContinuous <- function(obj, metric.name, lbls, name, reduction) { #DimPlot with continious colors by metric
  name <- paste0(name, "_", reduction)
  data <- GetDimPlotPoints(obj, reduction, metric.name) 
  plot <- ggplot(data, aes(x=axis1, y=axis2, color=eval(parse(text=eval(metric.name))))) + geom_point(size = 0.5) + scale_colour_gradientn(colours=rev(rainbow(4))) + labs(color=metric.name) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + ggtitle(name)
  for (cl in levels(obj$seurat_clusters)) { #add cluster labels
    cluster.red <- subset(data, cluster == cl)
    plot <- plot + annotate("text", x = mean(cluster.red$axis1), y = mean(cluster.red$axis2), label = lbls[as.numeric(cl) + 1], size = 3, fontface=2)
  } 
  return(plot)
}

DimPlotCluster <- function(obj, lbls, name, reduction) { #DimPlot colored by cluster
  name <- paste0(name, "_", reduction)
  data <- GetDimPlotPoints(obj, reduction)
  plot <- ggplot(data, aes(x=axis1, y=axis2, color=cluster)) + geom_point(size = 0.5) + guides(colour = guide_legend(override.aes = list(size=2))) + ggtitle(name) + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  for (cl in levels(obj$seurat_clusters)) { #add cluster labels
    cluster.red <- subset(data, cluster == cl)
    plot <- plot + annotate("text", x = mean(cluster.red$axis1), y = mean(cluster.red$axis2), label = lbls[as.numeric(cl) + 1], size = 3, fontface=2)
  }
  return(plot)
}

generatePlotsByMetric <- function(obj, name, lbls, metric.name.seurat, metric.name, name.suffix, add.dir="", save.log2=TRUE) { #different plots for QC metrics
  #themes and axis labels for plots
  t <- scale_x_discrete(labels=lbls)
  t1 <- scale_y_discrete(labels=lbls)
  t2 <- theme(axis.text.x = element_text(angle = 45, size=10, hjust=1, face="bold"), axis.text.y = element_text(size=10), legend.position="none", axis.title.x=element_blank())
  t3 <- ggtitle(name)
  t4 <- theme(legend.position="none")
  t5 <- facet_wrap(. ~ clusters, ncol=5, labeller = as_labeller(lbls))
  t6 <- stat_summary(fun.y=mean, geom="point", shape=23, fill="blue", size=3)
  t7 <- theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10, face="bold"), legend.position="none", axis.title.y=element_blank())
  l1 <- labs(y=metric.name)
  l2 <- labs(y=paste0("log2(", metric.name, ")"))
  l3 <- labs(x=metric.name)
  l4 <- labs(x=paste0("log2(", metric.name, ")"))
  l5 <- labs(y=paste0("Average ", metric.name))
  l6 <- labs(y=paste0("Average log2(", metric.name, ")"))
  if (metric.name.seurat == "percent.mt") {
    a1 <- scale_y_continuous(breaks=seq(0, 80, 5))
    a2 <- scale_x_continuous(breaks=seq(0, 80, 5))
  }
  else {
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
  data$clusters = with(data, reorder(clusters, -metric, median)) #order data by cluster median
  
  #all plots except combined density and tsne/umap are saved in to vesions: no trasformation and log2, if save.log2 is true
  
  name.prefix <- paste0(results.dir, add.dir)
  
  if (save.log2) {
    #barplot of cluster means
    ggsave1(filename=paste0(name.prefix, "bar_mean_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0) %>% group_by(clusters) %>% summarize(metric = mean(log2(metric))), aes(x=clusters, y=log2(metric))) + geom_bar(aes(fill=clusters), stat="identity") + t + t2 + t3 + c1 + c2 + l6, n.clusters = n.clusters)
    #boxplot by cluster
    ggsave1(filename=paste0(name.prefix, "box_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=clusters, y=log2(metric))) + geom_boxplot(aes(fill=clusters)) + t + t2 + t3 + t6 + c1 + c2 + l2, n.clusters = n.clusters) 
    #combined density plots for each cluster
    ggsave1(filename=paste0(name.prefix, "density_", name.suffix, ".pdf"), plot=ggplot(subset(data, metric > 0), aes(x=log2(metric))) + geom_density(aes(fill=clusters)) + t3 + t4 + t5 + c1 + c2 + l4, n.clusters = n.clusters)
    #joyplot by cluster
    ggsave1(filename=paste0(name.prefix, "density2_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=log2(metric), y=clusters)) + geom_density_ridges(aes(fill=clusters)) + t1 + t3 + t7 + c1 + c2 + l4, n.clusters = n.clusters) 
    #overall destiny plot
    ggsave1(filename=paste0(name.prefix, "density3_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=log2(metric))) + geom_density(aes(fill="red")) + t3 + l4, n.clusters = n.clusters) 
    #violin plot by cluster
    ggsave1(filename=paste0(name.prefix, "violin_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=clusters, y=log2(metric))) + geom_violin(aes(fill=clusters)) + t + t2 + t3 + t6 + c1 + c2 + l2, n.clusters = n.clusters)
  }
  
  #barplot of cluster means
  ggsave1(filename=paste0(name.prefix, "bar_mean_", name.suffix, ".pdf"), plot=ggplot(data %>% group_by(clusters) %>% summarize(metric = mean(metric)), aes(x=clusters, y=metric)) + geom_bar(aes(fill=clusters), stat="identity") + t + t2 + t3 + l5 + c1 + c2 + a1, n.clusters = n.clusters)
  #boxplot by cluster
  ggsave1(filename=paste0(name.prefix, "box_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=clusters, y=metric)) + geom_boxplot(aes(fill=clusters)) + t + t2 + t3 + t6 + c1 + c2 + l1 + a1, n.clusters = n.clusters) 
  #joyplot by cluster
  ggsave1(filename=paste0(name.prefix, "density2_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=metric, y=clusters)) + geom_density_ridges(aes(fill=clusters)) + t1 + t3 + t7 + c1 + c2 + l3 + a2, n.clusters = n.clusters) 
  #overall destiny plot
  ggsave1(filename=paste0(name.prefix, "density3_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=metric)) + geom_density(aes(fill="red")) + t3 + l3 + a2, n.clusters = n.clusters) 
  #violin plot by cluster
  ggsave1(filename=paste0(name.prefix, "violin_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=clusters, y=metric)) + geom_violin(aes(fill=clusters)) + t + t2 + t3 + t6 + c1 + c2 + l1 + a1, n.clusters = n.clusters)
  
  #tsne and umap continious dimplots
  ggsave1(filename=paste0(name.prefix, "tsne_", name.suffix, ".pdf"), plot=DimPlotContinuous(obj, metric.name.seurat, lbls, name, "tsne"))
  ggsave1(filename=paste0(name.prefix, "umap_", name.suffix, ".pdf"), plot=DimPlotContinuous(obj, metric.name.seurat, lbls, name, "umap"))
  
  #signature scatter plots
  if (grepl("cd", metric.name)) {
    cor1 <- round(cor(obj$percent.mt, data$metric, method = "pearson"), 2)
    ggsave1(filename=paste0(name.prefix, "/scatter-mito_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=obj$percent.mt, y=metric)) + geom_point() + geom_smooth() + l1 + ggtitle(paste(name, cor1)), n.clusters = n.clusters)
  }
}

generatePlots <- function(obj, name, cell.types, annotations, sig.plots) { #main plots function
  message("Making Plots")
  lbls <- NULL #create labels in the following format: cluster #, Panglao Cell Type \n annotated Cell Type
  for (i in 1:length(cell.types)) {
    lbls <- c(lbls, paste0((i - 1), " ", cell.types[i], "\n", annotations[i]))
  }
  names(lbls) <- 0:(length(lbls) - 1) #rename labels with cluster #
  
  generatePlotsByMetric(obj, name, lbls, "nCount_RNA", "Number of UMIS", "count") #nUMI plots
  generatePlotsByMetric(obj, name, lbls, "nFeature_RNA", "Number of Genes", "genes") #nGenes plots
  generatePlotsByMetric(obj, name, lbls, "percent.mt", "percent.mt", "mito") #%mito plots
  generatePlotsByMetric(obj, name, lbls, "percent.rb", "percent.rb", "ribo") #%ribo plots
  
  if (sig.plots) {
    #signatures plots
    generatePlotsByMetric(obj, name, lbls, "cd11", "cd1", "cd1", "signatures/", save.log2 = FALSE) 
    generatePlotsByMetric(obj, name, lbls, "cd21", "cd2", "cd2", "signatures/", save.log2 = FALSE)
    generatePlotsByMetric(obj, name, lbls, "cd31", "cd3", "cd3", "signatures/", save.log2 = FALSE)
  }
  
  #cluster colored dimplots
  ggsave1(filename=paste0(results.dir, "/tsne_clusters.pdf"), plot=DimPlotCluster(obj, lbls, name, "tsne"))
  ggsave1(filename=paste0(results.dir, "/umap_clusters.pdf"), plot=DimPlotCluster(obj, lbls, name, "umap"))
  
}

generateMarkerPlots <- function(obj, top5, name) { #generate plots of top marker genes. Does not work properly
  top5 <- obj.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
  dir.create(paste0(results.dir, name, "/marker_plots/"), showWarnings = FALSE)
  for (c in unique(top5$cluster)) {
    ggsave1(filename=paste0(results.dir, name, "/marker_plots/", c, "-violin.pdf"), plot=VlnPlot(obj, features = subset(top5, subset = cluster == c)$gene, pt.size = 0.25))
    ggsave1(filename=paste0(results.dir, name, "/marker_plots/", c, "-scatter.pdf"), plot=FeaturePlot(obj, features = subset(top5, subset = cluster == c)$gene, pt.size = 0.25))
  }
}


#clustering & finding markers 
clusterize <- function(obj, res, compute.reductions=TRUE, compute.markers=TRUE) { #function that will perform standart clustering procedures. If compute.markers == TRUE will compute DE genes. If compute.reductions = TRUE will calculate TSNE and UMAP
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(x = obj)
  obj <- ScaleData(obj, features = all.genes)
  obj <- RunPCA(obj, features = VariableFeatures(obj))
  
  obj <- FindNeighbors(obj, dims = 1:50)
  obj <- FindClusters(obj, resolution = res)
  
  if (compute.reductions) {
    obj <- RunUMAP(obj, dims = 1:50)
    obj <- tryCatch({RunTSNE(obj, dims = 1:50, check_duplicates = FALSE)}, error = function(e) {RunTSNE(obj, dims = 1:10, check_duplicates = FALSE)})
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


#add cell death module scores
addCDScores <- function(obj) {
  message("Calculating CD scores")
  
  signatures.path <- paste0(data.dir, "signatures/")
  if (is.human) {
    cd1 <- toupper(as.character(read.csv(paste0(signatures.path, "cd1_signatures.csv"), header = FALSE)$V1))
    cd2 <- toupper(as.character(read.csv(paste0(signatures.path, "cd2_signatures.csv"), header = FALSE)$V1))
    cd3 <- toupper(as.character(read.csv(paste0(signatures.path, "cd3_signatures.csv"), header = FALSE)$V1))
  }
  else {
    cd1 <- toTitleCase(tolower(as.character(read.csv(paste0(signatures.path, "cd1_signatures.csv"), header = FALSE)$V1)))
    cd2 <- toTitleCase(tolower(as.character(read.csv(paste0(signatures.path, "cd2_signatures.csv"), header = FALSE)$V1)))
    cd3 <- toTitleCase(tolower(as.character(read.csv(paste0(signatures.path, "cd3_signatures.csv"), header = FALSE)$V1)))
  }
  
  obj <- AddModuleScore(obj, features=list(cd1), name="cd1")
  obj <- AddModuleScore(obj, features=list(cd2), name="cd2")
  obj <- AddModuleScore(obj, features=list(cd3), name="cd3")
  
  return(obj)
}

#annotations & cell types assingment
formatMarkers <- function(lst) {
  st <- ""
  for (marker in unlist(lst)) {
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
  for (cl in levels(obj$seurat_clusters)) {
    cluster <- subset(obj, idents = cl)
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
      percents2 <- c(percents2, second.frequent / length(colnames(cluster$RNA)))
    }
    percents1 <- c(percents1, first.frequent / length(colnames(cluster$RNA)))
  }
  return(list("a1" = cluster.labels, "a2" = cluster.labels2, "p1" = percents1, "p2" = percents2))
}

assignCellTypes <- function(obj, markers, annotations, record.stats=FALSE) { #function that assigns cell types based on marker genes using gene to cell.type dictionary
  message("Assigning Cell Types")
  genes <- read.csv(paste0(data.dir, "markers.tsv"), sep="\t") #read cell type markers
  genes <- data.frame(genes)
  genes <- select(genes, official.gene.symbol, cell.type, organ)
  
  markers <- markers %>% group_by(cluster) %>% filter(avg_logFC > 0, p_val_adj < 0.05) #leave only markers that have positive fold change and pval < 0.05
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
    obj.cluster <- subset(obj, idents = cl)
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
      cells <- c(cells, length(colnames(obj.cluster$RNA)))
      genes.mean <- c(genes.mean, mean(obj.cluster$nFeature_RNA))
      genes.median <- c(genes.median, median(obj.cluster$nFeature_RNA))
      mito.mean <- c(mito.mean, mean(obj.cluster$percent.mt))
      mito.median <- c(mito.median, median(obj.cluster$percent.mt))
      
      formatted.markers <- c(formatted.markers, formatMarkers(cluster.markers$gene)) #format markers for csv recording by separating them with ;
    }
  }
  if (record.stats) {
    sm <- paste(sum(scores) / length(cluster.ids), sum(scores2) / length(cluster.ids), sum(cells), round(mean(genes.mean), 3), round(mean(mito.mean), 3), sep=",") #add score summary to the file
    clusters <- tibble("cluster" = cluster.ids, "annotation" = annotations[["a1"]], "annotation2" = annotations[["a2"]], "%annotation1" = round(annotations[["p1"]], 3), 
                  "%annotation2" = round(annotations[["p2"]], 3), "cell.type" = cell.types, "cell.type2" = cell.types2, "sum.avg_logFC" = round(scores, 3), 
                  "delta.score" = round(scores2, 3), "cells" = cells, "genes.mean" = round(genes.mean, 3), "genes.median" = round(genes.median, 3), 
                  "mito.mean" = round(mito.mean, 3), "mito.median" = round(mito.median, 3), "markers" = formatted.markers, "score.genes" = score.genes)
    return(list("clusters" = clusters, "sm" = sm, "markers" = markers))
  }
  else { #return only cell types and annotations
    res <- tibble("cluster" = cluster.ids, "annotation" = annotations[["a1"]], "cell.type" = cell.types)
    return(res)
  }
}


saveResults <- function(obj, clusters, obj.markers, mc_specific=FALSE, sm=NA) {
  message("Writing Results")

  write.csv(clusters, file=paste0(results.dir, "/!clusters.csv"))
  write.csv(obj.markers, file=paste0(results.dir, "/!markers.csv"))
  
  all.cells <- colnames(obj)
  umap <- as.data.frame(Embeddings(obj$umap)[all.cells, c(1, 2)])
  tsne <- as.data.frame(Embeddings(obj$tsne)[all.cells, c(1, 2)])
  
  cells <- data.frame("cell"=all.cells, "cluster"=obj$seurat_clusters, "cell.type"=clusters$cell.type[obj$seurat_clusters], 
                     "annotation" = clusters$annotation[obj$seurat_clusters], "UMAP_1"=umap$UMAP_1, "UMAP_2"=umap$UMAP_2,
                     "TSNE_1"=tsne$tSNE_1, "TSNE_2"=tsne$tSNE_2)
  write.csv(cells, file=paste0(results.dir, "/!cells.csv"))
  
  if (mc_specific) {
    if (res == 0.5 && method == "none") { 
      message("Recording Additional Stats")
      #record all cells with QC stats into csv for summary plots
      t <- tibble("tissue" = tissue, "cluster" = obj$seurat_clusters, "nCount_RNA" = obj$nCount_RNA, "nFeature_RNA" = obj$nFeature_RNA, "percent.mt" = obj$percent.mt, "percent.rb" = obj$percent.rb)
      write.table(t, paste0(results.dir, "../../stats_summary.csv"), sep=",", append=TRUE, col.names=FALSE)
      
      sm <- paste(task.name, sm, sep=",")
      write(sm, paste0(results.dir, "../score_summary.csv"), append=TRUE)
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


#main function
MCMain <- function() {
  tasks.per.res <- tasks.per.tiss #how many different methods per one resolution
  res <<- 1 #0.5 * (1 + ((task.id %% tasks.per.tiss) %/% tasks.per.res)) #clustering resolution
  method <<- switch(task.id %% tasks.per.res + 1, "none", "z_score", "cutoff", "outlier", "mad") #filtering method
  param <<- switch(task.id %% tasks.per.res + 1, 0, 2, 10, 0, 2) #filtering parameter
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
  
  tiss <<- addCDScores(tiss)
  
  tmp <- assignCellTypes(tiss, obj.markers, getAnnotations(tiss), record.stats = TRUE) #assign cell types
  #unpack returned object
  clusters <<- tmp$clusters
  sm <- tmp$sm
  obj.markers <<- tmp$markers
  
  generatePlots(tiss, task.name, clusters$cell.type, clusters$annotation, sig.plots = TRUE) #make plots
  #generateMarkerPlots(obj, filename)
  
  saveResults(tiss, clusters, obj.markers, mc_specific=TRUE, sm=sm)
  
  message(paste0("Finished task.id:", info.msg))
}
