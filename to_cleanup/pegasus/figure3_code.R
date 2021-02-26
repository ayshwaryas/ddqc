data.from.pg <- TRUE

source("../scripts/settings.R")
source("../scripts/local_settings.R")


#plots
ggsave1 <- function(filename, plot, n.clusters=30, type="h", custom.size=TRUE) {
  if (type == "h") {
    height = 10
    width = 14 / 30 * max(n.clusters, 30)
  }
  if (type == "v") {
    height = 10 / 30 * max(n.clusters, 30)
    width = 14
  }
  no_bkg <- theme(axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank()) 
  if (custom.size) {
    ggsave(filename = filename, plot = plot + no_bkg, width = width, height = height) #saves plot with custom dimensions 
  } else {
    ggsave(filename = filename, plot = plot + no_bkg) #saves plot withput custom dimensions 
  }
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
  t <- theme(axis.text.x = element_text(size=15), axis.title.x = element_text(size=15), axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), plot.title = element_text(size = 20, face = "bold"))
  cols <- scale_colour_gradientn(colours=rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")))
  if (log2) {
    plot <- ggplot(data, aes(x=axis1, y=axis2, color=log2(color))) + geom_point(size = 1.5) + labs(color=paste0("log2(", metric.name, ")")) + t + cols + ggtitle(name) + labs(x = "UMAP1", y = "UMAP2")
  } else {
    plot <- ggplot(data, aes(x=axis1, y=axis2, color=color)) + geom_point(size = 1.5) + labs(color=metric.name) + t + cols + ggtitle(name)
  }
  for (cl in levels(obj$seurat_clusters)) { #add cluster labels
    cluster.red <- subset(data, cluster == cl)
    plot <- plot + annotate("text", x = mean(cluster.red$axis1), y = mean(cluster.red$axis2), label = lbls[as.numeric(cl) + 1], size = 5, fontface=2)
  } 
  return(plot)
}

DimPlotCluster <- function(obj, lbls, name, reduction) { #DimPlot colored by cluster
  name <- paste0(name, "_", reduction)
  data <- GetDimPlotPoints(obj, reduction, "seurat_clusters")
  t <- theme(axis.text.x = element_text(size=15), axis.title.x = element_text(size=15), axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), plot.title = element_text(size = 20, face = "bold"))
  plot <- ggplot(data, aes(x=axis1, y=axis2, color="Red")) + geom_point(size = 1.5) + ggtitle(name) + t + labs(x = "UMAP1", y = "UMAP2")
  for (cl in levels(obj$seurat_clusters)) { #add cluster labels
    cluster.red <- subset(data, cluster == cl)
    plot <- plot + annotate("text", x = mean(cluster.red$axis1), y = mean(cluster.red$axis2), label = lbls[as.numeric(cl) + 1], size = 5, fontface=2)
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
    ggsave1(filename=paste0(name.prefix, "density2_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=log2(metric), y=clusters)) + geom_density_ridges(aes()) + t1 + t3 + t6h + t7 + c1 + c2 + l4 + vl, n.clusters = n.clusters) 
    #violin plot by cluster
    ggsave1(filename=paste0(name.prefix, "violin_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=clusters, y=log2(metric))) + geom_violin(aes()) + t + t2 + t3 + t6 + c1 + c2 + l2 + hl, n.clusters = n.clusters)
    #umap
    ggsave1(filename=paste0(name.prefix, "umap_", name.suffix, ".pdf"), plot=DimPlotContinuous(obj, metric.name.seurat, lbls, name, "umap", log2 = TRUE))
  } else {
    data$clusters = with(data, reorder(clusters, -metric, median)) #order data by cluster mean
    
    #boxplot by cluster
    ggsave1(filename=paste0(name.prefix, "box_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=clusters, y=metric)) + geom_boxplot(aes()) + t + t2 + t3 + t6 + c1 + c2 + l1 + a1 + hl, n.clusters = n.clusters) 
    #joyplot by cluster
    ggsave1(filename=paste0(name.prefix, "density2_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=metric, y=clusters)) + geom_density_ridges(aes()) + t1 + t3 + t6h + t7 + c1 + c2 + l3 + a2 + vl, n.clusters = n.clusters) 
    #violin plot by cluster
    ggsave1(filename=paste0(name.prefix, "violin_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=clusters, y=metric)) + geom_violin(aes()) + t + t2 + t3 + t6 + c1 + c2 + l1 + a1 + hl, n.clusters = n.clusters)
    #umap
    ggsave1(filename=paste0(name.prefix, "umap_", name.suffix, ".pdf"), plot=DimPlotContinuous(obj, metric.name.seurat, lbls, name, "umap"))
  }
}

generatePlots <- function(obj, name, cell.types, annotations) { #main plots function
  message("Making Plots")
  lbls <- NULL #create labels in the following format: cluster #, Panglao Cell Type \n annotated Cell Type
  for (i in 1:length(clusters$cell.type)) {
    lbls <- c(lbls, paste0((i - 1), " ", clusters$cell.type[i]))#, "\n", clusters$annotation[i]))
  }
  names(lbls) <- 0:(length(lbls) - 1)  #rename labels with cluster #
  
  generatePlotsByMetric(obj, name, lbls, "nCount_RNA", "Number of UMIS", "count") #nUMI plots
  generatePlotsByMetric(obj, name, lbls, "nFeature_RNA", "Number of Genes", "genes") #nGenes plots
  generatePlotsByMetric(obj, name, lbls, "percent.mt", "percent.mt", "mito") #%mito plots
  generatePlotsByMetric(obj, name, lbls, "percent.rb", "percent.rb", "ribo") #%ribo plots
  
  #cluster colored dimplots
  ggsave1(filename=paste0(results.dir, "/umap_clusters.pdf"), plot=DimPlotCluster(obj, lbls, name, "umap"))
}


#function that makes filtered cells plots
generateFCPlots <- function(obj, name, clusters) {
  message("Making FC Plots")
  lbls <- NULL #create labels in the following format: cluster #, Panglao Cell Type \n annotated Cell Type
  for (i in 1:length(clusters$cell.type)) {
    lbls <- c(lbls, paste0((i - 1), " ", clusters$cell.type[i]))#, "\n", clusters$annotation[i]))
  }
  names(lbls) <- 0:(length(lbls) - 1)  #rename labels with cluster #
  
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
  
  plot.cols <- c("C10 only" = "#DB6400", "MAD2 only" = "#16697A","All" = "#FFA62B" , "Neither" = "#99A8B2") #to keep consistent plot colors
  
  
  data <- data.frame(UMAP_1 = data$UMAP_1, UMAP_2 = data$UMAP_2, cluster = obj$seurat_clusters, color=obj$color, annotation=obj$annotations)
  
  t1 <- theme(axis.text.x = element_text(size=15), axis.title.x = element_text(size=15), axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), plot.title = element_text(size = 20, face = "bold"))
  t2 <- guides(colour = guide_legend(override.aes = list(size=2)))
  
  plot1 <- ggplot(data, aes(x=UMAP_1, y=UMAP_2, color=color)) + geom_point(size = 1.5) + t1 + t2 + scale_fill_manual(values = plot.cols) + scale_color_manual(values = plot.cols) #umap colored by filtering category
  for (cl in levels(obj$seurat_clusters)) { #add labels to plots
    cluster.umap <- subset(data, cluster == cl)
    plot1 <- plot1 + annotate("text", x = mean(cluster.umap$UMAP_1), y = mean(cluster.umap$UMAP_2), label = lbls[(as.numeric(cl) + 1)], size = 5, fontface=2)
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
  ggsave1(filename = paste0(results.dir, res, "-!filterplot.pdf"), plot=plot1)
  ggsave1(filename = paste0(results.dir, res, "-!barplot.pdf"), plot=plot3, n.clusters=n.clusters)
}


PATH <- "~/Downloads/figure3_plots/"

for (tissue in list.files(PATH)) {
  
  res <<- 1.4
  metric <<- "none"
  
  results.dir <<-  paste0(PATH, tissue, "/")
  
  tiss <- read.csv(paste0(results.dir, "!cells.csv"))
  clusters <- read.csv(paste0(results.dir, "!clusters.csv"))
  tiss <- tiss %>% rename(nFeature_RNA = n_genes, nCount_RNA = n_counts, percent.mt = percent_mito, percent.rb = percent_ribo, seurat_clusters = louvain_labels)
  tiss$seurat_clusters <- factor(tiss$seurat_clusters - 1)
  
  if (grepl("joint", "lung_joint", fixed = TRUE)) {
    generateFCPlots(tiss, tissue, clusters)
  } else {
    generatePlots(tiss, tissue, clusters$cell.type, clusters$annotation)
  }
}
