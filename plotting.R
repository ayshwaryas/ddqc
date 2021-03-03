library(ggplot2)
library(ggridges)
library(cowplot)

theme_horizontal <- theme(axis.text.x = element_text(angle = 45, size=15, hjust=1, face="bold"), 
                          axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), 
                          legend.position="none", axis.title.x=element_blank(), plot.title = element_text(size = 20, face = "bold"))
theme_horizontal_with_legend <- theme(axis.text.x = element_text(angle = 45, size=15, hjust=1, face="bold"), 
                          axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), 
                          axis.title.x=element_blank(), plot.title = element_text(size = 20, face = "bold"))
theme_vertical <- theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15, face="bold"), 
                        axis.title.x = element_text(size=15), legend.position="none", axis.title.y=element_blank(), 
                        plot.title = element_text(size = 20, face = "bold"))
theme_umap <- theme(axis.text.x = element_text(size=15), axis.title.x = element_text(size=15),
           axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), 
           plot.title = element_text(size = 20, face = "bold"), legend.title = element_text(size = 15), 
           legend.text = element_text(size = 10))



ggsave1 <- function(filename, plot, n.clusters=30, type="h") {  # custom ggsave function
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

DimPlotContinuous <- function(obj, metric.name.pg, metric.name, lbls, log2=FALSE) { #DimPlot with continious colors by metric
  cols <- scale_colour_gradientn(colours=rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", 
                                               "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", 
                                               "#5E4FA2")))
  data <- data.frame(UMAP1 = obj$umap1, UMAP2 = obj$umap2, color = obj[[metric.name.pg]], cluster = obj$louvain_labels)
  if (log2) {
    plot <- ggplot(data, aes(x=UMAP1, y=UMAP2, color=log(color), fill=cluster)) + geom_point(size = 1) + theme_umap + ttl + cols + scale_fill_discrete(labels = lbls)
    plot <- plot + labs(color=paste0("log2(", metric.name, ")")) 
  } else {
    plot <- ggplot(data, aes(x=UMAP1, y=UMAP2, color=color, fill=cluster)) + geom_point(size = 1) + theme_umap + ttl + cols + scale_fill_discrete(labels = lbls)
    plot <- plot + labs(color=metric.name)
  }
  for (cl in levels(obj$louvain_labels)) { #add cluster labels
    cluster.data <- subset(data, cluster == cl)
    plot <- plot + annotate("text", x = mean(cluster.data$UMAP1), y = mean(cluster.data$UMAP2), label=cl, size = 7, fontface=2)
  }
  return(plot)
}

DimPlotCluster <- function(obj, lbls) { #DimPlot colored by cluster
  data <- data.frame(UMAP1 = obj$umap1, UMAP2 = obj$umap2, cluster = obj$louvain_labels)
  plot <- ggplot(data, aes(x=UMAP1, y=UMAP2, color=cluster)) + geom_point(size = 1) + theme_umap + ttl
  for (cl in levels(obj$louvain_labels)) { #add cluster labels
    cluster.data <- subset(data, cluster == cl)
    plot <- plot + annotate("text", x = mean(cluster.data$UMAP1), y = mean(cluster.data$UMAP2), label=cl, size = 7, fontface=2)
  }
  return(plot)
}

generatePlotsByMetric <- function(obj, lbls, metric.name.pg, metric.name, name.suffix) { #different plots for QC metrics
  labels_horizontal <- scale_x_discrete(labels=lbls)
  labels_vertical <- scale_y_discrete(labels=lbls)
  
  mean_plot <- stat_summary(fun=mean, geom="point", shape=23, fill="blue", size=3)
  
  if (metric.name.pg == "percent_mito") {
    axis_breaks_vertical <- scale_x_continuous(breaks=seq(0, 80, 5))
    axis_breaks_horizontal <- scale_y_continuous(breaks=seq(0, 80, 5))
    horizontal_line <- geom_hline(yintercept=10, color="red", size=0.5)
    vertical_line <- geom_vline(xintercept=10, color="red", size=0.5)
  } else {
    if (metric.name.pg == "n_genes") {
      horizontal_line <- geom_hline(yintercept=log2(200), color="red", size=0.5)
      vertical_line <- geom_vline(xintercept=log2(200), color="red", size=0.5)
    } else { 
      horizontal_line <- NULL
      vertical_line <- NULL
    }
    axis_breaks_vertical <- NULL
    axis_breaks_horizontal <- NULL
  }
  
  n.clusters <- length(levels(obj$louvain_labels))
  data <- data.frame(metric=obj[[metric.name.pg]], clusters=obj$louvain_labels) 
  colnames(data) <- c("metric", "clusters") #rename data columns
  
  if (metric.name.pg == "n_counts" || metric.name.pg == "n_genes") {
    data$clusters = with(data, reorder(clusters, -log2(metric), mean)) #order data by cluster mean
    data$metric <- log2(data$metric)
    axis_labels_horizontal <- labs(y=paste0("log2(", metric.name, ")"))
    axis_labels_vertical <- labs(x=paste0("log2(", metric.name, ")"))
    name.suffix <- paste0(name.suffix, "_log.pdf")
    is.log2 = TRUE
  } else {
    data$clusters = with(data, reorder(clusters, -metric, mean)) #order data by cluster mean
    axis_labels_horizontal <- labs(y=metric.name)
    axis_labels_vertical <- labs(x= metric.name)
    name.suffix <- paste0(name.suffix, ".pdf")
    is.log2 = FALSE
  }
  
  #boxplot by cluster
  boxplot <- ggplot(subset(data, metric > 0), aes(x=clusters, y=metric)) + geom_boxplot() + theme_horizontal + mean_plot + labels_horizontal + horizontal_line + axis_breaks_horizontal + ttl
  #joyplot by cluster
  joyplot <- ggplot(subset(data, metric > 0), aes(x=metric, y=clusters)) + geom_density_ridges(aes()) + theme_vertical + mean_plot + labels_vertical + vertical_line + axis_breaks_vertical + ttl
  #violin plot by cluster
  vnlplot <- ggplot(subset(data, metric > 0), aes(x=clusters, y=metric)) + geom_violin(aes()) + theme_horizontal + mean_plot + labels_horizontal + horizontal_line + axis_breaks_horizontal + ttl
  #umap
  umapplot <- DimPlotContinuous(obj, metric.name.pg, metric.name, lbls, log2 = is.log2)
  
  #save plots 
  ggsave1(filename=paste0(results.dir, "box_", name.suffix), plot=boxplot + axis_labels_horizontal, n.clusters = n.clusters, type = "h") 
  ggsave1(filename=paste0(results.dir, "density_", name.suffix), plot=joyplot + axis_labels_vertical, n.clusters = n.clusters, type = "v") 
  ggsave1(filename=paste0(results.dir, "violin_", name.suffix), plot=vnlplot + axis_labels_horizontal, n.clusters = n.clusters, type = "h")
  ggsave1(filename=paste0(results.dir, "umap_", name.suffix), plot=umapplot, n.clusters = n.clusters, type = "u")
}

generateFCPlots <- function(obj, lbls) {
  plot.cols <- c("Cutoff only" = "#4A3933", "MAD2 only" = "#16697A",
                 "All" = "#FFA62B" , "Neither" = "#AAAAAA") #to keep consistent plot colors
  data <- data.frame(UMAP1 = obj$umap1, UMAP2 = obj$umap2, cluster = obj$louvain_labels, color=obj$color, annotation=obj$annotations)
  
  
  fltplot <- ggplot(data, aes(x=UMAP1, y=UMAP2, color=color, fill=cluster)) + geom_point(size = 1) + 
    theme_umap + guides(colour = guide_legend(override.aes = list(size=2))) + 
    scale_fill_discrete(labels = lbls) + scale_color_manual(values = plot.cols) + ttl
   
  for (cl in levels(obj$louvain_labels)) { #add cluster labels
    cluster.data <- subset(data, cluster == cl)
    fltplot <- fltplot + annotate("text", x = mean(cluster.data$UMAP1), y = mean(cluster.data$UMAP2), label=cl, size = 7, fontface=2)
  }
  
  #create barplot which shows category distribution
  table.color <- NULL
  table.cluster <- NULL
  table.freq <- NULL
  for (cl in levels(data$cluster)) {
    data.cluster <- subset(data, cluster == cl)
    table.tmp <- as.data.frame(table(data.frame(color=factor(data.cluster$color), cluster = as.character(data.cluster$cluster))))
    table.tmp$Freq <- table.tmp$Freq / sum(table.tmp$Freq)
    table.color <- c(table.color, as.character(table.tmp$color))
    table.cluster <- c(table.cluster, as.integer(as.character(table.tmp$cluster)))
    table.freq <- c(table.freq, as.character(table.tmp$Freq))
  }
  data1 <- data.frame(color=table.color, cluster = as.factor(table.cluster), freq=as.double(as.character(table.freq)) * 100)
  
  freqplot <- ggplot(data1, aes(x=cluster, y=freq, fill=color)) + geom_bar(stat="identity") + theme_horizontal_with_legend + ttl
  
  #write plots
  n.clusters <- length(levels(obj$louvain_labels))
  ggsave1(filename = paste0(results.dir, "!p_filterplot.pdf"), plot=fltplot, n.clusters = n.clusters, type = "u")
  ggsave1(filename = paste0(results.dir, "!p_barplot.pdf"), plot=freqplot, n.clusters=n.clusters, type = "h")
}

generatePlots <- function(obj, cell.types, joint=FALSE) { #main plots function
  message("Making Plots")
  lbls <- NULL #create labels in the following format: cluster #, Panglao Cell Type \n annotated Cell Type
  for (i in 1:length(cell.types)) {
    lbls <- c(lbls, paste0(i, " ", cell.types[i]))
  }
  ttl <<- ggtitle(task.name)
  
  generatePlotsByMetric(obj, lbls, "n_counts", "Number of UMIS", "count") #nUMI plots
  generatePlotsByMetric(obj, lbls, "n_genes", "Number of Genes", "genes") #nGenes plots
  generatePlotsByMetric(obj, lbls, "percent_mito", "percent_mito", "mito") #%mito plots
  generatePlotsByMetric(obj, lbls, "percent_ribo", "percent_ribo", "ribo") #%ribo plots
  
  #cluster colored dimplots
  ggsave1(filename=paste0(results.dir, "umap_clusters.pdf"), plot=DimPlotCluster(obj, lbls), n.clusters = length(levels(obj$louvain_labels)), type = "u")
  
  if (joint) { # make joint clustering plots
    generateFCPlots(obj, lbls)
  }
}

task.name <- commandArgs(trailingOnly = TRUE)[1]
results.dir <- commandArgs(trailingOnly = TRUE)[2]
script.type <- commandArgs(trailingOnly = TRUE)[3]

message("Starting R script to generate plots")

tiss <- read.csv(paste0(results.dir, "!cells.csv"))
tiss$louvain_labels <- as.factor(tiss$louvain_labels)
clusters <- read.csv(paste0(results.dir, "!clusters.csv"))

generatePlots(tiss, clusters$cell_type, script.type == "joint") #make plots
