library(ggplot2)
library(ggridges)
library(cowplot)
library(dplyr)

PATH <- "/Volumes/easystore/primes_storage/"
DATA.DIR <- paste0(PATH, "output_pg/")

theme_horizontal <- theme(axis.text.x = element_text(angle = 45, size=15, hjust=1, face="bold"), 
                          axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), 
                          legend.position="none", axis.title.x=element_blank(), plot.title = element_text(size = 20, face = "bold"))
theme_horizontal_with_legend <- theme(axis.text.x = element_text(angle = 45, size=15, hjust=1, face="bold"), 
                                      axis.text.y = element_text(size=15), axis.title.y = element_text(size=15), 
                                      axis.title.x=element_blank(), plot.title = element_text(size = 20, face = "bold"))
theme_vertical <- theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15, face="bold"), 
                        axis.title.x = element_text(size=15), legend.position="none", axis.title.y=element_blank(), 
                        plot.title = element_text(size = 20, face = "bold"))


ggsave1 <- function(filename, plot, n.tissues=30, type="h") {  # custom ggsave function
  if (type == "h") {
    height = 10
    width = 14 / 30 * max(n.tissues, 30)
  }
  if (type == "v") {
    height = 10 / 30 * max(n.tissues, 30)
    width = 14
  }
  if (type == "u") {
    height = 10
    width = 10 + 2 * ceiling(n.tissues / 13)
  }
  no_bkg <- theme(axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank()) 
  ggsave(filename = filename, plot = plot + no_bkg, width = width, height = height) #saves plot with custom dimensions 
}

generatePlotsByMetric <- function(obj, metric.name.pg, metric.name, name.suffix) { #different plots for QC metrics
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
  
  n.tissues <- length(levels(obj$tissue))
  data <- data.frame(metric=obj[[metric.name.pg]], tissues=obj$tissue) 
  colnames(data) <- c("metric", "tissues") #rename data columns
  
  if (metric.name.pg == "n_counts" || metric.name.pg == "n_genes") {
    data$tissues = with(data, reorder(tissues, -log2(metric), mean)) #order data by tissue mean
    data$metric <- log2(data$metric)
    axis_labels_horizontal <- labs(y=paste0("log2(", metric.name, ")"))
    axis_labels_vertical <- labs(x=paste0("log2(", metric.name, ")"))
    name.suffix <- paste0(name.suffix, "_log.pdf")
    is.log2 = TRUE
  } else {
    data$tissues = with(data, reorder(tissues, -metric, mean)) #order data by tissue mean
    axis_labels_horizontal <- labs(y=metric.name)
    axis_labels_vertical <- labs(x= metric.name)
    name.suffix <- paste0(name.suffix, ".pdf")
    is.log2 = FALSE
  }
  
  #boxplot by tissue
  boxplot <- ggplot(subset(data, metric > 0), aes(x=tissues, y=metric)) + geom_boxplot() + theme_horizontal + mean_plot + horizontal_line + axis_breaks_horizontal + ttl
  #joyplot by tissue
  joyplot <- ggplot(subset(data, metric > 0), aes(x=metric, y=tissues)) + geom_density_ridges(jittered_points = TRUE,
                                                                                              position = position_points_jitter(width = 0.05, height = 0),
                                                                                              point_shape = '|', point_size = 1, point_alpha = 0.5, alpha = 0.7) + theme_vertical + mean_plot + vertical_line + axis_breaks_vertical + ttl
  #violin plot by tissue
  vnlplot <- ggplot(subset(data, metric > 0), aes(x=tissues, y=metric)) + geom_violin() + theme_horizontal + mean_plot + horizontal_line + axis_breaks_horizontal + ttl

  
  #save plots 
  ggsave1(filename=paste0(results.dir, "box_", name.suffix), plot=boxplot + axis_labels_horizontal, n.tissues = n.tissues, type = "h") 
  ggsave1(filename=paste0(results.dir, "density_", name.suffix), plot=joyplot + axis_labels_vertical, n.tissues = n.tissues, type = "v") 
  ggsave1(filename=paste0(results.dir, "violin_", name.suffix), plot=vnlplot + axis_labels_horizontal, n.tissues = n.tissues, type = "h")
}

generatePlots <- function(obj) { #main plots function
  message("Making Plots")
  
  generatePlotsByMetric(obj, "n_counts", "Number of UMIS", "count") #nUMI plots
  generatePlotsByMetric(obj, "n_genes", "Number of Genes", "genes") #nGenes plots
  generatePlotsByMetric(obj, "percent_mito", "percent_mito", "mito") #%mito plots
  generatePlotsByMetric(obj, "percent_ribo", "percent_ribo", "ribo") #%ribo plots
}

generatePlotsPG <- function(project) {
  data.path <- paste0(DATA.DIR, project, "/")
  results.dir <<- paste0(PATH, "figure1_plots/", project, "/")
  
  if (!file.exists(paste0(results.dir, "!cells.csv"))) {
    dataset <- NULL
    for (directory in list.dirs(data.path, full.names = FALSE, recursive = FALSE)) {
      if (!file.exists(paste0(data.path, directory, "/1.4-none-0/!cells.csv"))) {
        next
      }
      d <- read.csv(paste0(data.path, directory, "/1.4-none-0/!cells.csv"))
      d$tissue <- directory
      if (is.null(dataset)) {
        dataset <- d
      } else {
        dataset <- full_join(dataset, d)
      }
    }
    dataset <- select(dataset, "tissue", "barcodekey", "Channel", "annotations", "n_genes", "n_counts", "percent_mito", "percent_ribo")
    dir.create(paste0(PATH, "figure1_plots/"), showWarnings = FALSE)
    dir.create(paste0(results.dir), showWarnings = FALSE)
    write.csv(dataset, paste0(results.dir, "!cells.csv"))
  } else {
    dataset <- read.csv(paste0(results.dir, "!cells.csv"))
  }
  
  
  ttl <<- ggtitle(project)
  generatePlots(dataset)
}


for (prj in list.dirs(DATA.DIR, full.names = FALSE, recursive = FALSE)) {
  generatePlotsPG(prj)
}



