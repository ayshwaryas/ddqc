ggsave1 <- function(filename, plot, n.tissues=30, type="h") {
  if (type == "h") {
    height = 10
    width = 14 / 30 * max(n.tissues, 30)
  }
  if (type == "v") {
    height = 10 / 30 * max(n.tissues, 20)
    width = 14
  }
  no_bkg <- theme(axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank()) 
  ggsave(filename = filename, plot = plot + no_bkg, width = width, height = height) #saves plot with custom dimensions 
}

generatePlots <- function(results.dir, dataset, project) {
  for (metric in c("nCount_RNA,Number of UMIS,count", "nFeature_RNA,Number of Genes,genes", 
                   "percent.mt,percent.mt,mito", "percent.rb,percent.rb,ribo")) {
    metric.name.seurat <- strsplit(metric, ",")[[1]][1]
    metric.name <- strsplit(metric, ",")[[1]][2]
    name.suffix <- strsplit(metric, ",")[[1]][3]
    
    name <- project
    add.dir <- ""
    
    #themes and axis labels for plots
    t2 <- theme(axis.text.x = element_text(angle = 45, size=10, hjust=1, face="bold"), axis.text.y = element_text(size=10), legend.position="none", axis.title.x=element_blank())
    t3 <- ggtitle(name)
    t4 <- theme(legend.position="none")
    t6 <- stat_summary(fun.y=mean, geom="point", shape=23, fill="blue", size=3)
    t6h <- stat_summary(fun.x=mean, geom="point", shape=23, fill="blue", size=3)
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
    
    
    n.tissues <- length(unique(dataset$tissue))
    plot.cols <- scales::hue_pal(c = 100, l = 65, h.start = 0)(n.tissues)
    names(plot.cols) <- unique(dataset$tissue)
    c1 <- scale_fill_manual(values = plot.cols) 
    c2 <- scale_color_manual(values = plot.cols)
    
    data <- data.frame(metric=dataset[[metric.name.seurat]], tissues=dataset$tissue) 
    colnames(data) <- c("metric", "tissues") #rename data columns
    
    #all plots except combined density and tsne/umap are saved in to vesions: no trasformation and log2, if save.log2 is true
    
    name.prefix <- paste0(results.dir, add.dir)
    
    if (metric.name.seurat == "nCount_RNA" || metric.name.seurat == "nFeature_RNA") {
      data$tissues = with(data, reorder(tissues, -log2(metric), mean)) #order data by tissue log mean
      
      #boxplot by tissue
      ggsave1(filename=paste0(name.prefix, "box_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=tissues, y=log2(metric))) + geom_boxplot() + t2 + t3 + t6 + c1 + c2 + l2 + hl, n.tissues = n.tissues) 
      #joyplot by tissue
      ggsave1(filename=paste0(name.prefix, "density2_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=log2(metric), y=tissues)) + geom_density_ridges() + t3 + t6h + t7 + c1 + c2 + l4 + vl, n.tissues = n.tissues) 
      #violin plot by tissue
      ggsave1(filename=paste0(name.prefix, "violin_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=tissues, y=log2(metric))) + geom_violin() + t2 + t3 + t6 + c1 + c2 + l2 + hl, n.tissues = n.tissues)
    } else {
      data$tissues = with(data, reorder(tissues, -metric, mean)) #order data by tissue mean
      
      #boxplot by tissue
      ggsave1(filename=paste0(name.prefix, "box_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=tissues, y=metric)) + geom_boxplot() + t2 + t3 + t6 + c1 + c2 + l1 + a1 + hl, n.tissues = n.tissues) 
      #joyplot by tissue
      ggsave1(filename=paste0(name.prefix, "density2_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=metric, y=tissues)) + geom_density_ridges() + t3 + t6h + t7 + c1 + c2 + l3 + a2 + vl, n.tissues = n.tissues) 
      #violin plot by tissue
      ggsave1(filename=paste0(name.prefix, "violin_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=tissues, y=metric)) + geom_violin() + t2 + t3 + t6 + c1 + c2 + l1 + a1 + hl, n.tissues = n.tissues)
    }
  }
}

generatePlotsPG <- function(project) {
  data.from.pg <<- TRUE
  #PATH <- "/ahg/regevdata/projects/scqc/"
  PATH <- "~/Downloads/"
  source("../scripts/settings.R")
  source("../scripts/local_settings.R")
  data.path <- paste0(PATH, "figure1_data/", project, "/")
  dataset <- NULL
  for (directory in list.files(data.path)) {
    d <- read.csv(paste0(data.path, directory, "/!cells.csv"))
    d$tissue <- directory
    if (is.null(dataset)) {
      dataset <- d
    }
    else {
      dataset <- full_join(dataset, d)
    }
  }
  dataset <- dataset %>% rename(nFeature_RNA = n_genes, nCount_RNA = n_counts, percent.mt = percent_mito, percent.rb = percent_ribo)
  
  results.dir <- paste0(PATH, "figure1_plots/", project, "/")
  dir.create(paste0(PATH, "figure1_plots/"), showWarnings = FALSE)
  dir.create(paste0(results.dir), showWarnings = FALSE)
  write.csv(dataset, paste0(results.dir, "!cells.csv"))
  
  generatePlots(results.dir, dataset, project)
}

generatePlotsPG("mc_tm")



