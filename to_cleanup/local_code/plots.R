ggsave1 <- function(filename, plot, n.tissues=30) {
  n.tissues <- max(10, n.tissues)
  no_bkg <- theme(axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank()) 
  ggsave(filename = filename, plot = plot + no_bkg, width = 14 / 30 * n.tissues, height = 10) #saves plot with custom dimensions 
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
      hl <- geom_hline(yintercept=10, color="red", size=1)
      vl <- geom_vline(xintercept=10, color="red", size=1)
    }
    else {
      if (metric.name.seurat == "nFeature_RNA") {
        hl <- geom_hline(yintercept=log2(200), color="red", size=1)
        vl <- geom_vline(xintercept=log2(200), color="red", size=1)
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
    data$tissues = with(data, reorder(tissues, -metric, mean)) #order data by tissue median
    
    #all plots except combined density and tsne/umap are saved in to vesions: no trasformation and log2, if save.log2 is true
    
    name.prefix <- paste0(results.dir, add.dir)
    
    if (metric.name.seurat == "nCount_RNA" || metric.name.seurat == "nFeature_RNA") {
    
      #boxplot by tissue
      ggsave1(filename=paste0(name.prefix, "box_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=tissues, y=log2(metric))) + geom_boxplot() + t2 + t3 + t6 + c1 + c2 + l2 + hl, n.tissues = n.tissues) 
      #joyplot by tissue
      ggsave1(filename=paste0(name.prefix, "density2_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=log2(metric), y=tissues)) + geom_density_ridges() + t3 + t7 + c1 + c2 + l4 + vl, n.tissues = n.tissues) 
      #violin plot by tissue
      ggsave1(filename=paste0(name.prefix, "violin_", name.suffix, "_log.pdf"), plot=ggplot(subset(data, metric > 0), aes(x=tissues, y=log2(metric))) + geom_violin() + t2 + t3 + t6 + c1 + c2 + l2 + hl, n.tissues = n.tissues)
    } else {
      #boxplot by tissue
      ggsave1(filename=paste0(name.prefix, "box_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=tissues, y=metric)) + geom_boxplot() + t2 + t3 + t6 + c1 + c2 + l1 + a1 + hl, n.tissues = n.tissues) 
      #joyplot by tissue
      ggsave1(filename=paste0(name.prefix, "density2_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=metric, y=tissues)) + geom_density_ridges() + t3 + t7 + c1 + c2 + l3 + a2 + vl, n.tissues = n.tissues) 
      #violin plot by tissue
      ggsave1(filename=paste0(name.prefix, "violin_", name.suffix, ".pdf"), plot=ggplot(data, aes(x=tissues, y=metric)) + geom_violin() + t2 + t3 + t6 + c1 + c2 + l1 + a1 + hl, n.tissues = n.tissues)
    }
  }
}

generatePlotsPG <- function(project) {
  data.from.pg <<- TRUE
  source("../scripts/settings.R")
  source("../scripts/local_settings.R")
  data.path <- paste0(source.dir.prefix, project, "/")
  dataset <- NULL
  for (directory in list.files(data.path)) {
    d <- read.csv(paste0(data.path, directory, "/1.4-none-0/!cells.csv"))
    d$tissue <- directory
    if (is.null(dataset)) {
      dataset <- d
    }
    else {
      dataset <- full_join(dataset, d)
    }
  }
  dataset <- dataset %>% rename(nFeature_RNA = n_genes, nCount_RNA = n_counts, percent.mt = percent_mito, percent.rb = percent_ribo, seurat_clusters = louvain_labels)
  
  results.dir <- paste0(source.dir.prefix, "summary_plots/", project, "/")
  dir.create(paste0(source.dir.prefix, "summary_plots/"), showWarnings = FALSE)
  dir.create(paste0(results.dir), showWarnings = FALSE)
  write.csv(dataset, paste0(results.dir, "!cells.csv"))
  
  generatePlots(results.dir, dataset, project)
}

generatePlotsSeurat <- function(project) {
  source("../scripts/settings.R")
  source("../scripts/local_settings.R")
  data.path <- paste0(source.dir.prefix, project, "/")
  dataset <- NULL
  for (directory in list.files(data.path)) {
    if (grepl(".csv", directory)) {
      next
    }
    d <- read.csv(paste0(data.path, directory, "/1-none-0/!cells.csv"))
    d$tissue <- directory
    if (is.null(dataset)) {
      dataset <- d
    }
    else {
      dataset <- full_join(dataset, d)
    }
  }
  # dataset <- dataset %>% rename(nFeature_RNA = n_genes, nCount_RNA = n_counts, percent.mt = percent_mito, percent.rb = percent_ribo, seurat_clusters = louvain_labels)
  
  results.dir <- paste0("~/Documents/primes_storage/output_pg/", "summary_plots/", project, "/")
  # dir.create(paste0(source.dir.prefix, "summary_plots/"), showWarnings = FALSE)
  dir.create(paste0(results.dir), showWarnings = FALSE)
  write.csv(dataset, paste0(results.dir, "!cells.csv"))
  
  generatePlots(results.dir, dataset, project)
}

generatePlotsInit <- function() {
  for (project in c("human", "mc_ebi", "mc_ebi_tm", "mc_mca", "mc_tm", "mc_ts30", "mc_PanglaoDB")) {
    data.path <- paste0(source.dir, project, "/stats_summary.csv")
    results.dir <- paste0(source.dir, "summary_plots/", project, "/")
    dir.create(paste0(source.dir, "summary_plots/"), showWarnings = FALSE)
    dir.create(paste0(results.dir), showWarnings = FALSE)
    
    if (project == "human") {
      d1 <- as.data.frame(read.csv(paste0(source.dir, "/mc_other/stats_summary.csv"), header = FALSE))
      d2 <- as.data.frame(read.csv(paste0(source.dir, "/mc_other_10X/stats_summary.csv"), header = FALSE))
      colnames(d1) <- c("#", "tissue", "cell", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
      colnames(d2) <- c("#", "tissue", "cell", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
      dataset <- full_join(d1, d2)
      dataset <- subset(dataset, !grepl("mouse", tissue))
      write.csv(dataset, "human.csv")
      return()
      
    } else {
      dataset <- as.data.frame(read.csv(data.path, header = FALSE))
      colnames(dataset) <- c("#", "tissue", "cell", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
    }
    generatePlots(results.dir, dataset, project)
  }
}

