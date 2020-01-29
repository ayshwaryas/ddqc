ggsave1 <- function(filename, plot) {
  no_bkg <- theme(axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank()) 
  ggsave(filename = filename, plot = plot + no_bkg, width = 14, height = 10) #saves plot with custom dimensions 
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
    t5 <- facet_wrap(. ~ tissues, ncol=5)
    t6 <- stat_summary(fun.y=mean, geom="point", shape=23, fill="blue", size=3)
    t7 <- theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10, face="bold"), legend.position="none", axis.title.y=element_blank())
    l1 <- labs(y=metric.name)
    l2 <- labs(y=paste0("log2(", metric.name, ")"))
    l3 <- labs(x=metric.name)
    l4 <- labs(x=paste0("log2(", metric.name, ")"))
    l5 <- labs(y=paste0("Average ", metric.name))
    l6 <- labs(y=paste0("Average log2(", metric.name, ")"))
    if (metric.name.seurat == "percent.mt") {
      a1 <- scale_y_continuous(breaks=seq(0, 100, 5))
      a2 <- scale_x_continuous(breaks=seq(0, 100, 5))
    }
    else {
      a1 <- NULL
      a2 <- NULL
    }
    
    data <- data.frame(metric=dataset[[metric.name.seurat]], tissues=dataset$tissue) 
    colnames(data) <- c("metric", "tissues") #rename data columns
    data$tissues = with(data, reorder(tissues, -metric, median)) #order data by tissue median
    
    #all plots except combined density and tsne/umap are saved in to vesions: no trasformation and log2, if save.log2 is true
    
    #barplot of tissue means
    ggsave1(filename=paste0(results.dir, add.dir, paste0("bar_mean_", name.suffix, "_log.pdf")), plot=ggplot(subset(data, metric > 0) %>% group_by(tissues) %>% summarize(metric = mean(log2(metric))), aes(x=tissues, y=log2(metric))) + geom_bar(aes(fill=tissues), stat="identity") + t2 + t3 + l6)
    #boxplot by tissue
    ggsave1(filename=paste0(results.dir, add.dir, paste0("box_", name.suffix, "_log.pdf")), plot=ggplot(subset(data, metric > 0), aes(x=tissues, y=log2(metric))) + geom_boxplot(aes(fill=tissues)) + t2 + t3 + t6 + l2) 
    #combined density plots for each tissue
    ggsave1(filename=paste0(results.dir, add.dir, paste0("density_", name.suffix, ".pdf")), plot=ggplot(subset(data, metric > 0), aes(x=log2(metric))) + geom_density(aes(fill=tissues)) + t3 + t4 + t5 + l4)
    #joyplot by tissue
    ggsave1(filename=paste0(results.dir, add.dir, paste0("density2_", name.suffix, "_log.pdf")), plot=ggplot(subset(data, metric > 0), aes(x=log2(metric), y=tissues)) + geom_density_ridges(aes(fill=tissues)) + t3 + t7 + l4) 
    #violin plot by tissue
    ggsave1(filename=paste0(results.dir, add.dir, paste0("violin_", name.suffix, "_log.pdf")), plot=ggplot(subset(data, metric > 0), aes(x=tissues, y=log2(metric))) + geom_violin(aes(fill=tissues)) + t2 + t3 + t6 + l2)

    
    #barplot of tissue means
    ggsave1(filename=paste0(results.dir, add.dir, paste0("bar_mean_", name.suffix, ".pdf")), plot=ggplot(data %>% group_by(tissues) %>% summarize(metric = mean(metric)), aes(x=tissues, y=metric)) + geom_bar(aes(fill=tissues), stat="identity") + t2 + t3 + l5 + a1)
    #boxplot by tissue
    ggsave1(filename=paste0(results.dir, add.dir, paste0("box_", name.suffix, ".pdf")), plot=ggplot(data, aes(x=tissues, y=metric)) + geom_boxplot(aes(fill=tissues)) + t2 + t3 + t6 + l1 + a1) 
    #joyplot by tissue
    ggsave1(filename=paste0(results.dir, add.dir, paste0("density2_", name.suffix, ".pdf")), plot=ggplot(data, aes(x=metric, y=tissues)) + geom_density_ridges(aes(fill=tissues)) + t3 + t7 + l3 + a2) 
    #violin plot by tissue
    ggsave1(filename=paste0(results.dir, add.dir, paste0("violin_", name.suffix, ".pdf")), plot=ggplot(data, aes(x=tissues, y=metric)) + geom_violin(aes(fill=tissues)) + t2 + t3 + t6 + l1 + a1)
  }
}

generatePlotsInit <- function() {
  for (project in c("human", "mc_ebi", "mc_ebi_tm", "mc_mca", "mc_tm", "mc_ts30")) {
    for(mito.filter in  c(80, 100)) {
      data.path <- paste0(source.dir, project, "/stats_summary.csv")
      results.dir <- paste0(source.dir, "summary_plots/", project, "/mito-", mito.filter, "/")
      dir.create(paste0(source.dir, "summary_plots/"), showWarnings = FALSE)
      dir.create(paste0(source.dir, "summary_plots/", project), showWarnings = FALSE)
      dir.create(results.dir, showWarnings = FALSE)
      
      if (project == "human") {
        d1 <- as.data.frame(read.csv(paste0(source.dir, "../ss_other.csv"), header = FALSE))
        d2 <- as.data.frame(read.csv(paste0(source.dir, "../ss_other_10X.csv"), header = FALSE))
        colnames(d1) <- c("#", "tissue", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
        colnames(d2) <- c("#", "tissue", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
        dataset <- full_join(d1, d2)
        dataset <- subset(dataset, !grepl("mouse", tissue))
        
      }
      else {
        dataset <- as.data.frame(read.csv(data.path, header = FALSE))
        colnames(dataset) <- c("#", "tissue", "cluster", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
      }
      dataset <- subset(dataset, percent.mt < mito.filter)
      generatePlots(results.dir, dataset, paste0(project, "-", mito.filter))
    }
  }
}

