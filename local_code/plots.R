source("~/Documents/primes/code/local_code/local_settings.R")

ggsave1 <- function(filename, plot) {
  ggsave(filename = filename, plot = plot, width = 12, height = 8) #saves plot with custom dimensions 
}

generatePlots <- function(results.dir, dataset, project) {
  for (metric in c("nCount_RNA count", "nFeature_RNA genes", "percent.mt mito", "percent.rb ribo")) {
    metric.name <- strsplit(metric, " ")[[1]][1]
    name.suffix <- strsplit(metric, " ")[[1]][2]
    
    #themes and axis labels for plots
    t2 <- theme(axis.text.x = element_text(angle = 45, size=10, hjust=1, face="bold"), axis.text.y = element_text(size=10), legend.position="none", axis.title.x=element_blank())
    t3 <- ggtitle(project)
    t4 <- theme(legend.position="none")
    t5 <- facet_wrap(. ~ tissues, ncol=5)
    t6 <- stat_summary(fun.y=mean, geom="point", shape=23, fill="blue", size=3)
    t7 <- theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10, face="bold"), legend.position="none", axis.title.y=element_blank())
    l1 <- labs(y=metric.name)
    l2 <- labs(y=paste0("log2(", metric.name, ")"))
    l3 <- labs(x=metric.name)
    l4 <- labs(x=paste0("log2(", metric.name, ")"))
    
    data <- data.frame(metric=dataset[[metric.name]], tissues=dataset$tissue) 
    colnames(data) <- c("metric", "tissues") #rename data columns
    data$tissues = with(data, reorder(tissues, -metric, mean)) #order data by tissue mean
    
    #all plots except combined density and tsne/umap are saved in to vesions: no trasformation and log2
    
    #barplot of tissue means
    ggsave1(filename=paste0(results.dir, paste0("bar_", name.suffix, ".pdf")), plot=ggplot(data %>% group_by(tissues) %>% summarize(metric = mean(metric)), aes(x=tissues, y=metric)) + geom_bar(aes(fill=tissues), stat="identity") + t2 + t3 + l1)
    ggsave1(filename=paste0(results.dir, paste0("bar_", name.suffix, "_log.pdf")), plot=ggplot(subset(data, metric > 0) %>% group_by(tissues) %>% summarize(metric = mean(metric)), aes(x=tissues, y=log2(metric))) + geom_bar(aes(fill=tissues), stat="identity") + t2 + t3 + l2)
    
    #boxplot by tissue
    ggsave1(filename=paste0(results.dir, paste0("box_", name.suffix, ".pdf")), plot=ggplot(data, aes(x=tissues, y=metric)) + geom_boxplot(aes(fill=tissues)) + t2 + t3 + t6 + l1) 
    ggsave1(filename=paste0(results.dir, paste0("box_", name.suffix, "_log.pdf")), plot=ggplot(subset(data, metric > 0), aes(x=tissues, y=log2(metric))) + geom_boxplot(aes(fill=tissues)) + t2 + t3 + t6 + l2) 
    
    #combined density plots for each tissue
    ggsave1(filename=paste0(results.dir, paste0("density_", name.suffix, ".pdf")), plot=ggplot(subset(data, metric > 0), aes(x=log2(metric))) + geom_density(aes(fill=tissues)) + t3 + t4 + t5 + l4)
    
    #joyplot by tissue
    ggsave1(filename=paste0(results.dir, paste0("density2_", name.suffix, ".pdf")), plot=ggplot(data, aes(x=metric, y=tissues)) + geom_density_ridges(aes(fill=tissues)) + t3 + t7 + l3) 
    ggsave1(filename=paste0(results.dir, paste0("density2_", name.suffix, "_log.pdf")), plot=ggplot(subset(data, metric > 0), aes(x=log2(metric), y=tissues)) + geom_density_ridges(aes(fill=tissues)) + t3 + t7 + l4) 
    
    #overall destiny plot
    ggsave1(filename=paste0(results.dir, paste0("density3_", name.suffix, ".pdf")), plot=ggplot(data, aes(x=metric)) + geom_density(aes(fill="red")) + t3 + l3) 
    ggsave1(filename=paste0(results.dir, paste0("density3_", name.suffix, "_log.pdf")), plot=ggplot(subset(data, metric > 0), aes(x=log2(metric))) + geom_density(aes(fill="red")) + t3 + l4) 
    
    #violin plot by tissue
    ggsave1(filename=paste0(results.dir, paste0("violin_", name.suffix, ".pdf")), plot=ggplot(data, aes(x=tissues, y=metric)) + geom_violin(aes(fill=tissues)) + t2 + t3 + t6 + l1)
    ggsave1(filename=paste0(results.dir, paste0("violin_", name.suffix, "_log.pdf")), plot=ggplot(subset(data, metric > 0), aes(x=tissues, y=log2(metric))) + geom_violin(aes(fill=tissues)) + t2 + t3 + t6 + l2)
  }
}

for (project in c("mc_ebi", "mc_mca", "mc_tm", "mc_ts30")) {
  for(mito.filter in  c(80, 100)) {
    data.path <- paste0(source.dir, project, "/stats_summary.csv")
    results.dir <- paste0(source.dir, "summary_plots/", project, "/mito-", mito.filter, "/")
    dir.create(paste0(source.dir, "summary_plots/"), showWarnings = FALSE)
    dir.create(paste0(source.dir, "summary_plots/", project), showWarnings = FALSE)
    dir.create(results.dir, showWarnings = FALSE)
    
    dataset <- as.data.frame(read.csv(data.path, header = FALSE))
    colnames(dataset) <- c("#", "tissue", "cluster", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
    dataset <- subset(dataset, percent.mt < mito.filter)
    generatePlots(results.dir, dataset, paste0(project, "-", mito.filter))
  }
}

