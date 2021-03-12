results.dir <- "/ahg/regevdata/projects/scqc/output_pg/"
fout <- "figure3_table.csv"


HEADER <- "Directory,Tissue,n_cells after basic QC,n_cells after ddqc,pct,n_cells after cutoff,pct,n_cells retained by both methods,pct,pct joint,n_cells retained only by ddqc,pct,pct joint,n_cells retained only by cutoff,pct,pct joint,n_cells that failed ddqc,pct,n_cells that failed ddqc due to n_counts,pct,n_cells that failed ddqc due to n_genes,pct,n_cells that failed ddqc due to percent.mito,pct,n_cells that failed ddqc due to n_genes and percent_mito,pct,n_cells that failed ddqc and were retained by cutoff,pct,n_cells that failed ddqc due to n_counts and were retained by cutoff,pct,n_cells that failed ddqc due to n_genes and were retained by cutoff,pct,n_cells that failed ddqc due to percent.mito and were retained by cutoff,pct,n_cells that failed ddqc due to n_genes and percent_mito and were retained by cutoff,pct"
write(HEADER, fout)

for (project in dir(results.dir, recursive=FALSE)) {
  for (tissue in dir(paste0(results.dir, project, "/"), recursive=FALSE)) {
    path <- paste0(results.dir, project, "/", tissue, "/")
    
    if (!file.exists(paste0(path,  "1.4-joint_clustering_old/!cells.csv")) || !file.exists(paste0(path, "1.4-mad-2/!filtered_counts.csv"))) {
      print(path)
      next
    }
    all.cells <- read.csv(paste0(path, "1.4-none-0/!cells.csv"))$barcodekey
    cutoff.cells <- read.csv(paste0(path, "1.4-cutoff-10/!cells.csv"))$barcodekey
    mad.cells <- read.csv(paste0(path, "1.4-mad-2/!cells.csv"))$barcodekey
    joint.cells <- read.csv(paste0(path, "1.4-joint_clustering_old/!cells.csv"))$barcodekey
    joint.cells.color <- read.csv(paste0(path, "1.4-joint_clustering_old/!cells.csv"))$color
    
    mad.failed.counts <- read.csv(paste0(path, "1.4-mad-2/!filtered_counts.csv"))$barcodekey
    mad.failed.genes <- read.csv(paste0(path, "1.4-mad-2/!filtered_genes.csv"))$barcodekey
    mad.failed.mito <- read.csv(paste0(path, "1.4-mad-2/!filtered_mito.csv"))$barcodekey
    
    n.cells.basic.qc <- length(all.cells)
    n.cells.ddqc <- length(mad.cells)
    n.cells.cutoff <- length(cutoff.cells)
    
    n.cells.both <- table(joint.cells.color)["All"]
    n.cells.mad.only <- table(joint.cells.color)["MAD2 only"]
    n.cells.cutoff.only <- table(joint.cells.color)["Cutoff only"]
    
    n.cells.fail.ddqc <- length(union(union(mad.failed.genes, mad.failed.mito), mad.failed.counts))
    n.cells.fail.ddqc.counts <- length(mad.failed.counts)
    n.cells.fail.ddqc.genes <- length(mad.failed.genes)
    n.cells.fail.ddqc.mito <- length(mad.failed.mito)
    n.cells.fail.ddqc.genes_and_mito <- length(intersect(mad.failed.genes, mad.failed.mito))
    
    n.cells.cutoff.exclusive <- length(intersect(union(union(mad.failed.genes, mad.failed.mito), mad.failed.counts), cutoff.cells))
    n.cells.fail.cutoff.exclusive.counts <- length(intersect(cutoff.cells, mad.failed.counts))
    n.cells.fail.cutoff.exclusive.genes <- length(intersect(cutoff.cells, mad.failed.genes))
    n.cells.fail.cutoff.exclusive.mito <- length(intersect(cutoff.cells, mad.failed.mito)) 
    n.cells.fail.cutoff.exclusive.genes_and_mito <- length(intersect(cutoff.cells, intersect(mad.failed.genes, mad.failed.mito)))
    
    write(paste(path, tissue, n.cells.basic.qc, 
                n.cells.ddqc, round(n.cells.ddqc/n.cells.basic.qc, 3)*100, 
                n.cells.cutoff, round(n.cells.cutoff/n.cells.basic.qc, 3)*100,
                n.cells.both, round(n.cells.both/length(joint.cells), 3)*100, round(n.cells.both/n.cells.basic.qc, 3)*100,
                n.cells.mad.only, round(n.cells.mad.only/length(joint.cells), 3)*100, round(n.cells.mad.only/n.cells.basic.qc, 3)*100,
                n.cells.cutoff.only, round(n.cells.cutoff.only/length(joint.cells), 3)*100, round(n.cells.cutoff.only/n.cells.basic.qc, 3)*100,
                n.cells.fail.ddqc, round(n.cells.fail.ddqc/n.cells.basic.qc, 3)*100, 
                n.cells.fail.ddqc.counts, round(n.cells.fail.ddqc.counts/n.cells.fail.ddqc, 3)*100,
                n.cells.fail.ddqc.genes, round(n.cells.fail.ddqc.genes/n.cells.fail.ddqc, 3)*100,
                n.cells.fail.ddqc.mito, round(n.cells.fail.ddqc.mito/n.cells.fail.ddqc, 3)*100,
                n.cells.fail.ddqc.genes_and_mito, round(n.cells.fail.ddqc.genes_and_mito/n.cells.fail.ddqc, 3)*100,
                n.cells.cutoff.exclusive, round(n.cells.cutoff.exclusive/n.cells.basic.qc, 3)*100, 
                n.cells.fail.cutoff.exclusive.counts, round(n.cells.fail.cutoff.exclusive.counts/n.cells.cutoff.exclusive, 3)*100, 
                n.cells.fail.cutoff.exclusive.genes, round(n.cells.fail.cutoff.exclusive.genes/n.cells.cutoff.exclusive, 3)*100, 
                n.cells.fail.cutoff.exclusive.mito, round(n.cells.fail.cutoff.exclusive.mito/n.cells.cutoff.exclusive, 3)*100,
                n.cells.fail.cutoff.exclusive.genes_and_mito, round(n.cells.fail.cutoff.exclusive.genes_and_mito/n.cells.cutoff.exclusive, 3)*100,
                sep=","), file=fout, append = TRUE)
  }
}