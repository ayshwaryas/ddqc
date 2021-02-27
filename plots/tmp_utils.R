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
  for (cl in unique(obj$louvain_labels)) {
    cluster <- subset(obj, louvain_labels == cl)
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
      percents2 <- c(percents2, second.frequent / length(cluster$louvain_labels))
    }
    percents1 <- c(percents1, first.frequent / length(cluster$louvain_labels))
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
    obj.cluster <- subset(obj, louvain_labels == cl)
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
      cells <- c(cells, length(obj.cluster$louvain_labels))
      genes.mean <- c(genes.mean, mean(obj.cluster$n_genes))
      genes.median <- c(genes.median, median(obj.cluster$n_genes))
      mito.mean <- c(mito.mean, mean(obj.cluster$percent_mito))
      mito.median <- c(mito.median, median(obj.cluster$percent_mito))
      
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




