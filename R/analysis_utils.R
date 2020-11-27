library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(GenomicRanges)
library(dplyr)

ggClean <- function(rotate_axis=FALSE){
  tm <- theme_bw() + 
    theme(text = element_text(size=18),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1))
  if(rotate_axis){
    tm <- tm + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  
  return(tm)
  
}

LegendOff <- function(){
  theme(legend.position = "none")
}

qcVlnPlot <- function(df, x,y,fill){
  ggplot(df, aes_string(x=x, y=y, fill=fill)) + geom_violin() + geom_jitter(shape=16, size=0.5, position=position_jitter(0.4)) + ggClean(rotate_axis=T) + xlab("") + LegendOff()
}

pullGeneScoreMatrix <- function(archr_project){
  gene.score <- getMatrixFromProject(ArchRProj = archr_project, useMatrix = 'GeneScoreMatrix')
  gene.score.mat <- as.matrix(assays(gene.score)[[1]])
  row.names(gene.score.mat) <- rowData(gene.score)$name
  gene.score.mat.norm <- imputeMatrix(gene.score.mat, imputeWeights = getImputeWeights(archr_project))
  return(gene.score.mat.norm)
}

custom_archr_umap <- function(archr_project, group.by="regions", alpha=1, pt.size=1, palette=NULL, legend=TRUE, label=FALSE){
  
  cellColData <- as.data.frame(archr_project@cellColData)
  umap.df <- as.data.frame(archr_project@embeddings$UMAP$df)
  colnames(umap.df) <- c("UMAP_1","UMAP_2")
  
  umap.df[,"meta"] <- cellColData[,group.by, drop=T]
  p <- ggplot(umap.df, aes(x=UMAP_1, y=UMAP_2, color=meta)) + geom_point(size=pt.size, alpha=alpha) + theme_set(theme_grey()) + ggClean()
  
  if(!is.null(palette)){
    p <- p + scale_color_manual(values = palette)
  }
  
  if(!legend){
    p <- p + LegendOff()
  }
  
  if(label){
    text.pos <- suppressMessages(umap.df %>% group_by(meta) %>% summarise(pos_x = mean(UMAP_1), pos_y = mean(UMAP_2)))
    p <- p + ggrepel::geom_text_repel(data = text.pos, aes(x = pos_x, y = pos_y, label = meta), color="black")
  }
  
  return(p)
}

getQCPlots <- function(archr_project){
  cellColData <- as.data.frame(archr_project@cellColData)
  p1 <- qcVlnPlot(df = cellColData, x = "regions", y = "nFrags", fill = "regions")
  p2 <- qcVlnPlot(df = cellColData, x = "regions", y = "BlacklistRatio", fill = "regions")
  p3 <- qcVlnPlot(df = cellColData, x = "regions", y = "TSSEnrichment", fill="regions")
  p4 <- qcVlnPlot(df = cellColData, x = "regions", y = "NucleosomeRatio", fill="regions")
  p <- p1 + p2 + p3 + p4 + plot_layout(nrow=1)
}

make_freq_plot <- function(freq, palette){
  type.freq <- data.frame(freq=as.numeric(freq), celltype=factor(names(freq)))
  p <- ggplot(type.freq, aes(x=celltype, y=freq, fill=celltype)) + geom_bar(position = "dodge", stat="identity", color="black") + 
    ggClean(rotate_axis = T) + scale_fill_manual(values = palette) + xlab("") + ylab("Prop. Cells")
}

RenameIdentity <- function(idents, from, to){
  new.idents <- plyr::mapvalues(idents, from = from, to = to)
  return(new.idents)
}

GRToString <- function(gr){
  paste0(as.character(seqnames(gr)),':',start(gr),'-',end(gr))
}

StringToGR <- function(st){
  chr <- sub(pattern = ':.*', replacement = "", st)
  st.end <- sub(pattern = '.*:', replacement = "", st)
  s <- as.numeric(sub(pattern = '-.*', replacement = "", st.end))
  e <- as.numeric(sub(pattern = '.*-', replacement = "", st.end))
  gr <- GenomicRanges::GRanges(seqnames = chr, IRanges::IRanges(start = s, end = e))
  return(gr)
}

qOverlapS <- function(q, s, minPoverlap){
  
  hits <- GenomicRanges::findOverlaps(query = q, subject = s)
  overlaps <- pintersect(q[queryHits(hits)], s[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(q[queryHits(hits)])
  hits <- hits[percentOverlap > minPoverlap]
  
  inQuery <- q[queryHits(hits)]
  propIn <- length(unique(GRToString(inQuery)))/length(q)
  return(propIn)
}


getSharingMat <- function(gr.list){
  l <- length(gr.list)
  m <- matrix(0, nrow=l, ncol=l)
  for(i in 1:l){
    curr_sum <- 0
    for(j in 1:l){
      m[i,j] <- qOverlapS(q = gr.list[[i]], s = gr.list[[j]], 0.5)
    }
  }
  m
}

generateBits <- function(n){
  max_n <- 2^n - 1
  bitList <- list()
  for(i in 1:max_n){
    bitList[[i]] <- as.integer(intToBits(i)[1:n])
  }
  return(bitList)
}

getIdealAccess <- function(cell_type_vec){
  
  cellTypes <- unique(cell_type_vec)
  nCellTypes <- length(cellTypes)
  profiles <- generateBits(nCellTypes)
  nProfiles <- length(profiles)
  
  ideal.mat <- matrix(0, nrow = nProfiles, ncol = length(cell_type_vec))

  for(i in 1:nProfiles){
    curr.profile <- profiles[[i]]
    if(sum(curr.profile) >= 3)
    toInclude <- which(curr.profile == 1)
    cellTypesIn <- cellTypes[toInclude]
    for(type in cellTypesIn){
      ideal.mat[i,] <- ideal.mat[i,] + 1*(cell_type_vec==type)
    }
  }
  return(ideal.mat)
  
}

multi.jaccard.index <- function(X, y){
  stopifnot(is.matrix(X))
  stopifnot(dim(X)[1] == length(y))
  intsct <- X * y
  yunion <- X | y
  colSums(intsct) / colSums(yunion)
}

peakToCluster <- function(peak.mat, ideal.mat){
  
  classify.mat <- matrix(0, nrow = nrow(peak.mat), ncol = nrow(ideal.mat))
  ideal.mat.t <- t(ideal.mat)
  
  for(i in 1:nrow(peak.mat)){
    if(i %% 1000 == 0){
      print(i)
    }
    x <- peak.mat[i,,drop=T]
    dists <- multi.jaccard.index(ideal.mat.t, x)
    classify.mat[i,] <- dists
  }
  return(classify.mat)
}


peakToClusterBatch <- function(peak.mat, ideal.mat, chunk.size=1e5){
  N <- nrow(peak.mat)
  k <- floor(N/1e5)
  classify.mat.list <- list()
  for(i in 1:k){
    print(paste0('Getting batch ',i,'...'))
    curr.mat <- as.matrix(peak.mat[(1+(i-1)*chunk.size):(chunk.size*i),])
    classify.mat.list[[i]] <- peakToCluster(peak.mat = curr.mat, ideal.mat = ideal.mat)
    rm(curr.mat)
    gc()
  }
  if((k*chunk.size) < N){
    curr.mat <- as.matrix(peak.mat[(1+k*chunk.size):N,])
    classify.mat.list[[i+1]] <- peakToCluster(peak.mat = curr.mat, ideal.mat = ideal.mat)
    rm(curr.mat)
    gc()
  }
  classify.mat <- Reduce(rbind, classify.mat.list)
  return(classify.mat)
}

