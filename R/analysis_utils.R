library(ggplot2)
library(patchwork)
library(RColorBrewer)
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
  ggplot(df, aes_string(x=x, y=y, fill=fill)) + geom_violin() + geom_jitter(shape=16, size=0.5, position=position_jitter(0.4)) + ggClean() + xlab("") + LegendOff()
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
    ggClean() + scale_fill_manual(values = palette) + xlab("") + ylab("Prop. Cells")
}

RenameIdentity <- function(idents, from, to){
  new.idents <- plyr::mapvalues(idents, from = from, to = to)
  return(new.idents)
}


