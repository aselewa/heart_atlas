library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(GenomicRanges)
library(dplyr)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

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

custom_archr_plot <- function(archr_project, reduction = 'UMAP', group.by="regions", alpha=1, pt.size=1, palette=NULL, legend=TRUE, label=FALSE){
  
  cellColData <- as.data.frame(archr_project@cellColData)
  
  if(reduction == 'UMAP'){
    reduc.df <- as.data.frame(archr_project@embeddings$UMAP$df)
  }
  if(reduction == 'TSNE'){
    reduc.df <- as.data.frame(archr_project@embeddings$TSNE$df)
  }
  colnames(reduc.df) <- c("x1","x2")
  
  reduc.df[,"meta"] <- cellColData[,group.by, drop=T]
  p <- ggplot(reduc.df, aes(x=x1, y=x2, color=meta)) + ggrastr::rasterise(geom_point(size=pt.size, alpha=alpha), dpi = 200) + theme_set(theme_grey()) + ggClean() +
    xlab(paste0(reduction,'1')) + ylab(paste0(reduction,'2'))
  
  if(!is.null(palette)){
    p <- p + scale_color_manual(values = palette)
  }
  
  if(!legend){
    p <- p + LegendOff()
  }
  
  if(label){
    text.pos <- suppressMessages(reduc.df %>% group_by(meta) %>% summarise(pos_x = mean(x1), pos_y = mean(x2)))
    p <- p + ggrepel::geom_text_repel(data = text.pos, aes(x = pos_x, y = pos_y, label = meta), color="black")
  }
  
  return(p)
}


custom_dim_plot <- function(seurat, reduction = 'UMAP', group.by="cellTypes", alpha=1, pt.size=1, palette=NULL, legend=TRUE, label=FALSE){
  
  cellColData <- as.data.frame(seurat@meta.data)
  
  if(reduction == 'UMAP'){
    reduc.df <- as.data.frame(seurat@reductions$umap@cell.embeddings)
  }
  colnames(reduc.df) <- c("x1","x2")
  
  reduc.df[,"meta"] <- cellColData[,group.by, drop=T]
  p <- ggplot(reduc.df, aes(x=x1, y=x2, color=meta)) + ggrastr::rasterise(geom_point(size=pt.size, alpha=alpha), dpi = 200) + theme_set(theme_grey()) + ggClean() +
    xlab(paste0(reduction,'1')) + ylab(paste0(reduction,'2'))
  
  if(!is.null(palette)){
    p <- p + scale_color_manual(values = palette)
  }
  
  if(!legend){
    p <- p + LegendOff()
  }
  
  if(label){
    text.pos <- suppressMessages(reduc.df %>% group_by(meta) %>% summarise(pos_x = mean(x1), pos_y = mean(x2)))
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

snpIDtoGR <- function(id){
  chr <- sub(pattern = '_.*', replacement = "", x = id)
  tmp <- sub(pattern = 'chr[0-9]+_', replacement = "", x = id)
  pos <- sub(pattern = '_.+', replacement = "", x = tmp)
  
  id.gr <- StringToGR(paste0(chr, ':',pos,'-',pos))
  
  return(id.gr)
}

removeOverlaps <- function(X, to.remove){
  hits <- GenomicRanges::findOverlaps(query = X, subject = to.remove)
  if(length(hits) > 0){
    X <- X[-queryHits(hits),]
  }
  X
}

subsetByOverlapProp <- function(q, s, minPoverlap){
  
  hits <- GenomicRanges::findOverlaps(query = q, subject = s)
  overlaps <- pintersect(q[queryHits(hits)], s[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(q[queryHits(hits)])
  hits <- hits[percentOverlap >= minPoverlap]
  
  return(hits)
}

SNPGenomeDistrib <- function(snp.gr, genomic.annots) {
  
  annot <- plyranges::join_overlap_inner(snp.gr, genomic.annots)
  
  annot.freq <- as.data.frame(annot@elementMetadata) %>% 
    group_by(SNP) %>% 
    dplyr::count(type) %>% 
    mutate(n_bin = 1*(n>0)) %>% 
    group_by(SNP) %>% 
    mutate(f = n_bin / sum(n_bin))
  
  snpsIn <- length(unique(annot$SNP))
  nIntergenic <- length(snp.gr) - snpsIn
  snp.dist <- c("intergenic"=nIntergenic)
  snp.dist["exon"] <- sum(annot.freq$f[annot.freq$type=="exon"])
  snp.dist["UTR"] <- sum(annot.freq$f[annot.freq$type=="UTR"])
  snp.dist["intronic"] <- sum(annot.freq$f[annot.freq$type=="intron"])
  snp.dist["promoter"] <- sum(annot.freq$f[annot.freq$type=="promoter"])
  
  snp.dist.df <- data.frame(freq=100*snp.dist/length(snp.gr), category=names(snp.dist))
  snp.dist.df$category <- unfactor(snp.dist.df$category)
  
  snp.dist.df
  
}

join_overlap_list <- function(gr.list, X){
  res.list <- list()
  for(n in names(gr.list)){
    res.list[[n]] <- plyranges::join_overlap_inner(X, gr.list[[n]])
  }
  return(res.list)
}


get_elements_overlap_snps <- function(snp.gr, annotations) {
    for (f in annotations) {
        name <- paste0(basename(f), "_d")
        curr <- rtracklayer::import(f, format = "bed")
        curr <- GenomicRanges::reduce(curr)
        overlap.df <- plyranges::join_overlap_inner(curr, snp.gr) %>% 
            as_tibble() %>% mutate(enhancer = paste0('chr',seqnames, ':', start, '-', end)) %>% 
            dplyr::select(snp, enhancer)
        colnames(overlap.df) <- c("snp", sub('.bed','',basename(annotations)))
    }
    return(overlap.df)
}

get_snp_to_gene_dist <- function(snp.gr, gene.gr, genes=NULL){
    
    if(!is.null(genes)){
        gene.gr <- gene.gr[gene.gr$gene_name %in% genes,]
        
    }
    
}


