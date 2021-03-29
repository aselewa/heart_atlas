library(ggplot2)
library(Gviz)
library(lattice)
source('R/analysis_utils.R')

require(org.Hs.eg.db)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(cowplot)

## below are functios for plotting each component of the track plot
## Note that the code assumes column names that I specified in my dataframes
## you should change these column names accordingly

# Needs a position, p-value, and r2 (LD) column
plot.pval.track.y.axis <- function(pip.df, curr.locus.gr){
  
  pvl.plot <- ggplot(pip.df, aes(x=pos, y=pval, color=r2)) + 
    geom_point() + 
    theme_classic() + 
    LegendOff() + 
    ylab("-log10 pvalue") + 
    xlab("") + 
    scale_color_manual(values = c("black","darkblue","darkgreen","green","red")) + 
    xlim(c(start(curr.locus.gr), end(curr.locus.gr))) +
    theme(axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank())
  
  pvl.plot
}

plot.pval.track <- function(pip.df, curr.locus.gr){
    
    pvl.plot <- ggplot(pip.df, aes(x=pos, y=pval, color=r2)) + 
        geom_point() + 
        theme_void() + 
        LegendOff() + 
        scale_color_manual(values = c("black","darkblue","darkgreen","green","red")) + 
        xlim(c(start(curr.locus.gr), end(curr.locus.gr))) 
    pvl.plot
}

# Needs a position, PIP from fine-mapping, and a snp column with rsIDs
plot.pip.track.w.y.axis <- function(pip.df, curr.locus.gr){
  
  pip.plot <- ggplot(pip.df, aes(x=pos, y=susie_pip)) + 
    geom_point() + 
    ggrepel::geom_text_repel(data=subset(pip.df, isAnnotate==1), aes(label=snp), size=4, min.segment.length = 1.5) +
    theme_classic() + 
    LegendOff() + 
    ylab("PIP") + 
    xlab("") + 
    xlim(c(start(curr.locus.gr), end(curr.locus.gr))) +
    theme(axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank())
  
  pip.plot
}

plot.pip.track <- function(pip.df, curr.locus.gr){
    
    pip.plot <- ggplot(pip.df, aes(x=pos, y=susie_pip)) + 
        geom_point() + 
        ggrepel::geom_text_repel(data=subset(pip.df, isAnnotate==1), aes(label=snp), size=4, min.segment.length = 1.5) +
        theme_void() + 
        LegendOff() + 
        xlim(c(start(curr.locus.gr), end(curr.locus.gr)))
        
    pip.plot
}


# only requires the region of interest
knownGeneObject <- function(curr.locus.gr, genome){
  
  ga.track <- GenomeAxisTrack()
  if(genome == "hg19"){
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  }
  if(genome == "hg38"){
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  }
  
  grtrack <- GeneRegionTrack(range = txdb, 
                             genome = "hg19", 
                             chromosome = as.character(seqnames(curr.locus.gr)),
                             start = start(curr.locus.gr),
                             end = end(curr.locus.gr))
  
  symbol(grtrack) <- mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                            keys=sub("\\.\\d+$", "", gene(grtrack)), 
                            keytype="ENTREZID", column="SYMBOL")
  
  symbol(grtrack) <- ifelse(is.na(symbol(grtrack)), gene(grtrack), symbol(grtrack))
  return(grtrack)
}

geneTrackPlot <- function(curr.locus.gr, collapseTranscripts=T, genome="hg19"){
    
    knownGeneObject <- knownGeneObject(curr.locus.gr = curr.locus.gr, genome)
    
    gene.track <- lattice::xyplot(1 ~ 1 | 1, 
                                  strip = F, 
                                  panel = function(x){plotTracks(knownGeneObject, 
                                                                 panel.only = T, 
                                                                 transcriptAnnotation = "symbol", 
                                                                 collapseTranscripts= collapseTranscripts, 
                                                                 add = T, 
                                                                 chromosome = as.character(seqnames(curr.locus.gr)),
                                                                 from = start(curr.locus.gr), 
                                                                 to = end(curr.locus.gr))},
                                  xlab=NULL,
                                  ylab=NULL,
                                  scales = list(draw = FALSE),
                                  par.settings = list(axis.line = list(col = 0)))
    return(gene.track)
}

# needs a GRanges object with a score metadata column that contains the accessibility
plotATAC <- function(gr, region.focus){
  gr <- subsetByOverlaps(gr, region.focus)
  score <- gr$score
  plot.df <- data.frame(start = start(gr), end = end(gr), score = score)
  
  p <- ggplot(plot.df) + geom_rect(aes(xmin = start, xmax = end, ymin=0, ymax=score), color = "darkgreen", fill="darkgreen") + 
    theme_classic() + 
    ylab("") +
    xlab("") + 
    LegendOff() + 
    theme_void() +
    ylim(c(0, 150)) +
    coord_cartesian(xlim=c(start(region.focus), end(region.focus))) 
  
  return(p)
}

plotRect <- function(gr, region.focus, col){
    
    seqlevelsStyle(gr) <- "UCSC"
    seqlevelsStyle(region.focus) <- "UCSC"
    
    gr <- subsetByOverlaps(gr, region.focus)
    
    if(length(gr) > 0){
        plot.df <- data.frame(start = start(gr), end = end(gr))
        
        p <- ggplot(plot.df) + geom_rect(aes(xmin = start, xmax = end, ymin=0, ymax=1), color = col, fill=col) + 
            theme_classic() + 
            ylab("") +
            xlab("") + 
            LegendOff() + 
            xlim(c(start(region.focus), end(region.focus))) +
            theme_void() +
            ylim(c(0, 1))
        
        return(p)   
    }
    else{
        return(NULL)
    }
}

# requires a dataframe with a chromosome column (chr), a position of the enhancer (enhancer_mid), and promoter (promoter_mid)
HiC.track <- function(HiC.Midpoints, curr.locus.gr, links.to.highlight=NULL, curv=0.5){
  
  HiC.Midpoints$highlight <- FALSE
  if(!is.null(links.to.highlight)){
    HiC.Midpoints[links.to.highlight,"highlight"] <- TRUE
  }
  
  curr.region <- HiC.Midpoints[HiC.Midpoints$chr == as.character(seqnames(curr.locus.gr)), ]
  
  hic.plot <- ggplot(mapping = aes(x = enhancer_mid, xend = promoter_mid, y = 0, yend = 0, color = highlight, size = highlight, alpha = highlight)) +
    
    geom_curve(data = subset(curr.region, enhancer_mid > promoter_mid),curvature = -curv, ncp = 1000) +
    geom_curve(data = subset(curr.region, enhancer_mid < promoter_mid),curvature = curv, ncp = 1000) + 
    
    xlim(c(start(curr.locus.gr), end(curr.locus.gr))) + 
    theme_void() + 
    LegendOff() + scale_color_manual(values = c("black","red")) + 
    scale_size_manual(values = c(0.3, 1)) + 
    scale_alpha_manual(values = c(0.2, 0.5)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.01, 0))
  
  hic.plot
}


snp.investigator <- function(snp.gr){
    
    bw.files <- c("ArchR/ArchR_heart_latest_noAtrium/GroupBigWigs/CellTypes/Cardiomyocyte-TileSize-100-normMethod-ReadsInTSS-ArchR.bw",
                  "ArchR/ArchR_heart_latest_noAtrium/GroupBigWigs/CellTypes/Fibroblast-TileSize-100-normMethod-ReadsInTSS-ArchR.bw")
    bw.list <- lapply(bw.files, function(x){rtracklayer::import(x)})
    
    bed.files <- c('GWAS/annotations_hg38/CM_specific_peaks.bed', 'GWAS/annotations_hg38/CM_shared_peaks.bed', 'GWAS/annotations_hg38/non_CM_peaks.bed')
    bed.list <- lapply(bed.files, function(x){rtracklayer::import(x)})
    
    
    #my.snp <- GRanges(seqnames = "chr1", ranges = IRanges(start = 170618199, end = 170618199))
    curr.locus.gr <- snp.gr + 10000
    
    plot.list <- list()
    #plot.list[["gene.track"]] <- geneTrackPlot(curr.locus.gr = curr.locus.gr, collapseTranscripts = "longest", genome = "hg38")
    
    plot.list[["cm.track"]] <- plotATAC(gr = bw.list[[1]], region.focus = curr.locus.gr) + geom_vline(xintercept = 170618199)
    plot.list[["fibro.track"]] <- plotATAC(gr = bw.list[[2]], region.focus = curr.locus.gr) + geom_vline(xintercept = 170618199)
    
    plot.list[["CM_specific_peaks"]] <- plotRect(gr = bed.list[[1]], region.focus = curr.locus.gr, col = 'blue')
    plot.list[["CM_shared_peaks"]] <- plotRect(gr = bed.list[[2]], region.focus = curr.locus.gr, col='red')
    plot.list[["non_CM_peaks"]] <- plotRect(gr = bed.list[[3]], region.focus = curr.locus.gr, col='green')
    
    cowplot::plot_grid(plotlist = plot.list, align = 'v', ncol = 1, axis = 'b')
}

