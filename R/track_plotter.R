library(ggplot2)
library(lattice)
source('/project2/gca/aselewa/heart_atlas_project/R/analysis_utils.R')

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
                             end = end(curr.locus.gr), name = "Genes")
  
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
                                                                 panel.only = F, 
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
plotATAC <- function(gr, region.focus, max.ylim=150, fill.col = "darkgreen"){
  gr <- subsetByOverlaps(gr, region.focus)
  score <- gr$score
  plot.df <- data.frame(start = start(gr), end = end(gr), score = score)
  
  p <- ggplot(plot.df) + geom_rect(aes(xmin = start, xmax = end, ymin=0, ymax=score), color = fill.col, fill=fill.col) + 
    theme_classic() + 
    ylab("") +
    xlab("") + 
    LegendOff() + 
    ylim(c(0, max.ylim)) +
    coord_cartesian(xlim=c(start(region.focus), end(region.focus))) +
      theme(axis.line.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.major.x=element_blank())
  
  return(p)
}

plotRect <- function(gr, region.focus, col){
    
    seqlevelsStyle(gr) <- "UCSC"
    seqlevelsStyle(region.focus) <- "UCSC"
    
    gr <- subsetByOverlaps(gr, region.focus)
    
    if(length(gr) > 0){
        plot.df <- data.frame(start = start(gr), end = end(gr))
    }
    else{
        plot.df <- data.frame(start = start(region.focus), end = end(region.focus))
        col <- "white"
    }
    p <- ggplot(plot.df) + geom_rect(aes(xmin = start, xmax = end, ymin=0, ymax=1), color = col, fill=col) + 
        theme_classic() + 
        ylab("") +
        xlab("") + 
        LegendOff() + 
        ylim(c(0, 1)) + 
        theme_void() +
        coord_cartesian(xlim = c(start(region.focus), end(region.focus)))
        
    
    return(p)   

}

plotRectGene <- function(gr, region.focus, col){
    
    seqlevelsStyle(gr) <- "UCSC"
    seqlevelsStyle(region.focus) <- "UCSC"
    
    gr <- subsetByOverlaps(gr, region.focus)
    
    if(length(gr) > 0){
        plot.df <- data.frame(start = start(gr), end = end(gr))
    }
    else{
        plot.df <- data.frame(start = start(region.focus), end = end(region.focus))
        col <- "white"
    }
    p <- ggplot(plot.df) + geom_rect(aes(xmin = start, xmax = end, ymin=0, ymax=1), color = col, fill=col) + 
        theme_classic() + 
        ylab("") +
        xlab("") + 
        LegendOff() + 
        ylim(c(0, 1)) + 
        theme_void() + 
        coord_cartesian(xlim = c(start(region.focus), end(region.focus)))
    
    
    return(p)   
    
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
    scale_color_manual(values = c("black","red")) + 
    scale_size_manual(values = c(0.3, 1)) + 
    scale_alpha_manual(values = c(0.2, 0.5)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.01, 0)) +
    theme_classic() + 
    theme(axis.line.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.major.x=element_blank()) +
    LegendOff()
  
  hic.plot
}




