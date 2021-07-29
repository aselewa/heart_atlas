
setwd('/project2/gca/aselewa/heart_atlas_project/')

library(Gviz)
library(lattice)
library(gridExtra)

palette <- readRDS('notebooks/palette.rds')
markers <- c("TNNT2","VWF","DCN","CD8A","CD14","PLP1","RGS5","TAGLN")
genomic.annots <- readRDS('hg38_genomic_annotations.gr.rds')
gene.annots <- genomic.annots$genes

markers.gr <- gene.annots[match(markers, gene.annots$gene_name),]
n <- length(markers)

#satac <- loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/')
#getGroupBW(satac, groupBy = "CellTypes")
bigwig.files <- list.files('ArchR/ArchR_heart_latest_noAtrium/GroupBigWigs/CellTypes/', pattern = '*.bw', full.names = T)
bigwig.names <- sub("-.*","",list.files('ArchR/ArchR_heart_latest_noAtrium/GroupBigWigs/CellTypes/', pattern = '*.bw', full.names = F))
bigwig.names <- sub('\\.', " ", bigwig.names)

bigwig.list <- lapply(bigwig.files, FUN = function(x){rtracklayer::import(x)})
names(bigwig.list) <- bigwig.names

my.gene.track <- GeneRegionTrack(range = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, genome = "hg38", fill = "black", col = "black")
symbol.from.entrez <- mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                             keys=sub("\\.\\d+$", "", gene(my.gene.track)), 
                             keytype="ENTREZID", column="SYMBOL")
symbol.from.entrez[sapply(symbol.from.entrez, is.null)] <- NA
symbol(my.gene.track) <- unlist(symbol.from.entrez)
symbol(my.gene.track) <- ifelse(is.na(symbol(my.gene.track)), gene(my.gene.track), symbol(my.gene.track)) 

atac.plots <- list()
gene.plots <- list()
for(i in 1:n){
   curr.gene <- markers[i]
   curr.locus.gr <- markers.gr[markers.gr$gene==curr.gene,] + 5000
   
   atac.plots[[curr.gene]] <- lapply(names(bigwig.list), FUN = function(x) { DataTrack(range = bigwig.list[[x]], type='histogram', genome = "hg38", 
                                                                                chromosome = seqnames(curr.locus.gr),
                                                                                fill.histogram = palette[x], 
                                                                                col.histogram = palette[x]) })
   gene.plots[[curr.gene]] <- my.gene.track
}


markers.df <- data.frame(markers = factor(markers, levels = markers))
p1 <- xyplot(1 ~ markers | markers, data = markers.df, layout=c(8,1), between = list(x = 0.1),
             panel = function(x) { plotTracks(atac.plots[[x]], 
                                               chromosome = as.character(seqnames(markers.gr[markers.gr$gene_name==x,] + 5000)),
                                               from = start(markers.gr[markers.gr$gene_name==x,] + 5000), 
                                               to = end(markers.gr[markers.gr$gene_name==x,] + 5000), 
                                               panel.only = T, 
                                               ylim = c(0, 150), 
                                               add = T) },
             scales = list(draw = FALSE), xlab = NULL, ylab = NULL)


p2 <- xyplot(1 ~ markers | markers, data = markers.df, layout=c(8,1), strip=FALSE, 
             panel = function(x) { plotTracks(gene.plots[[x]], 
                                             chromosome = as.character(seqnames(markers.gr[markers.gr$gene_name==x,] + 5000)),
                                             from = start(markers.gr[markers.gr$gene_name==x,] + 5000), 
                                             to = end(markers.gr[markers.gr$gene_name==x,] + 5000), 
                                             panel.only = T, 
                                             collapseTranscripts = "longest",
                                             ylim = c(0, 150), 
                                             add = T) },
             between = list(x = 0.1),
             par.settings = list(axis.line = list(col = 0)),
             scales = list(draw = FALSE), xlab = NULL, ylab = NULL)

pdf('manuscript_figures/figure1/gviz_celltype_tracks.pdf', width = 10, height=6)
gridExtra::grid.arrange(p1, p2, ncol=1, heights=c(8.5,2))
dev.off()




