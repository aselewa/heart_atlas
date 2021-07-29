library(ArchR)
library(tidyverse)
require(plyranges)

setwd('/project2/gca/aselewa/heart_atlas_project/')
palette <- readRDS('notebooks/palette.rds')
source('R/analysis_utils.R')


# load some basics before figures
satac <- suppressMessages(loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/'))

peaks <- getPeakSet(satac)
peaks$peakID <- GRToString(peaks)
peaks <- peaks[peaks$peakType!="Promoter",]
peak.markers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')

genomic.annots <- readRDS('genomic_annotations/hg38_gtf_genomic_annots.gr.rds')
gene.annots <- genomic.annots$genes

### annotate all peaks with DA test results
celltype_ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")
peak.markers <- lapply(celltype_ideal_order, function(x){peak.markers[[x]]})
names(peak.markers) <- celltype_ideal_order
peak.markers.str <- unlist(lapply(peak.markers, function(x){GRToString(x)}), use.names = F)
peaks.celltypes <- data.frame(peakID=peak.markers.str, cellType = factor(rep(names(peak.markers), lengths(peak.markers)), levels = names(peak.markers)))

### load co-accessibility results
enhancer.coacc <- readRDS('ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coacc_ENHANCERS_AllCellTypes_overlapCutoff50_k200_corr_cut_-1_maxDist_1Mb_hg38.gr.rds')
enhancer.coacc <- enhancer.coacc[enhancer.coacc$correlation > 0.5,]
 
enhancer.coacc.tbl <- enhancer.coacc %>% 
    as_tibble() %>%  
    dplyr::select(peakID, coacc_gene_name, correlation, distToCoaccGene) %>%
    group_by(peakID) %>% arrange(-correlation) %>% slice(1) %>% ungroup()

peak.info.coacc.df <- left_join(peak.info.df, enhancer.coacc.tbl, on = "peakID") %>% 
    filter(!is.na(coacc_gene_name)) %>% 
    distinct(peakID, coacc_gene_name, .keep_all = T)

npeaks.per.gene <- peak.info.coacc.df %>% 
    group_by(coacc_gene_name) %>%
    summarise(peakCount = dplyr::n()) %>% 
    mutate(peakFreq = peakCount/sum(peakCount))

bks <- c(0, 2, 10, 20, 30, 40, 50, 60, 10000)
labs <- c("1-2","3-10","11-20", "21-30","31-40", "41-50","51-60","61+")
peakCount.bins <- npeaks.per.gene %>% summarise(nGenes = as.vector(table(cut(peakCount, breaks = bks, labels = labs)))) %>% mutate(nPeaks = labs)
peakCount.bins$nPeaks <- factor(peakCount.bins$nPeaks, levels = labs)

# 3A is cartoon of co-accessibility

# 3B with inset, number of peaks per gene
pdf('manuscript_figures/figure3/Fig3_percent_peaks_distribution.pdf', width = 8, height=6)
ggplot(peakCount.bins, aes(x = nPeaks, y = nGenes/1e3)) + 
    geom_bar(stat='identity', position='dodge', fill="#238b45", width=0.7) + ggClean() +
    ylab('Number of Genes (x 1000)') + xlab('Number of Peaks') + scale_y_continuous(breaks = seq(0, 8, 1))
dev.off()

dist.df <- data.frame(dist = log10(peak.info.coacc.df$distToCoaccGene),
                      type = peak.info.coacc.df$correlation>0)

pdf('manuscript_figures/figure3/Fig3A_2_distance_dist.pdf', width = 8, height=6)
ggplot(dist.df, aes(x = dist)) + geom_density(adjust=2) + ggClean() + xlim(c(0,7))
dev.off()


# gene expression of peaks linked

srna.exp <- Seurat::AverageExpression(srna)
srna.exp.mat <- as.matrix(log2(srna.exp$RNA + 1))
rSDs <- matrixStats::rowSds(srna.exp.mat)
srna.exp.mat <- srna.exp.mat[rSDs > 0,]
srna.exp.mat.z <- sweep(x = srna.exp.mat - rowMeans(srna.exp.mat), MARGIN = 1, STATS = matrixStats::rowSds(srna.exp.mat), FUN = '/')

nCellType <- npeaks.per.gene.celltype %>% group_by(coacc_gene_name) %>% summarise(nType = dplyr::n()) 
unique.npeaks.per.gene <- npeaks.per.gene.celltype[npeaks.per.gene.celltype$coacc_gene_name %in% nCellType$coacc_gene_name,]
genes.expressed.df <- unique.npeaks.per.gene[unique.npeaks.per.gene$coacc_gene_name %in% rownames(srna.exp.mat.z),] 

cm.genes <- genes.expressed.df$coacc_gene_name[genes.expressed.df$cellType=="Cardiomyocyte"]
npeaks <- genes.expressed.df$peakCount[genes.expressed.df$cellType=="Cardiomyocyte"]

srna.egene.mat.z.sub <- srna.exp.mat.z[cm.genes, celltype_ideal_order]
expOrder <- order(srna.egene.mat.z.sub[,"Cardiomyocyte"], decreasing = T)
srna.egene.mat.z.sub <- srna.egene.mat.z.sub[expOrder,]
npeaks <- npeaks[expOrder]

pdf('manuscript_figures/figure3/Fig3F_Genes_Linked_Expression.pdf', width = 10, height=11)
Heatmap(srna.egene.mat.z.sub, 
        show_row_names = F, 
        cluster_rows = F, 
        cluster_columns = F, 
        left_annotation = rowAnnotation(npeaks = npeaks),
        col = circlize::colorRamp2(c(-2,0,2), c("lightblue","white","firebrick")), name = 'log2 Expression',
        use_raster = T)
dev.off()


# # Fig 3 supplement: correlation of # of peaks with transcriptional units
# 
gene.annots.tbl <- gene.annots %>% as_tibble() %>% select(seqnames, start, end, gene_name) %>% mutate(gene_length = end-start) %>% dplyr::rename(coacc_gene_name = gene_name)

bks <- c(0, 2, 50, 10000)
labs <- c("<=2", "3-50","50+")
npeaks.per.gene.length <- left_join(npeaks.per.gene, gene.annots.tbl, on="coacc_gene_name") %>%
    mutate(nGenes = cut(peakCount, breaks = bks, labels = labs))

pdf('manuscript_figures/figure3/Fig3_Supplementary_GeneLength_vs_nPeaks.pdf', width = 8, height=6)
ggplot(npeaks.per.gene.length, aes(x =  nGenes, y = log10(gene_length), color=nGenes)) + geom_violin(color = 'black') + geom_boxplot(width=0.1)  + ggClean() +
    xlab('Number of Peaks') + ylab('log10 Gene Length')
dev.off()

# 3E is GO enrichment, run script R/GO_enrichment.R

# 3F cell-type specificity of each gene linked
npeaks.linked <- peak.info.coacc.df %>% filter(cellType != "Shared") %>% group_by(coacc_gene_name) %>% summarise(nPeaks = dplyr::n())

ncelltypes.per.gene <- peak.info.coacc.df %>% 
    filter(cellType != "Shared") %>% 
    group_by(cellType, coacc_gene_name) %>% 
    summarise(cellTypeFreq = dplyr::n())

gene.celltype.prob <- ncelltypes.per.gene %>% 
    group_by(coacc_gene_name) %>% 
    mutate(geneProb = cellTypeFreq / sum(cellTypeFreq)) %>% 
    select(-cellTypeFreq) 

gene.celltype.prob.wide <- gene.celltype.prob %>% pivot_wider(values_from = geneProb, names_from=cellType)

gene.celltype.prob.wide <- as.data.frame(gene.celltype.prob.wide)
rownames(gene.celltype.prob.wide) <- gene.celltype.prob.wide$coacc_gene_name
gene.celltype.prob.wide$coacc_gene_name <- NULL
gene.celltype.prob.wide[is.na(gene.celltype.prob.wide)] <- 0
gene.celltype.prob.wide <- gene.celltype.prob.wide[,celltype_ideal_order]

gene.counts <- rowSums(srna@assays$RNA@counts)
srna.filt <- srna[gene.counts > 5000,]
srna.exp <- Seurat::AverageExpression(srna.filt)
mean.gene.exp <- log2(rowMeans(srna.exp$RNA))
srna.exp.prob <- srna.exp$RNA / rowSums(srna.exp$RNA)
srna.exp.prob <- srna.exp.prob[, celltype_ideal_order]
same_genes <- intersect(rownames(srna.exp.prob), rownames(gene.celltype.prob.wide))

gene.celltype.prob.wide <- gene.celltype.prob.wide[same_genes,]
srna.exp.prob <- srna.exp.prob[same_genes,]
mean.gene.exp <- mean.gene.exp[same_genes]
    
cor.res <- rep(0, length(same_genes))
for(i in 1:length(same_genes)){
    cor.res[i] <- cor(x = as.numeric(gene.celltype.prob.wide[i,]), y = as.numeric(srna.exp.prob[i,]))
}

cor.df <- data.frame(coacc_gene_name = same_genes, pearson = cor.res, meanExp = mean.gene.exp) %>% left_join(., npeaks.linked, on="coacc_gene_name")

bks <- c(0, 2, 10, 20, 30, 40, 50, 60, 10000)
labs <- c("1-2","3-10","11-20", "21-30","31-40", "41-50","51-60","61+")
cor.df <- cor.df %>% mutate(nPeaks = cut(cor.df$nPeaks, breaks = bks, labels = labs))
cor.df$nPeaks <- factor(cor.df$nPeaks, levels = labs)

data_summary <- function(x) {
    m <- mean(x)
    ymin <- max(min(x), m-sd(x))
    ymax <- min(max(x), m+sd(x))
    return(c(y=m,ymin=ymin,ymax=ymax))
}
markers <- c("TNNT2","MYBPC3","MYH7","NPPA","RGS5","ABCC9","MYH11","TAGLN","DCN","PDGFRA","PECAM1","VWF","PLP1","CD8A","LCK","CD14","FOLR2")
cor.df$highlight <- cor.df$coacc_gene_name %in% markers

pdf('manuscript_figures/figure3/Fig3F_CellTypeSpecificty_Pearson.pdf', width = 10, height=6)
p <- ggplot(cor.df, aes(x=nPeaks, y=pearson)) + 
    ggrastr::geom_jitter_rast(position=position_jitter(0.3), aes(color=highlight, alpha=highlight, size=highlight)) + ggClean() + LegendOff() + scale_color_manual(values = c("black","green")) +
    scale_alpha_manual(values = c(0.3, 1)) + coord_cartesian(ylim = c(-1, 1))
p + stat_summary(fun.data=data_summary, color="red")
dev.off()s

ggplot(cor.df, aes(x=nPeaks, y=meanExp)) + geom_boxplot() + ggClean() + ylab('log2 Expression') + xlab('Number of Peaks')


# 3D, overlap with RNA genes
rna.markers <- readRDS('seurat/diff_expr_markers.df.rds')
rna.markers$cluster <- as.character(rna.markers$cluster)

celltypes <- celltype_ideal_order
gene.overlap <- matrix(0, nrow=length(celltypes), ncol=length(celltypes))
for(i in 1:length(celltypes)){
    for(j in 1:length(celltypes)){
        gene.overlap[i,j] <- sum(rna.markers$gene[rna.markers$cluster==celltypes[i]] %in% peak.info.coacc.df$coacc_gene_name[peak.info.coacc.df$cellType==celltypes[j]])/sum(rna.markers$cluster==celltypes[i])
    }
}

rownames(gene.overlap) <- celltypes
colnames(gene.overlap) <- paste0("PeakGenes_", celltypes)

pdf('manuscript_figures/figure3/Fig3B_overlap_w_rna.pdf', width = 8, height=6)=
ComplexHeatmap::Heatmap(gene.overlap, 
                        cluster_rows = F, 
                        cluster_columns = F, 
                        col = circlize::colorRamp2(c(0, 0.3, 0.8), c("lightblue","white","firebrick")))
dev.off()


# 3C gene track plot

#interest.eqtl <- eqtl.snp.gr[eqtl.snp.gr$eGene=="ESRRB",]
#curr.locus.gr <- GRanges(seqnames = seqnames(interest.eqtl), ranges = IRanges(start = start(interest.eqtl)-14000, end = end(interest.eqtl)+14000))

offset <- 100000
my.gene = "MYBPC3"
gene.of.interest <- gene.annots[gene.annots$gene_name==my.gene]
gene.start <- ifelse(strand(gene.of.interest) == "+", start(gene.of.interest), end(gene.of.interest))
curr.locus.gr <- GRanges(seqnames = seqnames(gene.of.interest), ranges = IRanges(start = gene.start-offset, end = gene.start+offset))

bw.files <- list.files(path = 'ArchR/ArchR_heart_latest_noAtrium/GroupBigWigs/CellTypes/', pattern = '*.bw$', full.names = T)
bw.files <- bw.files[c(1,2,3)]
peaks <- satac@peakSet

bw.list <- lapply(bw.files, function(x){rtracklayer::import(x)})
names(bw.list) <- sub('\\.', ' ', sub('-.*', '', basename(bw.files)))

atac.tracks <- lapply(names(bw.list), function(x){
    DataTrack(range = bw.list[[x]], type = 'h', col = palette[x], name = x, showAxis=FALSE, ylim=c(0, 0.2))
 })

union.calls <- DataTrack(range = peaks, type = 'h', genome = "hg38", name = "Union", col = "grey", showAxis=FALSE, ylim = c(0,1))

bamNormFunction <- function(x,libsize){
    x/libsize
}
cm.rna.bam <- DataTrack(range = '/project2/gca/Heart_Atlas/Nuc_seq/SP-HE-HE200915RNA-397RNA/output/Cardiomyocyte.bam',
                     type = 'h', genome = 'hg38', name = 'RNA', col = palette["Cardiomyocyte"], 
                     transformation = function(x){bamNormFunction(x,11188127)}, ylim=c(1e-5, 2e-4),
                     col.baseline='white')

endo.rna.bam <- DataTrack(range = '/project2/gca/Heart_Atlas/Nuc_seq/SP-HE-HE200915RNA-397RNA/output/Endothelial.bam',
                        type = 'h', genome = 'hg38', name = 'RNA', col = palette["Endothelial"], 
                        transformation = function(x){bamNormFunction(x,20603239)}, ylim=c(1e-5, 2e-4))

fibro.rna.bam <- DataTrack(range = '/project2/gca/Heart_Atlas/Nuc_seq/SP-HE-HE200915RNA-397RNA/output/Fibroblast.bam',
                          type = 'h', genome = 'hg38', name = 'RNA', col = palette["Fibroblast"], 
                          transformation = function(x){bamNormFunction(x,16849499)}, ylim=c(1e-5, 2e-4))

rna.tracks <- list(cm.rna.bam, endo.rna.bam, fibro.rna.bam)

promoter.gr <- StringToGR(enhancer.coacc$promoterID)
enhancer.promoter <- GenomicInteractions::GenomicInteractions(anchor1 = enhancer.coacc, anchor2 = promoter.gr)
enhancer.promoter$counts <- round(100*enhancer.coacc$correlation)
enhancer.promoter <- enhancer.promoter[enhancer.promoter$anchor1.coacc_gene_name == my.gene,]

coacc.track <- InteractionTrack(enhancer.promoter, name = "Coaccess")
dpars <- list(col.interactions="red", 
              col.anchors.fill ="blue", 
              col.anchors.line = "black", 
              interaction.dimension = "height",  
              interaction.measure = "counts", 
              plot.trans = FALSE, 
              plot.outside = FALSE,  
              col.outside="lightblue",  
              anchor.height = 0.1)
displayPars(coacc.track) <- dpars

gene.track <- knownGeneObject(curr.locus.gr = curr.locus.gr, genome = "hg38")

pdf('manuscript_figures/figure3/MYBPC3_coaccessibility_track.pdf', width=8, height=4)
plotTracks(c(coacc.track, union.calls, atac.tracks, rna.tracks, gene.track),
           chromosome = as.character(seqnames(curr.locus.gr)),
           transcriptAnnotation = "symbol",
           collapseTranscripts= 'longest',
           from = start(curr.locus.gr),
           to = end(curr.locus.gr), 
           sizes = c(1, 0.2, rep(0.5, 3), rep(0.5, 3), 1),
           panel.only = F)
#plotTracks(atrack, chromosome = as.character(seqnames(curr.locus.gr)), from = start(curr.locus.gr), to = end(curr.locus.gr), ylim=c(0,10))
dev.off()


