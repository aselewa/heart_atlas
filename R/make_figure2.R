library(ArchR)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

# load and define basics
satac <- loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/')
peak.set <- satac@peakSet
satac$aggregate <- "agg"
celltype_ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")

### 2A - visualize accessibility
mean.peak <- getGroupSE(ArchRProj = satac, 
                        useMatrix = 'PeakMatrix', 
                        groupBy = 'CellTypes', 
                        divideN = T, 
                        scaleTo = 10^5)

set.seed(100)

mean.peak.mat <- mean.peak@assays@data$PeakMatrix
mean.peak.mat.norm <- sweep(x = mean.peak.mat - rowMeans(mean.peak.mat), MARGIN = 1, STATS = matrixStats::rowSds(mean.peak.mat), FUN = '/')
rData <- rowData(mean.peak)
rownames(mean.peak.mat) <- paste0(as.character(rData$seqnames), ':', rData$start, '-', rData$end)

nClust = 14
row.clust.res <- kmeans(mean.peak.mat.norm, centers = nClust)
labs <- row.clust.res$cluster
rAnno <- rowAnnotation(cluster = factor(sort(labs)))
Heatmap(mean.peak.mat.norm[names(sort(labs)),], cluster_rows = F, cluster_columns = F, show_row_names = F, left_annotation = rAnno)

ideal_order <- c("Cardiomyocyte1","Cardiomyocyte2","Smooth Muscle", "Pericyte","Endothelial1","Endothelial2",
                 "Fibroblast","Neuronal","Lymphoid","Myeloid1","mix1","mix2","mix3","mix4")
lab_map <- data.frame(num = 1:nClust,
                      celltype = c("Smooth Muscle", "Myeloid1", "Lymphoid", "Cardiomyocyte1", "mix1", "Endothelial1", "Fibroblast", "Pericyte",
                                   "mix2","mix3","Endothelial2","Cardiomyocyte2","Neuronal","mix4"),
                      stringsAsFactors = F)

labs.named <- plyr::mapvalues(x = labs, from = lab_map$num, to = lab_map$celltype)
labs.named <- factor(labs.named, levels = ideal_order)
labs.sort <- sort(labs.named)
mean.peak.mat.orderd <- mean.peak.mat.norm[names(labs.sort), celltype_ideal_order]

pdf('manuscript_figures/Fig2A_peak_heatmap.pdf', width=14, height=16)
Heatmap(matrix = mean.peak.mat.orderd, cluster_rows = F, 
        cluster_columns = F, show_row_names = F, 
        name = "Row-normalized accessibility",
        col = circlize::colorRamp2(c(-2, 0, 2), c("lightblue","white","firebrick")), 
        column_split = factor(colnames(mean.peak.mat), levels=colnames(mean.peak.mat)),
        row_title = NULL)
dev.off()

# 2B
# genomic distribution of peaks
markers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')
markers <- lapply(celltype_ideal_order, function(x){markers[[x]]})
names(markers) <- celltype_ideal_order
DA.peaks <- unlist(lapply(markers, function(x){GRToString(x)}), use.names = F)
peaks.celltypes <- data.frame(peakID=DA.peaks, cellType = factor(rep(names(markers), lengths(markers)), levels = names(markers)))
    
peak.set <- getPeakSet(satac)
peak.set.str <- GRToString(peak.set)
peak.set$peakID <- peak.set.str

# categories non DA peaks into number of cell-types shared
shared.peaks <- peak.set[!(peak.set$peakID %in% peaks.celltypes$peakID),]
mean.peak.mat.sub <- mean.peak.mat[shared.peaks$peakID,]

for(i in 1:ncol(mean.peak.mat.sub)){
  curr.thresh <- as.numeric(quantile(mean.peak.mat.sub[,i], 0.75))
  mean.peak.mat.sub[,i] <- 1*(mean.peak.mat.sub[,i] > curr.thresh)
}
peak.count.df <- data.frame(cellType = rowSums(mean.peak.mat.sub), peakID = rownames(mean.peak.mat.sub))

non.da.bol <- peak.count.df$cellType == 0 | peak.count.df$cellType == 1
two.three.bol <- peak.count.df$cellType >=2 & peak.count.df$cellType <=3
four.plus.bol <- peak.count.df$cellType >=4

peak.count.df$cellType[non.da.bol] <- "Non-DA"
peak.count.df$cellType[two.three.bol] <- "2-3"
peak.count.df$cellType[four.plus.bol] <- "4+"
rownames(peak.count.df) <- NULL

peaks.celltypes.shared <- bind_rows(peaks.celltypes, peak.count.df)

peak.info.df <- peak.set %>% as_tibble() %>% left_join(., peaks.celltypes.shared, on="peakID") 
peak.info.df$cellType <- factor(peak.info.df$cellType, levels = c(celltype_ideal_order, "Non-DA","2-3","4+"))

peak.dist.df <- peak.info.df %>% group_by(cellType, peakType) %>% 
    summarise(numberPeaks = n()) %>% 
    group_by(cellType) %>% 
    mutate(PropPeaks = numberPeaks/sum(numberPeaks)) %>%
    mutate(TotPeaks = sum(numberPeaks))

peak.dist.df$cellType <- factor(peak.dist.df$cellType, levels = c(celltype_ideal_order, "Non-DA","2-3","4+"))

p2 <- ggplot(peak.dist.df, aes(x = cellType, y = PropPeaks, fill = peakType)) + 
    geom_bar(stat='identity') + 
    ggClean(rotate_axis = T) + 
    scale_fill_manual(values = c("#aaaaaa","#bee6af","#74c476","#238b45"))

pdf('manuscript_figures/figure2/Fig2A_2_peaks_per_celltype.pdf', width=6, height=4)
p2
dev.off()

##### Fig 2C 
### distance to nearest gene

peak.set <- getPeakSet(satac)
peak.set$peakID <- GRToString(peak.set)

peaks.info.df <- peak.set %>% 
  as_tibble() %>% 
  left_join(., peaks.celltypes)

palette <- readRDS('notebooks/palette.rds')
pdf('manuscript_figures/figure2/Fig2A_3_distToNearestGene.pdf', width=8, height=5)

ggplot(peak.info.df, aes(x = log10(distToGeneStart + 1),  color = cellType)) + 
  geom_density(adjust = 1, size = 1.2) + 
  ggClean() + 
  xlab('Distance (bp)') + 
  scale_color_manual(values = c(palette[celltype_ideal_order],c("lightgrey","grey","black"))) 

dev.off()

# Fig Supplemental 2 - donor correlation

mean.donor.peak <- getGroupSE(ArchRProj = satac, 
                        useMatrix = 'PeakMatrix', 
                        groupBy = 'individual', 
                        divideN = T, 
                        scaleTo = 10^4)

mean.donor.peak.mat <- mean.donor.peak@assays@data$PeakMatrix
donor.corr.mat <- cor(mean.donor.peak.mat)

pdf('manuscript_figures/Fig2B_1_donor_corr.pdf', width=12, height=10)
Heatmap(matrix = donor.corr.mat, cluster_rows = F, 
        cluster_columns = F, show_row_names = T, row_names_side = "left",
        name = "Pearson",
        col = circlize::colorRamp2(c(0.8, 1), c("white","firebrick")), 
        row_title = NULL,
        column_title = NULL,
        row_gap = unit(1, "mm"),
        column_gap = unit(1, "mm"),
        na_col = "white",
        use_raster = T)
dev.off()


# Fig Supplemental 2 - region correlation

mean.peak.region <- getGroupSE(ArchRProj = satac, 
                               useMatrix = 'PeakMatrix', 
                               groupBy = 'regions', 
                               divideN = T,
                               scaleTo = 10^4)

mean.peak.region.mat <- mean.peak.region@assays@data$PeakMatrix
region.cor.mat <- cor(mean.peak.region.mat)

pdf('manuscript_figures/Fig2B_2_region_corr.pdf', width=12, height=10)
Heatmap(matrix = region.cor.mat, cluster_rows = F, 
        cluster_columns = F, show_row_names = T, row_names_side = "left",
        name = "Pearson",
        col = circlize::colorRamp2(c(0.8, 1), c("white","firebrick")), 
        row_title = NULL,
        column_title = NULL,
        row_gap = unit(1, "mm"),
        column_gap = unit(1, "mm"),
        na_col = "white",
        use_raster = T)
dev.off()

# Figure 2D correlation with ENCODE DNase
markers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')

#fnames <- c('ENCODE/H3K27ac/heart_left_ventricle.bed.gz', 'ENCODE/H3K27ac/heart_right_ventricle.bed.gz')
fnames <- c(list.files('ENCODE/DNase', pattern = "*.bed.gz$", full.names = T), 'ENCODE/H3K27ac/Heart_LVRV.bed.gz')
            #'ENCODE/H3K27ac/heart_left_ventricle.bed.gz', 'ENCODE/H3K27ac/heart_right_ventricle.bed.gz')
encode.file.df <- data.frame(path = fnames, type = sub('_hg38','',sub('.bed.gz','',basename(fnames))), assay = basename(dirname(fnames)), stringsAsFactors = F)
                    
encode.beds <- lapply(encode.file.df$path, function(x){ 
    encode <- vroom::vroom(x, col_names = F); 
    encode.gr <- GRanges(seqnames = encode$X1, ranges = IRanges(start = encode$X2, end = encode$X3));
    encode.gr
    }
)
names(encode.beds) <- paste(encode.file.df$type)

encode.GR <- unlist(GRangesList(encode.beds))
overlap.res <- subsetByOverlapProp(q = peak.set, s = encode.GR, minPoverlap = 0.1)
union.peaks.In.str <- GRToString(peak.set[unique(queryHits(overlap.res)),])

getPeakEnrichment <- function(gr, encode.gr){
    rangesIn <- subsetByOverlapProp(q = gr, s = encode.gr, minPoverlap = 0.1)
    length(unique(queryHits(rangesIn))) / length(gr)
}

markers$`Smooth Muscle` <- NULL
markers$Neuronal <- NULL
celltypes <- names(markers)
enrich.result <- list()
for( c in celltypes){
    enrich.result[[c]] <- sapply(encode.beds, function(x){ getPeakEnrichment(markers[[c]], x) })
}
enrich.df <- as.data.frame(enrich.result)
dnase.mat <- enrich.df[encode.file.df$assay=="DNase",]
h3k.mat <- enrich.df[encode.file.df$assay=="H3K27ac",]

dnase.mat <- dnase.mat[c("Heart_RVLV","CardiacMuscle","DermisEndothelialBlood","Tcell","Bcell","CD14Monocyte"), 
                       c("Cardiomyocyte","Pericyte","Endothelial","Fibroblast","Lymphoid","Myeloid")]

h3k.mat <- h3k.mat[,c("Cardiomyocyte","Pericyte","Endothelial","Fibroblast","Lymphoid","Myeloid")]

pdf('manuscript_figures/figure2/Overlap_encode_dnase_enrichment.pdf', width=12, height=10)

Heatmap(matrix = dnase.mat, cluster_rows = F, 
        cluster_columns = F, show_row_names = T, row_names_side = "left",
        name = "Proportion Overlap",
        col = circlize::colorRamp2(c(0,0.25,0.5), c("lightblue","white","firebrick")), 
        row_title = "ENCODE DNase",
        column_title = "Cell Types",
        row_gap = unit(1, "mm"),
        column_gap = unit(1, "mm"),
        na_col = "white",
        rect_gp = gpar(col = "black", lwd = 0.5),
        use_raster = F)

dev.off()


h3k.mat <- h3k.mat %>% 
    pivot_longer(cols=  everything()) %>% 
    rename(propOverlap = value, celltypes=name)

pdf('manuscript_figures/figure2/Overlap_encode_h3k27ac_enrichment.pdf', width=4, height=4)
ggplot(h3k.mat, aes(x=celltypes, y=propOverlap)) + 
    geom_bar(stat='identity', width=0.7) + 
    ggClean(rotate_axis = T) + 
    ylab('Proportion Overlap') + 
    xlab('')
dev.off()


# 2E
# motif enrichment plot
archrPeaks <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_markerPeaks.rds')
archr_motifUp <- suppressMessages(peakAnnoEnrichment(seMarker = archrPeaks, 
                                                     ArchRProj = satac, 
                                                     peakAnnotation = "Motif", 
                                                     cutOff = "FDR <= 0.01 & Log2FC >= 1"))
TF.enrich.mat <- archr_motifUp@assays@data$Enrichment
TF.fdr.mat <- archr_motifUp@assays@data$mlog10Padj
pvalue.mat <- archr_motifUp@assays@data$mlog10p

tf.info <- data.frame(tf.names = rownames(TF.fdr.mat), stringsAsFactors = F) %>% 
    mutate(gene_name = sub('_.*', '', tf.names)) %>% 
    mutate(max_enrichment = apply(X = TF.fdr.mat, MARGIN = 1, FUN = function(x){max(x)})) %>%
    group_by(gene_name) %>% filter(max_enrichment == max(max_enrichment))

TF.fdr.mat <- TF.fdr.mat[tf.info$tf.names,]
pvalue.mat <- pvalue.mat[tf.info$tf.names,]

rownames(TF.enrich.mat) <- tf.info$gene_name
rownames(TF.fdr.mat) <- tf.info$gene_name
rownames(pvalue.mat) <- tf.info$gene_name

ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")
TF.fdr.mat$TFs <- rownames(TF.fdr.mat)
TF.fdr.mat.long <- TF.fdr.mat %>% pivot_longer(cols = all_of(ideal_order))

TF.enrich.mat$TFs <- rownames(TF.enrich.mat)
TF.enrich.mat.long <- TF.enrich.mat %>% pivot_longer(cols = all_of(ideal_order))
TF.fdr.mat.long$Enrichment <- TF.enrich.mat.long$value
TF.fdr.mat.long <- TF.fdr.mat.long %>% dplyr::rename(Celltypes = name, mlog10_fdr=value)

TF.fdr.mat.long %>% arrange(Celltypes, -mlog10_fdr, -Enrichment) %>% write_tsv('misc/All_TFs_celltype_enrichment.tsv')

geneScore.cor.res <- read_tsv('TF_Correlation_Results/All_TFMotif_GeneScore_Correlations_CisBP.tsv')

ct.specific.tf <- list()
for(t in ideal_order){
    curr <- TF.fdr.mat.long[TF.fdr.mat.long$Celltypes == t,]
    ct.specific.tf[[t]] <- curr$TFs[order(curr$mlog10_fdr, decreasing = T)]
}
ct.motifs <- unlist(ct.specific.tf, use.names = F)
ct.motifs <- ct.motifs[!duplicated(ct.motifs)]
enrich.df <- pvalue.mat[ct.motifs, ideal_order]

geneScore.cor.res <- geneScore.cor.res[abs(geneScore.cor.res$GeneScore_Correlation) > 0.5,]
enrich.df <- enrich.df[rownames(enrich.df) %in% geneScore.cor.res$gene_name, ]

# save final summary of TFs
TF.fdr.mat.long %>% 
  filter(TFs %in% rownames(enrich.df)) %>% 
  dplyr::rename(gene_name = TFs) %>% 
  left_join(., geneScore.cor.res) %>%
  write_tsv('TF_Correlation_Results/Final_TFs_summary.tsv')


#tfs.to.highlight <- unlist(lapply(ct.specific.tf, function(x){sub('_.*','',x[1:3])}), use.names = F)
tfs.to.highlight <- c("TBX5","GATA4","MEF2A","MEF2B","MEF2C","TEAD1","HAND2","PRRX1","NKX25","ESRRB","ESRRA","NR2F2","PITX2")

pdf('manuscript_figures/figure2/TF_Motif_enriched.pdf', width = 7, height=6)

rAnno <- rowAnnotation(n = anno_mark(at = match(tfs.to.highlight, rownames(enrich.df)), labels = tfs.to.highlight, labels_gp = gpar(fontsize = 22)))
Heatmap(matrix = enrich.df, cluster_rows = F, cluster_columns = F, show_column_names = T, show_row_names = F,
        col = circlize::colorRamp2(c(0,80), c("white","firebrick")), right_annotation = rAnno,
        row_title = NULL, border = T, column_split = 1:ncol(enrich.df), column_gap = unit(0, "mm"),
        column_title = NULL,
        name = "Motif Enrichment (-log10 pval)", 
        use_raster = T)

dev.off()


## 2F -- archr overlays

genescore.tbx5 <- plotEmbedding(satac, colorBy = 'GeneScoreMatrix', name='ID3')
p1 <- genescore.tbx5 + ggClean() + ggtitle('Gene Scores: ID3') +xlab('') + ylab('') + theme(legend.position = 'right')
motifaccess.tbx5 <- plotEmbedding(satac, colorBy = 'MotifMatrix', name='z:ID3_38')
p2 <- motifaccess.tbx5 + ggClean() + ggtitle('Motif Access.: ID3') + xlab('') + ylab('') + theme(legend.position = 'right')
p1+p2

genescore.tbx5 <- plotEmbedding(satac, colorBy = 'GeneScoreMatrix', name='ETV1')
p1 <- genescore.tbx5 + ggClean() + ggtitle('Gene Scores: ETV1') +xlab('') + ylab('') + theme(legend.position = 'right')
motifaccess.tbx5 <- plotEmbedding(satac, colorBy = 'MotifMatrix', name='z:ETV1_320')
p2 <- motifaccess.tbx5 + ggClean() + ggtitle('Motif Access.: ETV1') + xlab('') + ylab('') + theme(legend.position = 'right')

p1+p2

