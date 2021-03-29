library(ArchR)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
source('R/analysis_utils.R')
setwd('/project2/gca/aselewa/heart_atlas_project/')

satac <- loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/')
satac$aggregate <- "agg"

mean.peak <- getGroupSE(ArchRProj = satac, 
                        useMatrix = 'PeakMatrix', 
                        groupBy = 'CellTypes', 
                        divideN = T, 
                        scaleTo = 10^4)

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

celltype_ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")

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

# 2A 1
pdf('manuscript_figures/Fig2A_peak_heatmap.pdf', width=14, height=16)
Heatmap(matrix = mean.peak.mat.orderd, cluster_rows = F, 
        cluster_columns = F, show_row_names = F, 
        name = "Row-normalized accessibility",
        col = circlize::colorRamp2(c(-2, 0, 2), c("lightblue","white","firebrick")), 
        column_split = factor(colnames(mean.peak.mat), levels=colnames(mean.peak.mat)),
        row_title = NULL,
        column_title = NULL,
        row_gap = unit(1, "mm"),
        column_gap = unit(1, "mm"),
        use_raster = T)
dev.off()

# 2A 2
# genomic distribution of peaks
markers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')
markers <- lapply(celltype_ideal_order, function(x){markers[[x]]})
names(markers) <- celltype_ideal_order
DA.peaks <- unlist(lapply(markers, function(x){GRToString(x)}), use.names = F)
peaks.celltypes <- data.frame(peakID=DA.peaks, cellType = factor(rep(names(markers), lengths(markers)), levels = names(markers)))
    
peak.set <- getPeakSet(satac)
peak.set.str <- GRToString(peak.set)
peak.set$peakID <- peak.set.str

peak.info.df <- peak.set %>% as_tibble() %>% left_join(., peaks.celltypes, on="peakID")
levels(peak.info.df$cellType) <- c(levels(peak.info.df$cellType), "Shared")

peak.info.df$cellType[is.na(peak.info.df$cellType)] <- "Shared"

peak.dist.df <- peak.info.df %>% group_by(cellType, peakType) %>% 
    summarise(numberPeaks = n()) %>% 
    group_by(cellType) %>% 
    mutate(PropPeaks = numberPeaks/sum(numberPeaks)) %>%
    mutate(TotPeaks = sum(numberPeaks))

peak.dist.df.tot <- peak.dist.df %>% distinct(cellType, TotPeaks)

p2 <- ggplot(peak.dist.df, aes(x = cellType, y = PropPeaks, fill = peakType)) + 
    geom_bar(stat='identity') + 
    ggClean(rotate_axis = T) + 
    scale_fill_manual(values = c("#aaaaaa","#bee6af","#74c476","#238b45"))

pdf('manuscript_figures/figure2/Fig2A_2_peaks_per_celltype.pdf', width=6, height=4)
p2
dev.off()

# Fig 2A 3

mean.acc.bulk <- getGroupSE(ArchRProj = satac, useMatrix = "PeakMatrix", groupBy = "aggregate", divideN = T, scaleTo = 10^4)
mean.acc.bulk.mat <- as.data.frame(mean.acc.bulk@assays@data$PeakMatrix)
rData <- rowData(mean.acc.bulk)
mean.acc.bulk.mat$peakID <- paste0(as.character(rData$seqnames), ':', rData$start, '-', rData$end)

mean.acc.bulk.mat$isCellTypeSpec <- mean.acc.bulk.mat$peakID %in% peaks.celltypes$peakID

mean.acc.bulk.mat$cellType <- "aggregate"
mean.acc.bulk.mat <- mean.acc.bulk.mat[order(mean.acc.bulk.mat$agg, decreasing = T),] %>% rename(mean.acc = agg)
mean.acc.bulk.mat$idx <- 1:nrow(mean.acc.bulk.mat)

xvec <- seq(0, 0.4, by = 0.001)
shared.cdf <- sapply(xvec, function(x){mean(mean.acc.bulk.mat$mean.acc[mean.acc.bulk.mat$isCellTypeSpec==F] < x)})
spec.cdf <- sapply(xvec, function(x){mean(mean.acc.bulk.mat$mean.acc[mean.acc.bulk.mat$isCellTypeSpec==T] < x)})

pdf('manuscript_figures/figure2/Fig2A_3_shared_specific_cdf.pdf', width=6, height=5)
plot(x = xvec, y = shared.cdf, col='black', type = 'lines',xlab = 'Scaled Mean Accessibility', ylab = 'CDF')
lines(x = xvec, y = spec.cdf, col='red')
legend(x = 0.2, y = 0.6, legend = c("Shared", "Cell Type Specific"), col = c("black","red"), lty = 1)
dev.off()


peak.set <- getPeakSet(satac)
peak.set$peakID <- GRToString(peak.set)

peaks.info.df <- peak.set %>% as_tibble() %>% left_join(., peaks.celltypes)
levels(peaks.info.df$cellType) <- c(levels(peaks.info.df$cellType), "Shared")
peaks.info.df$cellType[is.na(peaks.info.df$cellType)] <- "Shared"

palette <- readRDS('notebooks/palette.rds')
pdf('manuscript_figures/figure2/Fig2A_3_distToNearestGene.pdf', width=8, height=5)

ggplot(peaks.info.df, aes(x = log10(distToGeneStart + 1),  color = cellType)) + geom_density(adjust = 1, size = 1.2) + ggClean() + xlab('Distance (bp)') +
    scale_color_manual(values = c(palette[celltype_ideal_order],"black")) 

dev.off()

# 2B 

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
        col = circlize::colorRamp2(c(0, 0.65, 1), c("lightblue","white","firebrick")), 
        row_title = NULL,
        column_title = NULL,
        row_gap = unit(1, "mm"),
        column_gap = unit(1, "mm"),
        na_col = "white",
        use_raster = T)
dev.off()


#2B 2
mean.peak.region <- getGroupSE(ArchRProj = satac, 
                               useMatrix = 'PeakMatrix', 
                               groupBy = 'regions', 
                               divideN = T,
                               scaleTo = 10^4)

mean.peak.region.mat <- mean.peak.region@assays@data$PeakMatrix

region.cor.mat <- cor(mean.peak.region.mat)
region.cor.mat[upper.tri(region.cor.mat)] <- NA

pdf('manuscript_figures/Fig2B_2_region_corr.pdf', width=12, height=10)
Heatmap(matrix = region.cor.mat, cluster_rows = F, 
        cluster_columns = F, show_row_names = T, row_names_side = "left",
        name = "Pearson",
        col = circlize::colorRamp2(c(0, 0.65, 1), c("lightblue","white","firebrick")), 
        row_title = NULL,
        column_title = NULL,
        row_gap = unit(1, "mm"),
        column_gap = unit(1, "mm"),
        na_col = "white",
        use_raster = T)
dev.off()

# 2B 3
encode.files <- list.files('ENCODE/DNase', pattern = "*.bed.gz", full.names = T)
encode.beds <- lapply(encode.files, function(x){ 
    encode <- vroom::vroom(x, col_names = F); 
    encode.gr <- GRanges(seqnames = encode$X1, ranges = IRanges(start = encode$X2, end = encode$X3));
    encode.gr
    }
)
names(encode.beds) <- sub("_.*", "", basename(encode.files))

getPeakEnrichment <- function(gr, encode.gr){
    
    encode_size <- sum(width(encode.gr)) 
    peak.set.size <- sum(width(gr))
    peak.overlaps.size <- sum(width(intersect(gr, encode.gr)))
    genome.size <- 3e9
    
    enrich <- log2((peak.overlaps.size / peak.set.size) / (encode_size / genome.size))
    enrich
}

celltypes <- names(markers)
enrich.result <- list()
for( c in celltypes){
    enrich.result[[c]] <- sapply(encode.beds, function(x){ getPeakEnrichment(markers[[c]], x) })
}
enrich.df <- as.data.frame(enrich.result)
enrich.df <- enrich.df[c("HeartLV","HeartRV","CardiacMuscle","DermisEndothelialBlood","CardiacFibroblast","NeuralStemCell","Tcell","Bcell","CD14Monocyte"),]
enrich.df <- enrich.df[,c("Cardiomyocyte","Smooth.Muscle","Pericyte","Endothelial","Fibroblast","Neuronal","Lymphoid","Myeloid")]
enrich.df.zscore <- enrich.df - rowMeans(enrich.df)

pdf('manuscript_figures/Fig2b_3_encode_enrichment.pdf', width=12, height=10)

Heatmap(matrix = enrich.df, cluster_rows = F, 
        cluster_columns = F, show_row_names = T, row_names_side = "left",
        name = "log2 Enrichment",
        col = circlize::colorRamp2(c(0,3.5,5.5), c("lightblue","white","firebrick")), 
        row_title = "ENCODE DNase",
        column_title = "Cell Types",
        row_gap = unit(1, "mm"),
        column_gap = unit(1, "mm"),
        na_col = "white",
        rect_gp = gpar(col = "black", lwd = 0.5),
        use_raster = F)

dev.off()

# 2C
# motif enrichment plot
archrPeaks <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_markerPeaks.rds')
archr_motifUp <- suppressMessages(peakAnnoEnrichment(seMarker = archrPeaks, 
                                                     ArchRProj = satac, 
                                                     peakAnnotation = "Vierstra", 
                                                     cutOff = "FDR <= 0.01 & Log2FC >= 1"))
TF.fdr.mat <- archr_motifUp@assays@data$mlog10Padj
pvalue.mat <- archr_motifUp@assays@data$mlog10p
    
tf.info <- data.frame(tf.names = rownames(TF.fdr.mat)) %>% 
    mutate(tf.prefix = sub('#.*', '', tf.names)) %>% 
    mutate(tf.cluster = sub(':.*', '', sub('.*#', '', tf.names))) %>%
    mutate(max_enrichment = apply(X = TF.fdr.mat, MARGIN = 1, FUN = function(x){max(x)})) %>%
    group_by(tf.cluster) %>% filter(max_enrichment == max(max_enrichment))

TF.fdr.mat <- TF.fdr.mat[tf.info$tf.names,]
pvalue.mat <- pvalue.mat[tf.info$tf.names,]

rownames(TF.fdr.mat) <- tf.info$tf.prefix
rownames(pvalue.mat) <- tf.info$tf.prefix

ct.specific.tf <- list()
ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")
for(t in ideal_order){
    curr <- TF.fdr.mat[TF.fdr.mat[,t] > 2,]
    ct.specific.tf[[t]] <- rownames(curr)[order(curr[t], decreasing = T)]
}
ct.motifs <- unlist(ct.specific.tf, use.names = F)
ct.motifs <- ct.motifs[!duplicated(ct.motifs)]
enrich.df <- pvalue.mat[ct.motifs, ideal_order]

tfs.to.highlight <- unlist(lapply(ct.specific.tf, function(x){sub('_.*','',x[1:3])}), use.names = F)
rownames.clean <- sub('_.*', '', rownames(enrich.df))

pdf('manuscript_figures/Fig2C_TF_enrich.pdf', width = 12, height=10)
rAnno <- rowAnnotation(n = anno_mark(at = match(tfs.to.highlight, rownames.clean), labels = tfs.to.highlight, labels_gp = gpar(fontsize = 22)))
Heatmap(matrix = enrich.df, cluster_rows = F, cluster_columns = F, show_column_names = T, show_row_names = F,
        col = circlize::colorRamp2(c(0,100), c("white","firebrick")), right_annotation = rAnno,
        row_title = NULL, border = T, column_split = 1:ncol(enrich.df), column_gap = unit(0, "mm"),
        column_title = NULL,
        name = "Motif Enrichment (-log10 pval)", 
        use_raster = T)
dev.off()


# 2D
# Expression of TFs
ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")

srna <- readRDS('seurat/Heart_RNA_Processed_Combined_NoAtrium.rds')
satac <- loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/')

motif.mat.obj <- getGroupSE(ArchRProj = satac, useMatrix = "VierstraMatrix", groupBy = "CellTypes", divideN = T, scaleTo = NULL)
normal.exp <- Seurat::AverageExpression(object = srna)

motif.mat <- motif.mat.obj@assays@data$VierstraMatrix
rowd <- rowData(motif.mat.obj)
motif.mat <- motif.mat[as.character(rowd$seqnames) == "deviations",]
rownames(motif.mat) <- rowd$name[as.character(rowd$seqnames) == "deviations"]

tf.info <- data.frame(tf.names = rownames(motif.mat)) %>% 
    mutate(tf.prefix = sub('#.*', '', tf.names)) %>% 
    mutate(tf.cluster = sub('/.*','',sub(':.*', '', sub('.*#', '', tf.names)))) %>%
    mutate(gene_name = sub('_.*', '', tf.names))

normal.exp.mat <- log2(as.matrix(normal.exp$RNA) + 1)
normal.exp.mat <- sweep(x = normal.exp.mat - rowMeans(normal.exp.mat), MARGIN = 1, STATS = matrixStats::rowSds(normal.exp.mat), FUN = '/')

TFs.enriched.in.RNA <- intersect(tf.info$gene_name, rownames(normal.exp.mat))

tf.info$In <- tf.info$gene_name %in% TFs.enriched.in.RNA
tf.info <- tf.info[tf.info$In,] 

tf.normal.exp <- as.data.frame(normal.exp.mat[tf.info$gene_name, ideal_order])
motif.enrich.df <- motif.mat[tf.info$tf.names, ideal_order]
motif.enrich.df <- sweep(x = motif.enrich.df - rowMeans(motif.enrich.df), MARGIN = 1, STATS = matrixStats::rowSds(motif.enrich.df), FUN = '/')

tf.info$RNA_gene_name <- rownames(tf.normal.exp)

### CORRELATION
corr.res <- rep(0, nrow(tf.normal.exp))
for(i in 1:nrow(tf.normal.exp)){
    corr.res[i] <- cor(as.numeric(tf.normal.exp[i,]), as.numeric(motif.enrich.df[i,]))
}

tf.info$cor <- corr.res
corr.cut <- 0.5
tf.info <- tf.info[tf.info$cor > corr.cut,]
tf.info <- tf.info %>% group_by(tf.cluster) %>% arrange(desc(cor)) %>% slice(1)

motif.enrich.df.cut <- motif.enrich.df[tf.info$tf.names,]
tf.normal.exp.cut <- tf.normal.exp[tf.info$RNA_gene_name,]

set.seed(100)
nClust = 8
row.clust.res <- kmeans(motif.enrich.df.cut, centers = nClust)
labs <- row.clust.res$cluster
rAnno <- rowAnnotation(cluster = factor(sort(labs)))
Heatmap(motif.enrich.df.cut[names(sort(labs)),], cluster_rows = F, cluster_columns = F, left_annotation = rAnno, col = circlize::colorRamp2(c(-2,0,2), c("lightblue","white","firebrick")))

lab_map <- data.frame(num = 1:nClust,
                      celltype = c("mix","Lymphoid1", "Lymphoid2", "Neuronal", "Myeloid1", "Cardiomyocyte", "Myeloid2","Endothelial"),
                      stringsAsFactors = F)

labs.named <- plyr::mapvalues(x = labs, from = lab_map$num, to = lab_map$celltype)
labs.named <- factor(labs.named, levels = c("Cardiomyocyte","mix","Pericyte","Endothelial","Fibroblast","Neuronal","Lymphoid1","Lymphoid2","Myeloid2","Myeloid1"))
labs.order <- order(labs.named)
labs.named <- labs.named[labs.order]

motif.enrich.df.cut <- motif.enrich.df.cut[labs.order, ideal_order]
tf.normal.exp.cut <- tf.normal.exp.cut[labs.order, ideal_order]
tf.info <- tf.info[labs.order,]

#tf.info$cellTypeOrdered <- labs.named
#tf.info2 <- tf.info %>% arrange(cellTypeOrdered, tf.cluster)
#motif.enrich.df.cut <- motif.enrich.df.cut[tf.info2$tf.names, ideal_order]

#tf.midpts <- tf.info2 %>% group_by(cellTypeOrdered, tf.cluster) %>% summarise(midpt  = floor(mean(index))) %>% arrange(midpt)
#tf.info2 <- tf.info2 %>% group_by(cellTypeOrdered, tf.cluster) %>% mutate(tf.cluster.binary = cur_group_id() %% 2)

rAnno <- HeatmapAnnotation(which = "column",
                           cluster_names = anno_mark(at = 1:nrow(tf.info),
                                                     labels = tf.info$gene_name,
                                                     side = "top"),
                           r = anno_barplot(x = tf.info$cor),
                           show_legend = F)

p1 <- Heatmap(matrix = t(motif.enrich.df.cut),
              cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = T, 
              col = circlize::colorRamp2(c(-2,0.1,2), c("lightblue","white","firebrick")), show_row_dend = F, show_column_dend = F,
              row_title = 'Motif', rect_gp = gpar(col = "black", lwd = 0.5),
              column_title = NULL,
              name = "Row-Normalized Motif Access", show_heatmap_legend = F,
              use_raster = T)

p2 <- Heatmap(matrix = t(tf.normal.exp.cut),
              cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = T,
              col = circlize::colorRamp2(c(-2,0.1,2), c("lightblue","white","firebrick")), 
              top_annotation = rAnno,  rect_gp = gpar(col = "black", lwd = 0.5),
              row_title = 'Expression', 
              show_heatmap_legend = T,
              column_title = NULL,
              name = "Row-Normalized Expression",
              use_raster = T)

pdf('manuscript_figures/figure2/Fig2D_TFMotif_Expression_Heatmap.pdf', width=14, height=7)
p2 %v% p1
dev.off()

