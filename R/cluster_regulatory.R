library(ArchR)
library(ComplexHeatmap)
source('R/analysis_utils.R')
satac <- loadArchRProject('ArchR/ArchR_heart/')

cell.type.mat <- suppressMessages(getGroupSE(ArchRProj = satac, useMatrix = "PeakMatrix", groupBy = "CellTypes", scaleTo = 10^6, verbose = F))
type.mat <- cell.type.mat@assays@data$PeakMatrix
row.names(type.mat) <- GRToString(getPeakSet(satac))
type.mat.norm <- type.mat/rowSums(type.mat)

set.seed(6699)
K=9
kmeans.res <- kmeans(type.mat.norm, centers = K)
labs <- kmeans.res$cluster

ideal_order <- c("Cardiomyocyte","Fibroblast", "Pericyte", "Neuronal", "Myeloid","Endothelial", "Lymphoid", "All")
lab_map <- data.frame(num = 1:K,
                      celltype = c("Cardiomyocyte", "Fibroblast", "Myeloid", "Endothelial", "All", "Lymphoid", "Neuronal", "Pericyte"),
                      stringsAsFactors = F)

labs.named <- plyr::mapvalues(x = labs, from = lab_map$num, to = lab_map$celltype)
labs.named <- factor(labs.named, levels = ideal_order)

labs.sort <- sort(labs.named)
type.mat.sort <- type.mat[names(labs.sort), c("Vent. CM", "Fibroblast", "Pericyte", "Neuronal", "Myeloid", "Endothelial", "Lymphoid")]
type.mat.sort <- t(scale(t(type.mat.sort), center = T, scale = T))

markers <- readRDS('ArchR/ArchR_heart/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')
cm.peaks.str <- GRToString(markers$`Vent. CM`)

colslist <- RColorBrewer::brewer.pal(n = K, "Set3")
names(colslist) <- paste0("K",1:K)
ha <- rowAnnotation(cluster = paste0("K", as.numeric(labs.sort)), DA_peak=(names(labs.sort) %in% cm.peaks.str) , col = list(cluster = colslist))

pdf('kmeans_peaks_DA.pdf', width = 14, height=8)
ComplexHeatmap::Heatmap(matrix = type.mat.sort, 
                        cluster_rows = F, 
                        cluster_columns = F, 
                        show_row_names = F, 
                        name = "Row accessibility\nZ-score", left_annotation = ha)
dev.off()



