library(ArchR)
library(edgeR)

setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

satac <- loadArchRProject('ArchR/ArchR_heart/')
peak.info <- getPeakSet(satac)
peak.names <- GRToString(peak.info)

agg.insert.mat <- getGroupSE(ArchRProj = satac, useMatrix = "PeakMatrix", groupBy = "CellTypes", divideN = F)
agg.insert.mat <- agg.insert.mat@assays@data$PeakMatrix
row.names(agg.insert.mat) <- peak.names

cell.type.mat <- log2(10^6*t(agg.insert.mat)/colSums(agg.insert.mat) + 1)

enhancer.umap <- uwot::umap(X = t(cell.type.mat)[toViz,], n_neighbors = 50, n_components = 2, min_dist = 0.4)

kres <- stats::kmeans(x = cell.type.mat, centers = 5)
labels <- kres$cluster

plot(enhancer.umap[,1], enhancer.umap[,2], pch=16)
