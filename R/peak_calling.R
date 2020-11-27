library(ArchR)
library(ggplot2)
library(tidyverse)

addArchRThreads(threads = 8) 
setwd('/project2/gca/aselewa/heart_atlas_project/')
macs2 <- '/project2/gca/software/miniconda3/bin/macs2'

source('R/analysis_utils.R')

run_workflow <- function(archr_project_path){
  
  projHeart <- loadArchRProject(archr_project_path)
  
  projHeart <- addGroupCoverages(ArchRProj = projHeart, groupBy = "CellTypes")
  projHeart <- addReproduciblePeakSet(ArchRProj = projHeart, groupBy = "CellTypes", pathToMacs2 = macs2, cutOff = 0.01, verbose = T)
  awprojHeart <- addPeakMatrix(projHeart, force = T)
  
  markersPeaks <- getMarkerFeatures(ArchRProj = projHeart, 
                                    useMatrix = "PeakMatrix", 
                                    groupBy = "CellTypes", 
                                    bias = c("TSSEnrichment", "log10(nFrags)"), 
                                    test)
  
  saveRDS(markersPeaks, paste0(archr_project_path,'/PeakCalls/DA_markerPeaks.rds'))
  
  markers <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = T)
  saveRDS(markers, file = paste0(archr_project_path,'/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds'))
  
  saveArchRProject(projHeart) 
  
  # unsupervised clustering of peaks
  peak.mat <- getMatrixFromProject(projHeart, useMatrix = "PeakMatrix", binarize = T)
  peak.info <- getPeakSet(projHeart)
  peak.mat <- peak.mat@assays@data$PeakMatrix
  collapsedCellTypes <- projHeart$CellTypes
  
  if("Neuronal" %in% projHeart$CellTypes){ # not many cells in general
    peak.mat <- peak.mat[, projHeart$CellTypes != "Neuronal"]
    collapsedCellTypes <- collapsedCellTypes[collapsedCellTypes != "Neuronal"]
  }
  collapsedCellTypes[collapsedCellTypes == "Myeloid" | collapsedCellTypes == "Lymphoid"] <- "Immune"
  collapsedCellTypes[collapsedCellTypes == "Endothelial" | collapsedCellTypes == "Pericyte"] <- "EndoPericyte"
  
  ideal.access <- getIdealAccess(cell_type_vec = collapsedCellTypes)
  classify.mat <- peakToClusterBatch(peak.mat = peak.mat, ideal.mat = ideal.access, chunk.size = 1e5)

  saveRDS(classify.mat, file = paste0(archr_project_path,'/PeakCalls/peak_to_cluster.mtx.rds'))

  peak.profile <- apply(classify.mat, 1, which.max)
  gr.list <- list()
  types <- unique(collapsedCellTypes)
  profiles <- generateBits(n = length(types))
  for(i in 1:length(profiles)){
    curr <- profiles[[i]]
    profile_name <- paste0(types[which(curr==1)], collapse = '_and_')
    gr.list[[profile_name]] <- peak.info[peak.profile == i]
  }
  saveRDS(gr.list, file = paste0(archr_project_path,'/PeakCalls/cluster_peaks.gr.rds'))
  
}

create_torus_annotations <- function(snpmap, gr.list){
  
  celltypes <- names(gr.list)
  for(c in celltypes){
    overlap <- IRanges::subsetByOverlaps(snpmap, gr.list[[c]])
    snpsIn <- unique(overlap$snp)
    if(length(snpsIn) > 0){
      annot <- snpmap@elementMetadata
      n <- gsub(pattern = " ", replacement = "", c)
      n <- gsub(pattern = "/", replacement = "", n)
      annot[,paste0(n,'_d')] <- ifelse(annot$snp %in% snpsIn, 1, 0)
      vroom::vroom_write(as_tibble(annot), path = paste0('eQTL_enrich/annotations/DA_peaks_',n,'_edgeR_peaks.txt.gz')) 
    }
  }
}

projects <- list.files(path = 'ArchR', pattern = "ArchR_project_*", full.names = )
for(p in projects){
  print(paste0('Running workflow on ',p,'...'))
  run_workflow(p)
}

snpmap <- readRDS('eQTL_enrich/metadata/hg38_SNP_map.gr.rds')
create_torus_annotations(snpmap, gr.list = gr.list)


#### Supervised clustering of peaks
# classify.mat <- readRDS('../ArchR/ArchR_project_03231/PeakCalls/peak_to_cluster_jaccard.mtx.rds')
# peak.profile <- apply(classify.mat, 1, which.max)
# 
# collapsedCellTypes <- satac$CellTypes
# collapsedCellTypes <- collapsedCellTypes[collapsedCellTypes != "Neuronal"]
# collapsedCellTypes[collapsedCellTypes == "Myeloid" | collapsedCellTypes == "Lymphoid"] <- "Myeloid/Lymphoid"
# collapsedCellTypes[collapsedCellTypes == "Endothelial" | collapsedCellTypes == "Pericyte"] <- "Endo/Pericyte"
# 
# numPeaks <- as.numeric(table(peak.profile))
# types <- unique(collapsedCellTypes)
# profiles <- generateBits(n = length(types))
# profile_names <- list()
# for(i in 1:length(profiles)){
#   curr <- profiles[[i]]
#   profile_names[[i]] <- types[which(curr==1)]
# }
# set_mtx <- t(list_to_matrix(profile_names))
# observ_list <- list()
# for (i in 1:nrow(set_mtx)){
#   observ_list[[i]] <- matrix(rep(set_mtx[i, ], each = numPeaks[i]), nrow = numPeaks[i])
# }
# observ_mat <- Reduce(rbind, observ_list)
# colnames(observ_mat) <- colnames(set_mtx)
# combat_mtx <- make_comb_mat(observ_mat)
# 
# ComplexHeatmap::UpSet(combat_mtx, comb_col = c("red","blue","purple","black")[comb_degree(combat_mtx)],
#                       top_annotation = HeatmapAnnotation("Peaks" = anno_barplot(comb_size(combat_mtx),
#                                                                                 ylim = c(0, 1e5),
#                                                                                 border = FALSE, 
#                                                                                 gp = gpar(fill = "black"), 
#                                                                                 height = unit(4, "cm")),
#                                                          annotation_name_side = "left", 
#                                                          annotation_name_rot = 90),
#                       right_annotation = NULL)

