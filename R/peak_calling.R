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
  projHeart <- addPeakMatrix(projHeart, force = T)
  projHeart <- addMotifAnnotations(ArchRProj = projHeart, motifSet = "cisbp", name = "Motif")
  
  # cell-type specific peaks
  markersPeaks <- getMarkerFeatures(ArchRProj = projHeart, 
                                    useMatrix = "PeakMatrix", 
                                    groupBy = "CellTypes", 
                                    bias = c("TSSEnrichment", "log10(nFrags)"), 
                                    test)
  
  saveRDS(markersPeaks, paste0(archr_project_path,'/PeakCalls/DA_markerPeaks.rds'))
  
  markers <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = T)
  saveRDS(markers, file = paste0(archr_project_path,'/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds'))
  
  saveArchRProject(projHeart) 
  
  # shared peaks
  # ASSUMING all peaks have non-zero mean accessibility in each cell-type! Or else we do this for each cell-type iteratively
  peak.info <- getPeakSet(projHeart)
  types <- names(markers)
  for(t in types){
    curr.gr <- markers[[t]]
    hits <- GenomicRanges::findOverlaps(query = peak.info, subject = curr.gr)
    peak.info <- peak.info[-subjectHits(hits),]
  }
  
  saveRDS(peak.info, file = paste0(archr_project_path, '/PeakCalls/Vent_CM_Shared_Peaks.gr.rds'))
  
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

