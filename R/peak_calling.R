library(ArchR)
library(ggplot2)
library(tidyverse)

addArchRThreads(threads = 8) 
setwd('/project2/gca/aselewa/heart_atlas_project/')
macs2 <- '/project2/gca/software/miniconda3/bin/macs2'

run_workflow <- function(archr_project_path){
  
  projHeart <- loadArchRProject(archr_project_path)
  
  projHeart <- addGroupCoverages(ArchRProj = projHeart, groupBy = "CellTypes")
  projHeart <- addReproduciblePeakSet(ArchRProj = projHeart, groupBy = "CellTypes", pathToMacs2 = macs2)
  projHeart <- addPeakMatrix(projHeart, force = T)
  
  markersPeaks <- getMarkerFeatures(ArchRProj = projHeart, useMatrix = "PeakMatrix", groupBy = "CellTypes", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
  
  saveRDS(markersPeaks, paste0(archr_project_path,'/PeakCalls/DA_markerPeaks.rds'))
  
  markers <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = T)
  saveRDS(markers, file = paste0(archr_project_path,'/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds'))
  
  saveArchRProject(projHeart) 
}

create_torus_annotations <- function(snpmap, archr_project_path){
  
  DA_markers <- readRDS(paste0(archr_project_path,'/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds'))
  donor <- strsplit(archr_project_path, split = '_')[[1]][3]
  celltypes <- names(DA_markers)
  for(c in celltypes){
    overlap <- IRanges::subsetByOverlaps(snpmap, DA_markers[[c]])
    snpsIn <- unique(overlap$snp)
    if(length(snpsIn) > 0){
      annot <- snpmap@elementMetadata
      n <- gsub(pattern = " ", replacement = "", c)
      annot[,paste0(n,'_d')] <- ifelse(annot$snp %in% snpsIn, 1, 0)
      vroom::vroom_write(as_tibble(annot), path = paste0('eQTL_enrich/annotations/DA_peaks_',n,'_donor_',donor,'_logfc1_minScore10.txt.gz')) 
    }
  }
}

projects <- list.files(pattern = "ArchR_project_*")
for(p in projects){
  print(paste0('Running workflow on ',p,'...'))
  run_workflow(p)
}

snpmap <- readRDS('eQTL_enrich/metadata/hg38_SNP_map.gr.rds')
for(p in projects){
  create_torus_annotations(snpmap, p)
}
