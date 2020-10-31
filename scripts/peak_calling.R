
library(ArchR)
library(ggplot2)
library(tidyverse)

addArchRThreads(threads = 8) 
setwd('/project2/gca/aselewa/heart_atlas_project/')
macs2 <- '/project2/gca/software/miniconda3/bin/macs2'

projHeart <- loadArchRProject('ArchR_project_03231/')

projHeart <- addGroupCoverages(ArchRProj = projHeart, groupBy = "CellTypes")
projHeart <- addReproduciblePeakSet(ArchRProj = projHeart, groupBy = "CellTypes", pathToMacs2 = macs2)

files <- list.files(path = 'ArchR_project_03231/PeakCalls/og_peakCalls/', pattern = '*', full.names = T)
for(f in files){
  x <- readRDS(f)
  x <- x[x$score  > 10, ]
  saveRDS(x, paste0('ArchR_project_03231/PeakCalls/',basename(f)))
}

projHeart <- addPeakMatrix(projHeart, force = T)

markersPeaks <- getMarkerFeatures(
  ArchRProj = projHeart, 
  useMatrix = "PeakMatrix", 
  groupBy = "CellTypes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveRDS(markersPeaks, 'ArchR_project_03231/PeakCalls/DA_markerPeaks.rds')

markers <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = T)
saveRDS(markers, file = 'ArchR_project_03231/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')

saveArchRProject(projHeart)

# prepare annotations for Torus enrichment analysis

snpmap <- readRDS('eQTL_enrich/metadata/hg38_SNP_map.gr.rds')
celltypes <- names(markers)
for(c in celltypes){
  overlap <- IRanges::subsetByOverlaps(snpmap, markers[[c]])
  snpsIn <- unique(overlap$snp)
  if(length(snpsIn) > 0){
    annot <- snpmap@elementMetadata
    n <- gsub(pattern = " ", replacement = "", c)
    annot[,paste0(n,'_d')] <- ifelse(annot$snp %in% snpsIn, 1, 0)
    vroom::vroom_write(as_tibble(annot), path = paste0('eQTL_enrich/annotations/DA_peaks_',n,'_donor_03231_logfc1_minScore10.txt.gz')) 
  }
}


