library(ArchR)
library(ggplot2)
library(tidyverse)

addArchRThreads(threads = 8) 
setwd('/project2/gca/aselewa/heart_atlas_project/')
macs2 <- '/project2/gca/software/miniconda3/bin/macs2'

source('R/analysis_utils.R')

archr_project_path <- 'ArchR/ArchR_heart/'

projHeart <- loadArchRProject(archr_project_path)

projHeart <- addGroupCoverages(ArchRProj = projHeart, groupBy = "CellTypes", force = T, maxCells = 10000)
projHeart <- addReproduciblePeakSet(ArchRProj = projHeart, groupBy = "CellTypes", pathToMacs2 = macs2, cutOff = 0.01, verbose = T)
projHeart <- addPeakMatrix(projHeart, force = T)

# cell-type specific peaks
markersPeaks <- getMarkerFeatures(ArchRProj = projHeart, 
                                  useMatrix = "PeakMatrix", 
                                  groupBy = "CellTypes", 
                                  bias = c("TSSEnrichment", "log10(nFrags)"))

saveRDS(markersPeaks, paste0(archr_project_path,'/PeakCalls/DA_markerPeaks.rds'))

markers <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = T)
saveRDS(markers, file = paste0(archr_project_path,'/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds'))

saveArchRProject(projHeart) 
  
projHeart <- addMotifAnnotations(ArchRProj = projHeart, motifSet = "cisbp", name = "Motif", force = T)
projHeart <- addDeviationsMatrix(ArchRProj = projHeart, peakAnnotation = "Motif", force = T)
satac <- addGeneIntegrationMatrix(
  ArchRProj = satac, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "harmony",
  seRNA = srna,
  addToArrow = TRUE,
  force= TRUE,
  groupRNA = "cellTypes",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

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

run_workflow('ArchR/ArchR_heart/')

snpmap <- readRDS('eQTL_enrich/metadata/hg38_SNP_map.gr.rds')
create_torus_annotations(snpmap, gr.list = gr.list)

