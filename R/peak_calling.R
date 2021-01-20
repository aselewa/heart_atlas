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
projHeart <- addGeneIntegrationMatrix(
  ArchRProj = projHeart, 
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

saveArchRProject(projHeart) 


