library(ArchR)
library(ggplot2)
library(tidyverse)

addArchRThreads(threads = 12) 
setwd('/project2/gca/aselewa/heart_atlas_project/')
macs2 <- '/project2/gca/software/miniconda3/bin/macs2'

source('R/analysis_utils.R')

archr_project_path <- 'ArchR/ArchR_heart_latest_noAtrium/'

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
saveRDS(markers, file = paste0(archr_project_path,'/PeakCalls/DA_MARKERS_FDRP_10_log2FC_0.5.rds'))

saveArchRProject(projHeart) 

# Motif Enrichment
motifPWMs <- readRDS("Vierstra_Motifs/Vierstra-Human-Motifs.rds")
tf.prefix <-  sub('_.*', '', names(motifPWMs))
motifPWMs <- motifPWMs[tf.prefix == toupper(tf.prefix)] #human only

projHeart <- addMotifAnnotations(projHeart, motifPWMs = motifPWMs, name = "Vierstra")
projHeart <- addDeviationsMatrix(ArchRProj = projHeart, peakAnnotation = "Vierstra", force = T)

# Co-accessibility
satac <- addCoAccessibility(ArchRProj = satac, reducedDims = 'harmony', maxDist = 1e6)

# BigWigs by cell-type
getGroupBW(ArchRProj = satac, groupBy = "CellTypes")

saveArchRProject(satac) 




