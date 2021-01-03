setwd('/project2/gca/aselewa/heart_atlas_project/')

library(ArchR)
source('R/analysis_utils.R')
palette <- readRDS('notebooks/palette.rds')

ReduceAndCluster <- function(archr_project, var.features=20000, min.dist=0.5, resolution=0.8, batch_correct=NULL){
  
  archr_project <- addIterativeLSI(ArchRProj = archr_project, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2, varFeatures = var.features, force = TRUE)
  
  do_batch_correction <- !is.null(batch_correct)
  if(do_batch_correction){
    archr_project <- addHarmony(ArchRProj = archr_project, reducedDims = "IterativeLSI", name = "harmony", groupBy = "individual", sigma=0.01, force = TRUE)
  }
  archr_project <- addClusters(input = archr_project, reducedDims = ifelse(do_batch_correction, "harmony","IterativeLSI"), resolution = resolution, force = TRUE)
  archr_project <- addUMAP(ArchRProj = archr_project, reducedDims = ifelse(do_batch_correction, "harmony","IterativeLSI"), force=TRUE, minDist = min.dist)
  
  return(archr_project)
}

run_workflow <- function(project_path, batch_correct=NULL){

  projHeart <- suppressMessages(loadArchRProject(path = project_path))
  projHeart <- subsetCells(ArchRProj = projHeart, cellNames = projHeart$cellNames[projHeart$DoubletEnrichment <= 1])
  projHeart <- ReduceAndCluster(archr_project = projHeart, batch_correct = batch_correct)
  projHeart <- addImputeWeights(projHeart)
  
  projHeart <- saveArchRProject(projHeart, load = T)

  markers <- c("TNNT2","MYBPC3","MYH7","NPPA","RGS5","ABCC9","MYH11","TAGLN","DCN","PDGFRA","PECAM1","VWF","PLP1","CD8A","LCK","CD14","FOLR2")
  p <- plotEmbedding(ArchRProj = projHeart, colorBy = "GeneScoreMatrix", name = markers,embedding = "UMAP")
  plotPDF(plotList = p, name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", ArchRProj = projHeart, addDOC = FALSE, width = 5, height = 5)
  
  p <- custom_archr_umap(archr_project = projHeart, pt.size=0.3, alpha=0.7, label=F, legend = T)
  ggsave(filename = paste0(project_path,"/Plots/umap-regions.png"), plot = p, dpi=150, width=8, height=6)
  
  p <- custom_archr_umap(archr_project = projHeart, group.by="Clusters", pt.size=0.3, alpha=0.7, label=T, legend = T)
  ggsave(filename = paste0(project_path,"/Plots/umap-Clusters.png"), plot = p, dpi=150, width=8, height=6)
  
  p <- getQCPlots(projHeart)
  ggsave(filename = paste0(project_path,"/Plots/VlnQCPlots.png"), plot = p, dpi=150, width=14, height=10)
  
  return(projHeart)
  
}

# reduce dim and cluster for each project
projects <- list.files(path = 'ArchR', pattern = "*_3", full.names = T)
for(p in projects){
  print(paste0('Running workflow on ',p,'...'))
  run_workflow(p)
}

# Donor-specific cell type assignment - cant be made generic
projHeart <- loadArchRProject('ArchR/ArchR_project_03231_3/')
custom_archr_umap(archr_project = projHeart, group.by="Clusters", palette = NULL, pt.size=0.3, alpha=0.7, label=T, legend = T)
cluster.ids <- paste0("C",1:10)
new.ids <- c("Myeloid", "Lymphoid","Endothelial","Pericyte","Vent. CM","Fibroblast","Fibroblast","Neuronal", "Vent. CM","Vent. CM")
projHeart$CellTypes <- RenameIdentity(idents = projHeart$Clusters, from = cluster.ids, to = new.ids)
p <- custom_archr_umap(archr_project = projHeart, group.by="CellTypes", palette = palette, pt.size=0.3, alpha=0.7, label=T, legend = T)
ggsave(filename = "ArchR/ArchR_project_03231_3/Plots/umap-CellTypes.png", plot = p, dpi=150, width=8, height=6)
saveArchRProject(projHeart)


# Donor-specific cell type assignment - cant be made generic
projHeart <- loadArchRProject('ArchR/ArchR_project_02207_3/')
custom_archr_umap(archr_project = projHeart, group.by="Clusters", palette = NULL, pt.size=0.3, alpha=0.7, label=T, legend = T)
cluster.ids <- paste0("C",1:13)
new.ids <- c("Vent. CM","Vent. CM","Vent. CM", "Endothelial", "Neuronal", "Pericyte","Smooth Muscle", "Fibroblast","Fibroblast","Fibroblast","Myeloid","Lymphoid","Myeloid")
projHeart$CellTypes <- RenameIdentity(idents = projHeart$Clusters, from = cluster.ids, to = new.ids)
p <- custom_archr_umap(archr_project = projHeart, group.by="CellTypes", palette = palette, pt.size=0.3, alpha=0.7, label=T, legend = T)
ggsave(filename = "ArchR/ArchR_project_02207_3/Plots/umap-CellTypes.png", plot = p, dpi=150, width=8, height=6)
saveArchRProject(projHeart)


# Donor-specific cell type assignment - cant be made generic
projHeart <- loadArchRProject('ArchR/ArchR_project_02336_3/')
custom_archr_umap(archr_project = projHeart, group.by="Clusters", palette = NULL, pt.size=0.3, alpha=0.7, label=T, legend = T)
cluster.ids <- paste0("C",1:8)
new.ids <- c("Vent. CM", "Vent. CM","Myeloid","Lymphoid","Endothelial","Pericyte","Fibroblast","Fibroblast")
projHeart$CellTypes <- RenameIdentity(idents = projHeart$Clusters, from = cluster.ids, to = new.ids)
p <- custom_archr_umap(archr_project = projHeart, group.by="CellTypes", palette = palette, pt.size=0.3, alpha=0.7, label=T, legend = T)
ggsave(filename = "ArchR/ArchR_project_02336_3/Plots/umap-CellTypes.png", plot = p, dpi=150, width=8, height=6)
saveArchRProject(projHeart)

## Combined project
p <- 'ArchR/ArchR_heart/'
run_workflow(p, batch_correct = TRUE)

projHeart <- loadArchRProject('ArchR/ArchR_heart/')
custom_archr_umap(archr_project = projHeart, group.by="Clusters", palette = NULL, pt.size=0.3, alpha=0.7, label=T, legend = T)
cluster.ids <- paste0("C",1:11)
new.ids <- c("Myeloid","Lymphoid","Vent. CM","Vent. CM","Vent. CM","Endothelial","Pericyte","Pericyte","Neuronal","Fibroblast","Fibroblast")
projHeart$CellTypes <- RenameIdentity(idents = projHeart$Clusters, from = cluster.ids, to = new.ids)
p <- custom_archr_umap(archr_project = projHeart, group.by="CellTypes", palette = palette, pt.size=0.3, alpha=0.7, label=T, legend = T)
ggsave(filename = "ArchR/ArchR_heart/Plots/umap-CellTypes.png", plot = p, dpi=150, width=8, height=6)
saveArchRProject(projHeart)
