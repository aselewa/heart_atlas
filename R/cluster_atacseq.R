setwd('/project2/gca/aselewa/heart_atlas_project/')

library(ArchR)
source('scripts/analysis_utils.R')
palette <- readRDS('palette.rds')

ReduceAndCluster <- function(archr_project, var.features=25000, resolution=0.1, min.dist=0.3, batch_correct=NULL){
  archr_project <- addIterativeLSI(ArchRProj = archr_project, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2, varFeatures = var.features, force = TRUE)
  
  do_batch_correction <- !is.null(batch_correct)
  if(do_batch_correction){
    archr_project <- addHarmony(ArchRProj = archr_project, reducedDims = "IterativeLSI", name = "harmony", groupBy = batch_correct, sigma = 0.05, force = TRUE)
  }
  archr_project <- addClusters(input = archr_project, reducedDims = ifelse(do_batch_correction, "harmony","IterativeLSI"), resolution = resolution, force = TRUE)
  archr_project <- addUMAP(ArchRProj = archr_project, reducedDims = ifelse(do_batch_correction, "harmony","IterativeLSI"), force=TRUE, minDist = min.dist)
  
  return(archr_project)
}

run_workflow <- function(project_path){

  projHeart <- suppressMessages(loadArchRProject(path = project_path))
  projHeart <- addDoubletScores(input = projHeart, k = 10, knnMethod = "UMAP", LSIMethod = 1)
  projHeart <- ArchR::subsetCells(projHeart, cellNames = projHeart$cellNames[!is.na(projHeart$DoubletScore)])
  projHeart <- ArchR::subsetCells(projHeart, cellNames = projHeart$cellNames[projHeart$DoubletScore < 2])
  
  projHeart <- ReduceAndCluster(archr_project = projHeart, resolution = 0.45)
  projHeart <- addImputeWeights(projHeart)
  
  saveArchRProject(projHeart)
  
  markers <- c("TNNT2","MYBPC3","MYH7","NPPA","RGS5","ABCC9","MYH11","TAGLN","DCN","PDGFRA","PECAM1","VWF","PLP1","CD8A","LCK","CD14","FOLR2")
  p <- plotEmbedding(ArchRProj = projHeart, colorBy = "GeneScoreMatrix", name = markers,embedding = "UMAP")
  plotPDF(plotList = p, name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", ArchRProj = projHeart, addDOC = FALSE, width = 5, height = 5)
  
  p <- custom_archr_umap(archr_project = projHeart, group.by="Clusters", palette = brewer.pal(n=10, name="Set3"), pt.size=0.3, alpha=0.7, label=F, legend = T)
  ggsave(filename = paste0(project_path,"/Plots/umap-clusters.png"), plot = p, dpi=150, width=8, height=6)
  
  p <- getQCPlots(projHeart)
  ggsave(filename = paste0(project_path,"/Plots/VlnQCPlots.png"), plot = p, dpi=150, width=14, height=10)
}

# reduce dim and cluster for each project
projects <- list.files(pattern = "ArchR_project_*")
for(p in projects){
  print(paste0('Running workflow on ',p,'...'))
  run_workflow(p)
}

# Donor-specific cell type assignment - cant be made generic
projHeart <- loadArchRProject('ArchR_project_03231/')
custom_archr_umap(archr_project = projHeart, group.by="Clusters", palette = brewer.pal(n=10, name="Set3"), pt.size=0.3, alpha=0.7, label=T, legend = T)
cluster.ids <- paste0("C",1:9)
new.ids <- c("Vent. CM", "Vent. CM","Vent. CM","Myeloid","Neuronal","Lymphoid","Endothelial", "Fibroblast","Pericyte")
projHeart$CellTypes <- RenameIdentity(idents = projHeart$Clusters, from = cluster.ids, to = new.ids)
p <- custom_archr_umap(archr_project = projHeart, group.by="CellTypes", palette = palette, pt.size=0.3, alpha=0.7, label=T, legend = T)
ggsave(filename = "ArchR_project_03231/Plots/umap-CellTypes.png", plot = p, dpi=150, width=8, height=6)
saveArchRProject(projHeart)


# Donor-specific cell type assignment - cant be made generic
projHeart <- loadArchRProject('ArchR_project_02207/')
custom_archr_umap(archr_project = projHeart, group.by="Clusters", palette = brewer.pal(n=10, name="Set3"), pt.size=0.3, alpha=0.7, label=T, legend = T)
cluster.ids <- paste0("C",1:10)
new.ids <- c("Myeloid", "Lymphoid","Vent. CM","Vent. CM","Vent. CM","Endothelial","Fibroblast", "Fibroblast","Pericyte","Smooth Muscle")
projHeart$CellTypes <- RenameIdentity(idents = projHeart$Clusters, from = cluster.ids, to = new.ids)
p <- custom_archr_umap(archr_project = projHeart, group.by="CellTypes", palette = palette, pt.size=0.3, alpha=0.7, label=T, legend = T)
ggsave(filename = "ArchR_project_02207/Plots/umap-CellTypes.png", plot = p, dpi=150, width=8, height=6)
saveArchRProject(projHeart)


# Donor-specific cell type assignment - cant be made generic
projHeart <- loadArchRProject('ArchR_project_02336/')
custom_archr_umap(archr_project = projHeart, group.by="Clusters", palette = brewer.pal(n=10, name="Set3"), pt.size=0.3, alpha=0.7, label=T, legend = T)
cluster.ids <- paste0("C",1:7)
new.ids <- c("Myeloid", "Lymphoid","Vent. CM","Vent. CM","Endothelial","Fibroblast","Pericyte")
projHeart$CellTypes <- RenameIdentity(idents = projHeart$Clusters, from = cluster.ids, to = new.ids)
p <- custom_archr_umap(archr_project = projHeart, group.by="CellTypes", palette = palette, pt.size=0.3, alpha=0.7, label=T, legend = T)
ggsave(filename = "ArchR_project_02336/Plots/umap-CellTypes.png", plot = p, dpi=150, width=8, height=6)
saveArchRProject(projHeart)
