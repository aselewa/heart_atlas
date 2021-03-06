library(ArchR)
library(RColorBrewer)
setwd('/project2/gca/aselewa/heart_atlas_project/')
addArchRThreads(threads = 8) 

# Integration of all RNA donors
heart_rna_list <- SplitObject(heart_rna, split.by =  "individual")
for (i in 1:length(heart_rna_list)) {
  heart_rna_list[[i]] <- NormalizeData(heart_rna_list[[i]], verbose = T)
}

int_featurres <- SelectIntegrationFeatures(object.list = heart_rna_list, nfeatures = 3000)
reference_dataset <- which(names(heart_rna_list) == "03231")

heart.anchors <- FindIntegrationAnchors(object.list = heart_rna_list, normalization.method = "LogNormalize", 
                                        anchor.features = int_featurres, reference = reference_dataset)

saveRDS(heart.anchors, file = 'seurat/Heart_anchors_integration.rds')

heart_rna_int <- IntegrateData(anchorset = heart.anchors, normalization.method = "LogNormalize")
heart_rna_int <- ScaleData(heart_rna_int)
heart_rna_int <- RunPCA(object = heart_rna_int, verbose = FALSE)
heart_rna_int <- RunUMAP(object = heart_rna_int, dims = 1:50, min.dist = 0.2)
heart_rna_int <- FindNeighbors(heart_rna_int, reduction = "pca", dims = 1:50)
heart_rna_int <- FindClusters(object = heart_rna_int, resolution=0.15)

# integrate RNA and ATAC
s.rna <- readRDS(file = 'seurat/Heart_RNA_list_processed.rds')
seRNA <- s.rna$`03231`
seRNA$CellTypes <- Seurat::Idents(seRNA)

projHeart <- loadArchRProject('ArchR_project_03231/')
projHeart <- addGeneIntegrationMatrix(ArchRProj = projHeart, useMatrix = 'GeneScoreMatrix', matrixName = 'GeneIntegrationMatrix' ,seRNA = seRNA, groupRNA = 'CellTypes', force = T)

colLookup <- readRDS('color_lookup_table.rds')
palette <- colLookup[sort(unique(projHeart$predictedGroup)), "cols"]
p <- custom_archr_umap(archr_project = projHeart, group.by="predictedGroup", palette = palette, pt.size=0.3, alpha=0.7, label=T, legend = T)
ggsave(filename = "ArchR_project_03231/Plots/umap-predicted_RNA_label.png", plot = p, dpi=150, width=8, height=6)

cell_type_agree <- c(sum(projHeart$CellTypes == projHeart$predictedGroup), sum(projHeart$CellTypes != projHeart$predictedGroup))
names(cell_type_agree) <- c('Agree','Disagree')
scores_agree <- projHeart$predictedScore[projHeart$CellTypes == projHeart$predictedGroup]
scores_disagree <-  projHeart$predictedScore[projHeart$CellTypes != projHeart$predictedGroup]

par(mfrow=c(1,2))
barplot(cell_type_agree)
boxplot(scores_agree, scores_disagree, names=c('Agree','Disagree'))


hist(projHeart$predictedScore, xlab='Max Predicted Score', main = '')
