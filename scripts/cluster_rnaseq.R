library(Seurat)
library(future)
library(harmony)
library(ggplot2)
library(RColorBrewer)
palette <- brewer.pal(12, name="Set3")

ggClean <- function(){
  theme_bw() + 
    theme(text = element_text(size=18),
          legend.position = "right",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1))
}

setwd('/project2/gca/aselewa/heart_atlas_project/')

heart_rna <- readRDS('seurat/Heart_Seurat_RNA_all_samples_raw.RDS')
heart_rna <- base::subset(heart_rna, subset = nCount_RNA > 1000 & nCount_RNA < 20000 & percent.mito < 10)

p <- VlnPlot(heart_rna, features = c("nCount_RNA","nFeature_RNA","percent.mito"), group.by = "individual", pt.size = 0.05)
ggsave(filename = 'seurat/plots/rna_qc_vln.png', plot = p, width = 10, height = 8)

heart_rna_list <- SplitObject(heart_rna, split.by =  "individual")

for (i in 1:length(heart_rna_list)) {
  heart_rna_list[[i]] <- NormalizeData(heart_rna_list[[i]], verbose = T)
  heart_rna_list[[i]] <- FindVariableFeatures(heart_rna_list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = T)
  heart_rna_list[[i]] <- ScaleData(heart_rna_list[[i]])
  heart_rna_list[[i]] <- RunPCA(heart_rna_list[[i]], verbose = T)
}

for(i in 1:length(1:length(heart_rna_list))){
  heart_rna_list[[i]] <- FindNeighbors(heart_rna_list[[i]], reduction = "pca")
  heart_rna_list[[i]] <- RunUMAP(heart_rna_list[[i]], reduction = "pca", dims = 1:10, min.dist = 0.3)
  heart_rna_list[[i]] <- FindClusters(heart_rna_list[[i]], resolution = 0.2)
}

heart_rna_list$`03231` <- RenameIdents(heart_rna_list$`03231`, '0'='Fibroblast', '1'='Endothelial','2'='Vent. CM','3'='Pericyte','4'='Vent. CM','5'='Myeloid','6'='Endothelial','7'='Atrial CM','8'='Lymphoid','9'='Neuronal','10'='Vent. CM','11'='DCN+ Endothelial','12'='Vent. CM')

p1 <- DimPlot(heart_rna_list$`03231`, group.by = 'CellTypes', pt.size = 0.05, label = T, cols=palette) + ggtitle("scRNA-seq") + theme_set(theme_gray()) + ggClean() 
ggsave(filename = 'seurat/plots/umap_03231_celltypes.png', plot = p1, width = 8, height=6, dpi = 300)

saveRDS(heart_rna_list, file = "seurat/Heart_RNA_list_processed.rds")

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


# heart markers
markers <- c("TTN","MYBPC3","MYH7","NPPA","RGS5","ABCC9","MYH11","TAGLN","DCN","PDGFRA","PECAM1","VWF","PLP1","CD8A","LCK","CD14","FOLR2")

DefaultAssay(heart_rna_int) <- "RNA"
feature_plots <- FeaturePlot(heart_rna_int, features = markers)
ggsave(filename = 'seurat/plots/umap_CM_markers.png', plot = feature_plots, width = 20, height=15, dpi = 150)

# transfer labels
curr.idents <- as.factor(0:14)
new.idents <- c("Endothelial","Vent. CM","Fibroblast","Pericyte","Myeloid","Lymphoid","Neuronal","Pericyte","Endothelial","Vent. CM","Atrial CM","Other","Other 2","Myeloid","Other 3")
cluster.ids <- plyr::mapvalues(heart_rna_int$seurat_clusters, from = curr.idents, to = new.idents)
cluster.ids <- factor(cluster.ids, levels = unique(new.idents))
heart_rna_int$CellTypes <- cluster.ids

saveRDS(heart_rna_int, file = "seurat/Heart_RNA_integrated_processed.rds")

palette <- brewer.pal(11, name="Set3")
p1 <- DimPlot(heart_rna_int, group.by = 'CellTypes', pt.size = 0.05, label = T, cols=palette) + ggtitle("scRNA-seq") + theme_set(theme_gray()) + ggClean() 
ggsave(filename = 'seurat/plots/umap_celltypes.png', plot = p1, width = 8, height=6, dpi = 300)


# break down by cell type
freq <- table(heart_rna_int$CellTypes)/length(heart_rna$CellTypes)
type.freq <- data.frame(freq=as.numeric(freq), celltype=factor(names(freq), levels = levels(cluster.ids)))
p <- ggplot(type.freq, aes(x=celltype, y=freq, fill=celltype)) + 
  geom_bar(position = "dodge", stat="identity", color="black") + 
  ggClean() +
  scale_fill_manual(values = palette) +
  xlab("") +
  ylab("Prop. Cells")

ggsave(filename = 'seurat/plots/rna_prop_cell_type.png', p, dpi=300)



