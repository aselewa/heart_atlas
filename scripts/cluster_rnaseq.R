library(Seurat)
library(future)
library(harmony)
library(ggplot2)
library(RColorBrewer)

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
heart_rna <- base::subset(heart_rna, nCount_RNA > 500 & nCount_RNA < 10000 & percent.mito < 10)

p <- VlnPlot(heart_rna, features = c("nCount_RNA","nFeature_RNA","percent.mito"), group.by = "individual", pt.size = 0.05)
ggsave(filename = 'seurat/plots/rna_qc_vln.png', plot = p, width = 10, height = 8)

heart_rna_list <- SplitObject(heart_rna, split.by =  "individual")
min.feats <- c(500, 500, 500)
max.feats <- c(2500, 1500, 2500)
for(i in 1:length(heart_rna_list)) {
  
  heart_rna_list[[i]] <- base::subset(heart_rna_list[[i]], subset = nFeature_RNA > min.feats[i] & nFeature_RNA < max.feats[i])
  
  heart_rna_list[[i]] <- NormalizeData(heart_rna_list[[i]], verbose = T)
  heart_rna_list[[i]] <- FindVariableFeatures(heart_rna_list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = T)
  heart_rna_list[[i]] <- ScaleData(heart_rna_list[[i]])
  heart_rna_list[[i]] <- RunPCA(heart_rna_list[[i]], verbose = T)
}

markers <- c("TNNT2","MYBPC3","MYH7","NPPA","RGS5","ABCC9","MYH11","TAGLN","DCN","PDGFRA","PECAM1","VWF","PLP1","CD8A","LCK","CD14","FOLR2")
res <- c(0.2, 0.2, 0.25)
for(i in 1:length(1:length(heart_rna_list))){
  heart_rna_list[[i]] <- FindNeighbors(heart_rna_list[[i]], reduction = "pca")
  heart_rna_list[[i]] <- RunUMAP(heart_rna_list[[i]], reduction = "pca", dims = 1:10, min.dist = 0.3)
  heart_rna_list[[i]] <- FindClusters(heart_rna_list[[i]], resolution = res[i])
  
  p <- FeaturePlot(heart_rna_list[[i]], features = markers)
  ggsave(filename = paste0('seurat/plots/umap_',names(heart_rna_list[i]),'_celltypes.png'), plot = p, width = 20, height=12, dpi = 300)
}


heart_rna_list$`02207` <- RenameIdents(heart_rna_list$`02207`, 
                                       '0'='Endothelial', 
                                       '1'='Pericyte',
                                       '2'='Fibroblast',
                                       '3'='Myeloid',
                                       '4'='Vent. CM',
                                       '5'='Lymphoid',
                                       '6'='Smooth Muscle',
                                       '7'='Neuronal',
                                       '8'='Endothelial',
                                       '9'='Lymphoid')

heart_rna_list$`02336` <- RenameIdents(heart_rna_list$`02336`,
                                       '0'='Endothelial',
                                       '1'='Pericyte',
                                       '2'='Fibroblast',
                                       '3'='Myeloid',
                                       '4'='Vent. CM',
                                       '5'='Endothelial',
                                       '6'='Smooth Muscle',
                                       '7'='Lymphoid',
                                       '8'='Neuronal',
                                       '9'='Other')

heart_rna_list$`03231` <- RenameIdents(heart_rna_list$`03231`, 
                                       '0'='Fibroblast', 
                                       '1'='Endothelial',
                                       '2'='Vent. CM',
                                       '3'='Pericyte',
                                       '4'='Endothelial',
                                       '5'='Myeloid',
                                       '6'='Vent. CM',
                                       '7'='Lymphoid',
                                       '8'='Neuronal',
                                       '9'='Smooth Muscle',
                                       '10'='Atrial CM',
                                       '11'='Fibroblast')

palette <- setNames(brewer.pal(10, name="Set3"), c(levels(Idents(heart_rna_list$`03231`)),"Other"))

for(i in 1:length(heart_rna_list)){
  p <- DimPlot(heart_rna_list[[i]], pt.size = 0.05, label = T, cols=palette) + ggtitle(paste0("scRNA-seq - ",names(heart_rna_list)[i])) + theme_set(theme_gray()) + ggClean() 
  ggsave(filename = paste0('seurat/plots/umap_',names(heart_rna_list)[i],'_celltypes.png'), plot = p, width = 8, height=6, dpi = 300)
}

saveRDS(heart_rna_list, file = "seurat/Heart_RNA_list_processed.rds")





