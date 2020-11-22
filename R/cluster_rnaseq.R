library(Seurat)
library(harmony)
library(DoubletFinder)
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

# Load raw seurat object. Filter nUMI, and remove doublets per individual (little to no batch effects within individuals)
heart_rna <- readRDS('seurat/Heart_Seurat_RNA_all_samples_raw.RDS')
heart_rna <- base::subset(heart_rna, nCount_RNA > 1000 & nCount_RNA < 10000 & percent.mito < 10)
heart_rna_list <- SplitObject(heart_rna, split.by = 'individual')
for(i in 1:length(heart_rna_list)){
  heart_rna_list[[i]] <- NormalizeData(heart_rna_list[[i]], verbose = T)
  heart_rna_list[[i]] <- FindVariableFeatures(heart_rna_list[[i]], selection.method = "vst", 
                                    nfeatures = 2000, verbose = T)
  heart_rna_list[[i]] <- ScaleData(heart_rna_list[[i]])
  heart_rna_list[[i]] <- RunPCA(heart_rna_list[[i]], verbose = T, npcs = 100)
  
  heart_rna_list[[i]] <- doubletFinder_v3(seu = heart_rna_list[[i]], 
                                          PCs = 1:10, 
                                          pN = 0.015, 
                                          pK = 0.005, 
                                          nExp = round(0.075*nrow(heart_rna@meta.data)))
  n <- colnames(heart_rna_list[[i]]@meta.data)
  doubletCol <- n[startsWith("DF.classifications", x = n)]
  heart_rna_list[[i]]$doubletStatus <- heart_rna_list[[i]][[doubletCol]]
  heart_rna_list[[i]] <- base::subset(heart_rna_list[[i]], doubletStatus == "Singlet")
}

heart_rna_filtered <- merge(heart_rna_list$`03231`, c(heart_rna_list$`02207`, heart_rna_list$`02336`))

p <- VlnPlot(heart_rna_filtered, features = c("nCount_RNA","nFeature_RNA"), group.by = "individual", pt.size = 0.05)
ggsave(filename = 'seurat/plots/rna_qc_vln.png', plot = p, width = 10, height = 8)

saveRDS(heart_rna_filtered, file = 'seurat/Heart_Seurat_RNA_all_samples_Filtered.rds')

heart_rna_filtered <- NormalizeData(heart_rna_filtered, verbose = T)
heart_rna_filtered <- FindVariableFeatures(heart_rna_filtered, selection.method = "vst", 
                                             nfeatures = 2000, verbose = T)
heart_rna_filtered <- ScaleData(heart_rna_filtered)
heart_rna_filtered <- RunPCA(heart_rna_filtered, verbose = T)

heart_rna_filtered <- RunHarmony(heart_rna_filtered, group.by.vars = 'individual')

heart_rna_filtered <- RunUMAP(object = heart_rna_filtered, dims = 1:30, min.dist = 0.4, reduction = 'harmony')
heart_rna_filtered <- FindNeighbors(heart_rna_filtered, dims = 1:30, reduction ='harmony')
heart_rna_filtered <- FindClusters(object = heart_rna_filtered, resolution=0.2)

DimPlot(heart_rna_filtered, group.by = 'individual', pt.size = 0.01) 
DimPlot(heart_rna_filtered, label=T) 

markers <- c("TNNT2","MYBPC3","MYH7","NPPA","RGS5","ABCC9","MYH11","TAGLN","DCN","PDGFRA","PECAM1","VWF","PLP1","CD8A","LCK","CD14","FOLR2")
p <- FeaturePlot(heart_rna_filtered, features = markers)
ggsave(filename = paste0('seurat/plots/umap_combined_genes.png'), plot = p, width = 20, height=12, dpi = 300)

heart_rna_filtered <- RenameIdents(heart_rna_filtered, 
                                       '0'='Endothelial', 
                                       '1'='Fibroblast',
                                       '2'='Vent. CM',
                                       '3'='Pericyte',
                                       '4'='Myeloid',
                                       '5'='Smooth Muscle',
                                       '6'='Lymphoid',
                                       '7'='Endothelial',
                                       '8'='Neuronal',
                                       '9'='Endothelial',
                                       '10'='Neuronal',
                                       '11'='Atrial CM',
                                       '12'='Fibroblast')

palette <- setNames(brewer.pal(9, name="Set3"), sort(levels(Idents(heart_rna_filtered))))

p <- DimPlot(heart_rna_filtered, pt.size = 0.05, label = T, cols=palette) + ggtitle("scRNA-seq - 67114 Nuclei") + theme_set(theme_gray()) + ggClean() 
ggsave(filename = paste0('seurat/plots/umap_combined_celltypes.png'), plot = p, width = 8, height=6, dpi = 300)

saveRDS(heart_rna_filtered, file = "seurat/Heart_RNA_Processed_Combined.rds")
saveRDS(palette, 'notebooks/palette.rds')




