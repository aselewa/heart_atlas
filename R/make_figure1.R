library(Seurat)
library(ArchR)
setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

srna <- readRDS('seurat/Heart_RNA_Processed_Combined_NoAtrium.rds')
satac <- suppressMessages(loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/', showLogo = F))


#### 1B -- cell type UMAPs

palette <- readRDS('notebooks/palette.rds')

srna$cellTypes <- Idents(srna)

rna.umap <- custom_dim_plot(srna, group.by = "cellTypes", label = T, legend = F, palette = palette) + ggtitle('snRNA-seq') + xlim(c(-15,15)) 
atac.umap <- custom_archr_plot(archr_project = satac, group.by = 'CellTypes', palette = palette, legend = F, label = T) + ggtitle('snATAC-seq') + xlim(c(-15,15)) 

rna.umap + atac.umap

### 1B -- barplots on side of umaps

donors <- unique(srna$individual)
freq <- list()
for(i in 1:length(donors)){
  curr.obj <- subset(srna, individual == donors[i])
  freq[[i]] <- table(Idents(curr.obj))/length(Idents(srna))
}
freqs.unlist <- unlist(freq)
type.freq <- data.frame(freq=as.numeric(freqs.unlist), celltype=factor(names(freqs.unlist)), type = rep(donors, each = length(freq[[1]])))
p1 <- ggplot(type.freq, aes(x=freq, y=celltype, fill = celltype, group=type)) + geom_bar(position = "stack", stat="identity", color="black") + 
  ggClean(rotate_axis = F) + scale_fill_manual(values = palette) + xlab("Prop. Cells") + ylab("") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  LegendOff() + ggtitle('scnNA-seq')

freq <- list()
for(i in 1:length(donors)){
  curr.obj <- subsetCells(ArchRProj = satac, cellNames = satac$cellNames[satac$individual == donors[i]])
  freq[[i]] <- table(curr.obj$CellTypes)/length(satac$cellNames)
}
freqs.unlist <- unlist(freq)

type.freq <- data.frame(freq=as.numeric(freqs.unlist), celltype=factor(names(freqs.unlist)), type = rep(donors, each = length(freq[[1]])))
p2 <- ggplot(type.freq, aes(x=freq, y=celltype, fill = celltype, group=type)) + geom_bar(position = "stack", stat="identity", color="black") + 
  ggClean(rotate_axis = F) + scale_fill_manual(values = palette) + xlab("Prop. Cells") + ylab("") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  LegendOff() + ggtitle('snATAC-seq')


p1 + p2


### 1C -- barplots on side of umaps
### This figure is complex - so i made it its own script: make_figure1C.R



### 1D -- % nuclei expressing

cell.type.exp <- AverageExpression(object = srna, slot = 'data', verbose = F) # exponentiates and returns data in non-logged space
cell.type.exp.mat <- log(1 + cell.type.exp$RNA)

markers <- c("TNNT2","VWF","DCN","CD8A","CD14","PLP1","RGS5","TAGLN")

counts <- srna@assays$RNA@counts
counts <- counts[markers,] > 0
pct.express <- c()
srna$cellTypes <- Idents(srna)
for(type in unique(Idents(srna))){
  currobj <- subset(srna, cellTypes == type)
  counts <- currobj@assays$RNA@counts
  counts <- counts[markers,] > 0
  pct.express <- c(pct.express, rowMeans(counts))
}

cell.type.exp.mat <- cell.type.exp.mat[markers,]
pct.express.df <- data.frame('pExpressing' = unname(pct.express), 
                             genes = factor(names(pct.express), levels = rev(markers)),
                             celltypes = factor(rep(unique(Idents(srna)), each = length(markers)), levels = sort(levels(Idents(srna)))),
                             mean_exp = unlist(cell.type.exp.mat, use.names = F))

ggplot(pct.express.df, aes(x = celltypes, y = genes, color = mean_exp, size = pExpressing)) + 
  ggrastr::rasterise(geom_point(), dpi=300) + 
  ggClean(rotate_axis = T) +
  scale_color_gradientn(colors = c("lightblue","yellow","firebrick"))


