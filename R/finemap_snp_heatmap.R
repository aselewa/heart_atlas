library(tidyverse)
library(ComplexHeatmap)
library(ArchR)
setwd('/project2/gca/aselewa/heart_atlas_project/')

satac <- loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/')

celltype.bw <- list.files('ArchR/ArchR_heart_latest_noAtrium/GroupBigWigs/CellTypes', pattern = 'Hg19*', full.names = T)
bw.list <- lapply(celltype.bw, function(x){rtracklayer::import(x)})

names <- c("Cardiomyocyte","Endothelial","Fibroblast","Lymphoid","Myeloid","Neuronal","Pericyte","Smooth Muscle")

finemap.res <- readr::read_csv('GWAS/aFib_Finemapped_minPIP_0.2_07022021.csv')
finemap.res <- finemap.res[finemap.res$PIP>0.5,]

finemap.res.gr <- GRanges(seqnames = finemap.res$chr, ranges = IRanges(start = finemap.res$`b37 bp`, end = finemap.res$`b37 bp`), snp = finemap.res$SNP)
seqlevelsStyle(finemap.res.gr) <- "UCSC"

snp.access.mat <- list()
for(i in 1:length(bw.list)){
    snp.access.mat[[i]] <- subsetByOverlaps(bw.list[[i]], finemap.res.gr)$score
}

snp.access.mat.df <- data.frame(matrix(unlist(snp.access.mat), ncol = length(snp.access.mat)))
colnames(snp.access.mat.df) <- names

rAnno <- HeatmapAnnotation(which = "column",
                           `H3K27ac` = factor(1*(!is.na(finemap.res$AdultHeart_H3K27ac))),
                           `Fetal DHS` = factor(1*(!is.na(finemap.res$FetalHeart_DNase))),
                           `PC-HiC or Coaccess` = factor(ifelse(yes = "pcHiC or Coacc", no = "Other", test = finemap.res$`Link Method`=="coacc" |  finemap.res$`Link Method`=="HiC")),
                           `TBX5/NKX2-5/GATA4 ChIP` = factor(1*(!is.na(finemap.res$FGT_ChIPseq))),
                           `TF Motif Break`= factor(1*(!is.na(finemap.res$motif_break_strong))), 
                           col = list(`H3K27ac` = c("0" = "mistyrose", "1" = "magenta"),
                                      `TF Binding` = c("0" = "slategray1", "1"="dodgerblue"),
                                      `Fetal DHS` = c("0" = "palegreen", "1" = "seagreen"),
                                      `PC-HiC or Coaccess` = c("pcHiC or Coacc" = "darkgoldenrod1", "Other"="lightyellow"),
                                      `TF Motif Break` = c("0" = "lightgrey", "1"="black")))

snp.access.mat.df <- log2((snp.access.mat.df)*100 + 1)
snp.access.mat.df[snp.access.mat.df > 5] <- 5

pdf('snp_access_heatmap.pdf', width=8, height=5)
p <- Heatmap(t(snp.access.mat.df), 
        col = circlize::colorRamp2(c(0,2,5), c("lightblue","white","firebrick")), 
        bottom_annotation = rAnno, 
        name = 'log2 accessibility',
        cluster_rows = T, 
        cluster_columns = T,
        height = unit(4, "cm"), 
        show_row_dend = F, 
        show_column_dend = F,
        heatmap_legend_param = list(direction = "horizontal"))
draw(p, heatmap_legend_side = "top")
dev.off()

        