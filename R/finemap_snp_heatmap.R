library(tidyverse)
library(ComplexHeatmap)
library(ArchR)

satac <- loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/')

celltype.bw <- list.files('ArchR/ArchR_heart_latest_noAtrium/GroupBigWigs/CellTypes', pattern = 'Hg19*', full.names = T)
bw.list <- lapply(celltype.bw, function(x){rtracklayer::import(x)})


finemap.res <- readr::read_csv('GWAS/aFib_Finemapped_minPIP_0.2_06152021.csv')
finemap.res <- finemap.res[finemap.res$PIP>0.5,]

finemap.res.gr <- GRanges(seqnames = finemap.res$chr, ranges = IRanges(start = finemap.res$`b37 bp`, end = finemap.res$`b37 bp`), snp = finemap.res$SNP)
seqlevelsStyle(finemap.res.gr) <- "UCSC"

snp.access.mat <- list()
for(i in 1:length(bw.list)){
    snp.access.mat[[i]] <- subsetByOverlaps(bw.list[[i]], finemap.res.gr)$score
}

snp.access.mat.df <- data.frame(matrix(unlist(snp.access.mat), ncol = length(snp.access.mat)))

rAnno <- HeatmapAnnotation(which = "row",
                           h3k27ac = factor(1*(!is.na(finemap.res$H3K27ac_heart_concat))),
                           chipseq = factor(1*(!is.na(finemap.res$FGT_ChIP_lifted_from_mm10))))

snp.access.mat.df <- snp.access.mat.df/max(snp.access.mat.df)

Heatmap(snp.access.mat.df, 
        col = circlize::colorRamp2(c(0,0.05,0.1), c("lightblue","white","firebrick")), 
        right_annotation = rAnno, 
        cluster_rows = T, 
        cluster_columns = F)

        