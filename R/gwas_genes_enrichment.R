library(tidyverse)

setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

srna <- readRDS('seurat/Heart_RNA_Processed_Combined_NoAtrium.rds')
srna.exp <- Seurat::AverageExpression(srna)
avg.exp.mat <- log2(srna.exp$RNA + 1)

gene.annots <- readRDS('genomic_annotations/hg19_gtf_genomic_annots.gr.rds')
genes <- gene.annots$genes

diff.exp.res <- readRDS('seurat/diff_expr_markers.df.rds')
finemap.res <- readRDS('GWAS/finemapping/aFib_Finemapped_GeneMapped_06152021.tble.rds')

high.pip.snps <- finemap.res %>% distinct(snp, .keep_all=TRUE) %>% group_by(locus) %>% arrange(-pip) %>% slice(1)
high.pip.snp.gr <- GRanges(seqnames = paste0("chr",high.pip.snps$chr), ranges = IRanges(start = high.pip.snps$pos, end = high.pip.snps$pos))

top.pip.genes <- unique(finemap.res$gene_name[finemap.res$gene_pip >= 0.5])
hits <- findOverlaps(high.pip.snp.gr, genes, maxgap = 250000)
control.genes <- unique(genes$gene_name[subjectHits(hits)])

cell.types <- unique(diff.exp.res$cluster)
top.pip.genes.overlap <- sapply(cell.types, function(x){
    length(intersect(diff.exp.res$gene[diff.exp.res$cluster==x], top.pip.genes))/length(top.pip.genes)})

control.genes.overlap <- sapply(cell.types, function(x){
    length(intersect(diff.exp.res$gene[diff.exp.res$cluster==x], control.genes))/length(control.genes)})

compare.overlap.df <- data.frame(prop = c(top.pip.genes.overlap, control.genes.overlap), 
                                 celltype = as.character(cell.types),
                                 gene_type = rep(c("PIP>0.5","Control"),each=length(cell.types)))

compare.overlap.df <- compare.overlap.df[compare.overlap.df$celltype %in% c("Cardiomyocyte","Endothelial","Fibroblast"),]

pdf('manuscript_figures/GWAS_Genes_Control.pdf',width=8, height=6)
ggplot(compare.overlap.df, aes(x=celltype, y=prop, fill=gene_type)) + geom_bar(stat='identity', position='dodge', width = 0.7) + ggClean(rotate_axis = T) +
    xlab('') + ylab('Prop. Genes Differentially Expressed') + scale_fill_brewer(palette = "Paired")
dev.off()


avg.exp.mat.top.pip <- avg.exp.mat[intersect(rownames(avg.exp.mat), top.pip.genes),]
avg.exp.mat.top.pip <- pivot_longer(data = avg.exp.mat.top.pip, cols = everything())
avg.exp.mat.top.pip$gene_type <- "PIP>0.5"

avg.exp.mat.control <- avg.exp.mat[intersect(rownames(avg.exp.mat), control.genes),]
avg.exp.mat.control <- pivot_longer(data = avg.exp.mat.control, cols = everything())
avg.exp.mat.control$gene_type <- "Control"

exp.compare.df <- bind_rows(avg.exp.mat.top.pip,avg.exp.mat.control)
exp.compare.df$value[exp.compare.df$value > 3] <- 3

ggplot(exp.compare.df, aes(x=name, y=value, fill=gene_type)) + geom_boxplot(outlier.size = 0.5) + ggClean(rotate_axis = T) + 
    xlab('') + ylab('log2 TP10k')


