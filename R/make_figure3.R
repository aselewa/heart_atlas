library(ArchR)
library(tidyverse)
require(plyranges)
source('R/analysis_utils.R')
palette <- readRDS('notebooks/palette.rds')
setwd('/project2/gca/aselewa/heart_atlas_project/')

# load a bunch of stuff
satac <- suppressMessages(loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/'))

peaks <- getPeakSet(satac)
peaks$peakID <- GRToString(peaks)
peaks <- peaks[peaks$peakType!="Promoter",]
peak.markers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')

genomic.annots <- readRDS('hg38_gtf_genomic_annots.gr.rds')
gene.annots <- genomic.annots$genes

#annotate all peaks with DA test results
celltype_ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")
peak.markers <- lapply(celltype_ideal_order, function(x){peak.markers[[x]]})
names(peak.markers) <- celltype_ideal_order
peak.markers.str <- unlist(lapply(peak.markers, function(x){GRToString(x)}), use.names = F)
peaks.celltypes <- data.frame(peakID=peak.markers.str, cellType = factor(rep(names(peak.markers), lengths(peak.markers)), levels = names(peak.markers)))

peak.info.df <- peaks %>% as_tibble() %>% left_join(., peaks.celltypes, on = "peakID") %>% select(peakID, peakType, nearestGene, distToGeneStart, cellType)
levels(peak.info.df$cellType) <- c(levels(peak.info.df$cellType), "Shared")
peak.info.df$cellType[is.na(peak.info.df$cellType)] <- "Shared"

# accessibility with gene promoters
enhancer.coacc <- readRDS('ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coaccess_enhancers_corr_0.5_maxDist_1Mb_hg38.gr.rds')

enhancer.coacc.tbl <- enhancer.coacc %>% 
    as_tibble() %>%  
    select(peakID, coacc_gene_name, correlation, distToCoaccGene) %>%
    group_by(peakID) %>% arrange(-correlation) %>% slice(1) %>% ungroup()

peak.info.coacc.df <- inner_join(peak.info.df, enhancer.coacc.tbl, on = "peakID") %>% 
    mutate(coacc_gene_name = ifelse(is.na(coacc_gene_name), nearestGene, coacc_gene_name)) %>%
    mutate(correlation = ifelse(is.na(correlation), 0, correlation)) %>%
    mutate(distToCoaccGene = ifelse(is.na(distToCoaccGene), distToGeneStart, distToCoaccGene)) %>%
    distinct(peakID, coacc_gene_name, .keep_all = T)

npeaks.per.gene <- peak.info.coacc.df %>% 
    group_by(coacc_gene_name) %>%
    summarise(peakCount = dplyr::n()) %>% 
    mutate(peakFreq = peakCount/sum(peakCount))

bks <- c(0, 2, 10, 20, 30, 40, 50, 60, 10000)
labs <- c("1-2","3-10","11-20", "21-30","31-40", "41-50","51-60","61+")
peakCount.bins <- npeaks.per.gene %>% summarise(nGenes = as.vector(table(cut(peakCount, breaks = bks, labels = labs)))) %>% mutate(nPeaks = labs)
peakCount.bins$nPeaks <- factor(peakCount.bins$nPeaks, levels = labs)

# 3A is cartoon of co-accessibility

# 3B with inset, number of peaks per gene
pdf('manuscript_figures/figure3/Fig3_percent_peaks_distribution.pdf', width = 8, height=6)
ggplot(peakCount.bins, aes(x = nPeaks, y = nGenes/1e3)) + 
    geom_bar(stat='identity', position='dodge', fill="#238b45", width=0.7) + ggClean() +
    ylab('Number of Genes (x 1000)') + xlab('Number of Peaks') + scale_y_continuous(breaks = seq(0, 8, 1))
dev.off()

dist.df <- data.frame(dist = log10(c(peak.info.coacc.df$distToCoaccGene, peak.info.df$distToGeneStart)),
                      type = c(rep("coacc", nrow(peak.info.coacc.df)), rep("nearest", nrow(peak.info.df))))

pdf('manuscript_figures/figure3/Fig3A_2_distance_dist.pdf', width = 8, height=6)
ggplot(dist.df, aes(x = dist, color=type)) + geom_density(adjust=2) + ggClean() + xlim(c(0,7))
dev.off()

ggplot(npeaks.per.gene.celltype, aes(x = nPeaks, y = nGenes/1e3, fill = cellType)) + 
    geom_bar(stat='identity', position='dodge', width=0.7) + ggClean() +
    ylab('Number of Genes (x 1000)') + xlab('Number of Peaks') + scale_y_continuous(breaks = seq(0, 8, 1))

# break down by cell type

set.seed(100)
npeaks.per.gene.celltype <- peak.info.coacc.df %>% 
    group_by(cellType) %>% sample_n(size = 900) %>% ungroup() %>%
    group_by(cellType, coacc_gene_name) %>%
    summarise(peakCount = dplyr::n()) %>% 
    mutate(peakFreq = peakCount/sum(peakCount))

quants <- seq(1, 15, 1)
npeaks.per.celltype.cdf <- npeaks.per.gene.celltype %>% 
    group_by(cellType) %>% 
    arrange(peakFreq) %>% 
    summarise(ecdf = sapply(quants, function(x){mean(peakCount < x)}), quants=quants)

ggplot(npeaks.per.celltype.cdf, aes(x = quants, y = ecdf, color=cellType)) + geom_line() + ggClean() + xlab('Number of Peaks per Gene') + ylab('CDF') +
    scale_x_continuous(breaks = seq(0, 16, 4))

# # 3A 2 correlation of # of peaks with transcriptional units
# 
# gene.annots.tbl <- gene.annots %>% as_tibble() %>% select(seqnames, start, end, symbol) %>% mutate(gene_length = end-start) %>% dplyr::rename(coacc_gene_name = symbol)
# 
# bks <- c(0, 2, 20, 10000)
# labs <- c("<=2", "3-20","20+")
# npeaks.per.gene.length <- left_join(npeaks.per.gene, gene.annots.tbl, on="coacc_gene_name") %>% 
#     mutate(nGenes = cut(peakCount, breaks = bks, labels = labs))
# 
# ggplot(npeaks.per.gene.length, aes(x =  nGenes, y = log10(gene_length), color=nGenes)) + geom_violin(color = 'black') + geom_boxplot(width=0.1)  + ggClean()

# 3B is GO enrichment, run script R/GO_enrichment.R
    
# 3C, overlap with RNA genes
rna.markers <- readRDS('seurat/diff_expr_markers.df.rds')
rna.markers$cluster <- as.character(rna.markers$cluster)

celltypes <- celltype_ideal_order
gene.overlap <- matrix(0, nrow=length(celltypes), ncol=length(celltypes))
for(i in 1:length(celltypes)){
    for(j in 1:length(celltypes)){
        gene.overlap[i,j] <- sum(rna.markers$gene[rna.markers$cluster==celltypes[i]] %in% peak.info.coacc.df$coacc_gene_name[peak.info.coacc.df$cellType==celltypes[j]])/sum(rna.markers$cluster==celltypes[i])
    }
}

rownames(gene.overlap) <- celltypes
colnames(gene.overlap) <- paste0("PeakGenes_", celltypes)

pdf('manuscript_figures/figure3/Fig3B_overlap_w_rna.pdf', width = 8, height=6)

ComplexHeatmap::Heatmap(gene.overlap, cluster_rows = F, cluster_columns = F, col = circlize::colorRamp2(c(0, 0.6, 1), c("lightblue","white","firebrick")))

dev.off()

# 3D eQTL overlap with our peaks

finemap.res <- suppressMessages(readr::read_tsv('eQTL_enrich/gtex_finemapping/GTEx_v8_finemapping_DAPG/Heart_LV_Finemapping_CS95.txt', col_names = F))
colnames(finemap.res) <- c("tissue","gene","cluster_id","cluster_pip","variant_id","variant_pip")

finemap.res <- finemap.res[finemap.res$variant_pip > 0.8,]
finemap.res$genes_clean <- sub('[.].*', '', finemap.res$gene)
gene_to_symbol <- ensembldb::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, keys= finemap.res$genes_clean, keytype = "GENEID", columns = "SYMBOL")
colnames(gene_to_symbol) <- c("genes_clean","symbol")
finemap.res <- left_join(finemap.res, gene_to_symbol, on = 'genes_clean')
finemap.res <- finemap.res[!is.na(finemap.res$symbol),]

eqtl.snp.gr <- snpIDtoGR(unique(finemap.res$variant_id))

peak.markers.union <- peak.markers
peak.markers.union[['union']] <- peaks
eqtl.overlap.peaks <- lapply(peak.markers.union, function(x){plyranges::join_overlap_inner(eqtl.snp.gr, x)})
eqtl.overlap.peaks.prop <- sapply(eqtl.overlap.peaks, function(x){100*length(x)/length(eqtl.snp.gr)})

eqtl.overlap.peaks.df <- data.frame(eqtl.overlap.peaks.prop)
eqtl.overlap.peaks.df$celltypes <- factor(rownames(eqtl.overlap.peaks.df), levels = c(rev(celltype_ideal_order), "union"))

pdf('manuscript_figures/figure3/Fig3D_eqtl_overlap.pdf', width = 10, height=8)

ggplot(eqtl.overlap.peaks.df, aes(y = celltypes, x = eqtl.overlap.peaks.prop, fill = celltypes)) + 
    geom_bar(stat='identity', position='dodge') + ggClean() +
    scale_fill_manual(values = c(palette[rev(celltype_ideal_order)],"grey")) +
    ylab('Cell Types') + xlab('% eQTLs') + coord_cartesian(xlim=c(0, 15))

dev.off()

# 3E eGene expression

srna <- readRDS('seurat/Heart_RNA_Processed_Combined_NoAtrium.rds')
srna.exp <- Seurat::AverageExpression(srna)
srna.exp.mat <- as.matrix(log2(srna.exp$RNA + 1))
rSDs <- matrixStats::rowSds(srna.exp.mat)
srna.exp.mat <- srna.exp.mat[rSDs > 0,]
srna.exp.mat.z <- sweep(srna.exp.mat - rowMeans(srna.exp.mat), MARGIN = 1, STATS = matrixStats::rowSds(srna.exp.mat), FUN = '/')

egenes <- intersect(rownames(srna.exp.mat.z), finemap.res$symbol)
srna.egene.mat <- srna.exp.mat.z[egenes, ]

set.seed(100)
kres <- kmeans(srna.egene.mat, centers = 10)
labs <- kres$cluster
labs.sort <- sort(labs)

Heatmap(srna.exp.mat.z[names(labs.sort), celltype_ideal_order], show_row_names = F, cluster_rows = F, cluster_columns = F, left_annotation = rowAnnotation(cluster = factor(labs.sort)))

lab_map <- data.frame(num = 1:10,
                      celltype = c("mix1","mix2","Neuronal", "Myeloid", "Lymphoid", "mix3", "Endothelial", "Fibroblast", "Cardiomyocyte","mix4"),
                      stringsAsFactors = F)
ideal.order <- c("Cardiomyocyte","mix1","Endothelial","Fibroblast","Neuronal","mix3","mix2","mix4","Lymphoid","Myeloid")

labs.named <- plyr::mapvalues(x = labs, from = lab_map$num, to = lab_map$celltype)
labs.named <- factor(labs.named, levels = ideal.order)
labs.sort <- sort(labs.named)
srna.exp.mat.z <- srna.exp.mat.z[names(labs.sort), celltype_ideal_order]

pdf('manuscript_figures/figure3/Fig3F_eGeneExpression.pdf', width = 10, height=11)

Heatmap(srna.exp.mat.z, show_row_names = F, 
        cluster_rows = F, 
        cluster_columns = F, 
        col = circlize::colorRamp2(c(-4, 0, 4), c("lightblue","white","firebrick")), name = 'Expression Z-score',
        use_raster = T)


dev.off()

# 3F, eQTL eGenes co-accessibility

finemap.res <- suppressMessages(readr::read_tsv('eQTL_enrich/gtex_finemapping/GTEx_v8_finemapping_DAPG/Heart_LV_Finemapping_CS95.txt', col_names = F))
colnames(finemap.res) <- c("tissue","gene","cluster_id","cluster_pip","variant_id","variant_pip")

finemap.res <- finemap.res[finemap.res$variant_pip > 0.8,]
finemap.res$genes_clean <- sub('[.].*', '', finemap.res$gene)
gene_to_symbol <- ensembldb::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, keys= finemap.res$genes_clean, keytype = "GENEID", columns = "SYMBOL")
colnames(gene_to_symbol) <- c("genes_clean","symbol")
finemap.res <- left_join(finemap.res, gene_to_symbol, on = 'genes_clean')
finemap.res <- finemap.res[!is.na(finemap.res$symbol),] # some genes cant map to symbols, e.g. lncRNAs

eqtl.snp.gr <- snpIDtoGR(finemap.res$variant_id)
eqtl.snp.gr$eGene <- finemap.res$symbol
eqtl.snp.gr$snp <- finemap.res$variant_id

eqtl.gene.coacc <- plyranges::join_overlap_left(eqtl.snp.gr, enhancer.coacc)
eqtl.gene.coacc <- eqtl.gene.coacc[!is.na(eqtl.gene.coacc$gene_name),]
eqtl.gene.coacc <- eqtl.gene.coacc[eqtl.gene.coacc$eGene == eqtl.gene.coacc$gene_name,]
eqtl.gene.coacc.tbl <- eqtl.gene.coacc %>% as_tibble()

rand.snps <- readLines('matched_SNPs/eQTL_top_pip_hg19_5batches/snpsnap_match_hg38.txt') # load random SNPs
rand.gr <- StringToGR(rand.snps)
rand.gr <- rand.gr[seqnames(rand.gr) %in% paste0("chr",1:22),]
seqlevels(rand.gr) <-  paste0("chr",1:22)
rand.gr$snp <- paste0(seqnames(rand.gr),'_',start(rand.gr))

rand.gene.coacc <- plyranges::join_overlap_left(rand.gr, enhancer.coacc)
rand.gene.coacc.tbl <- rand.gene.coacc %>% as_tibble() %>% filter(!is.na(correlation))

plot.df <- data.frame(correlation = c(eqtl.gene.coacc.tbl$correlation, rand.gene.coacc.tbl$correlation, enhancer.coacc$correlation), 
                      type= c(rep("A_eQTLs",nrow(eqtl.gene.coacc.tbl)), rep("B_random", nrow(rand.gene.coacc.tbl)), rep("C_all", length(enhancer.coacc))))

pdf('manuscript_figures/figure3/Fig4D_coaccess_correlation.pdf', width = 7, height=8)

ggplot(plot.df, aes(x=type, y = correlation, fill=type)) + 
    geom_boxplot() + 
    ggClean() + 
    scale_fill_manual(values = c("#238b45", "#bae4b3", "#aaaaaa"))

dev.off()


# 3G track plot

interest.eqtl <- eqtl.snp.gr[eqtl.snp.gr$eGene=="ESRRB",]
curr.locus.gr <- GRanges(seqnames = seqnames(interest.eqtl), ranges = IRanges(start = start(interest.eqtl)-14000, end = end(interest.eqtl)+14000))

bw.files <- c("ArchR/ArchR_heart_latest_noAtrium/GroupBigWigs/CellTypes/Cardiomyocyte-TileSize-100-normMethod-ReadsInTSS-ArchR.bw",
              "ArchR/ArchR_heart_latest_noAtrium/GroupBigWigs/CellTypes/Endothelial-TileSize-100-normMethod-ReadsInTSS-ArchR.bw",
              "ArchR/ArchR_heart_latest_noAtrium/GroupBigWigs/CellTypes/Neuronal-TileSize-100-normMethod-ReadsInTSS-ArchR.bw",
              "ArchR/ArchR_heart_latest_noAtrium/GroupBigWigs/CellTypes/Lymphoid-TileSize-100-normMethod-ReadsInTSS-ArchR.bw")

bw.list <- lapply(bw.files, function(x){rtracklayer::import(x)})


plot.list <- list()
plot.list[["gene.track"]] <- geneTrackPlot(curr.locus.gr = curr.locus.gr, collapseTranscripts = "longest", genome = "hg38")

plot.list[["cm.track"]] <- plotATAC(gr = bw.list[[1]], region.focus = curr.locus.gr) 
plot.list[["endo.track"]] <- plotATAC(gr = bw.list[[2]], region.focus = curr.locus.gr) 
plot.list[["neuro.track"]] <- plotATAC(gr = bw.list[[3]], region.focus = curr.locus.gr)
plot.list[["lymph.track"]] <- plotATAC(gr = bw.list[[4]], region.focus = curr.locus.gr)

hicplot.df <- enhancer.coacc %>% as_tibble() %>% rename(chr = seqnames) %>% 
    filter(chr == as.character(seqnames(interest.eqtl))) %>% 
    mutate(dist.to.eqtl = abs(start - start(interest.eqtl))) %>%
    mutate(isClose = dist.to.eqtl == min(dist.to.eqtl))

plot.list[["HiCplot"]] <- HiC.track(HiC.Midpoints = hicplot.df, curr.locus.gr = curr.locus.gr, links.to.highlight = hicplot.df$isClose)

pdf('manuscript_figures/figure3/Fig3G_ESRRB_coaccess.pdf', width = 20, height=6)
cowplot::plot_grid(plotlist = plot.list, align = 'v', ncol = 1, axis = 'b')
dev.off()

gene.exp <- data.frame(esrrb = srna.exp.mat["ESRRB",])
gene.exp$cellType <- factor(rownames(gene.exp), levels = celltype_ideal_order)

pdf('manuscript_figures/figure3/Fig3G_ESRRB_exp.pdf', width = 8, height=6)
ggplot(gene.exp, aes(x=cellType, y = esrrb, fill = cellType)) + 
    geom_bar(stat='identity', width=0.7) + 
    ggClean() + 
    scale_fill_manual(values = palette[celltype_ideal_order]) +
    xlab('') + ylab('log2 Expression')
dev.off()
