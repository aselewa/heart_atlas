library(ArchR)
library(tidyverse)
library(Gviz)
library(GenomicInteractions)
library(ComplexHeatmap)

setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')
source('R/track_plotter.R')
palette <- readRDS('notebooks/palette.rds')

satac <- suppressMessages(loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/'))
peak.markers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')
celltype_ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")
peak.set <- satac@peakSet
peak.set$peakID <- GRToString(peak.set)

    
# Run notebooks/eqtl_enrichment before running this script 
# 4A is from notebooks/eqtl_enrichment.rmd

# 4B,C eQTL overlap with our peaks

eqtl.snp.gr <- readRDS('eQTL_enrich/gtex_finemapping/HighPIP_eQTLs_w_LFSR.gr.rds')

rand.snps <- readLines('matched_SNPs/eQTL_top_pip_hg19_5batches/snpsnap_match_hg38.txt') # load random SNPs
rand.gr <- StringToGR(rand.snps)
rand.gr <- rand.gr[seqnames(rand.gr) %in% paste0("chr",1:22),]
seqlevels(rand.gr) <-  paste0("chr",1:22)
rand.gr$snp <- paste0(seqnames(rand.gr),'_',start(rand.gr))

peak.markers.union <- split(peak.set, peak.set$peakType)
eqtl.overlap.peaks <- lapply(peak.markers.union, function(x){plyranges::join_overlap_inner(eqtl.snp.gr, x)})
eqtl.overlap.peaks.prop <- sapply(eqtl.overlap.peaks, function(x){length(x)/length(eqtl.snp.gr)})

random.overlap.peaks <- lapply(peak.markers.union, function(x){plyranges::join_overlap_inner(rand.gr, x)})
random.overlap.peaks.prop <- sapply(random.overlap.peaks, function(x){length(x)/length(rand.gr)})

overlap.peaks.df <- data.frame(overlap.peaks.prop = c(eqtl.overlap.peaks.prop, random.overlap.peaks.prop))
overlap.peaks.df$celltypes <- factor(rep(names(peak.markers.union),2))
                                     #, levels = c("union",celltype_ideal_order))
overlap.peaks.df$snpType <- rep(c("eQTLs","Random"), each = length(peak.markers.union))

pdf('manuscript_figures/figure4/Fig4B_union_peaktype_overlap.pdf', width = 6, height=8)
ggplot(overlap.peaks.df, aes(x = celltypes, y = overlap.peaks.prop, fill = snpType)) + 
    geom_bar(stat='identity', position='dodge') + ggClean(rotate_axis = T) +
    scale_fill_manual(values=c("black","grey")) + 
    xlab('Cell Types') + ylab('% eQTLs')
dev.off()

peak.markers.union <- peak.markers
eqtl.overlap.peaks <- lapply(peak.markers.union, function(x){plyranges::join_overlap_inner(eqtl.snp.gr, x)})
eqtl.overlap.peaks.prop <- sapply(eqtl.overlap.peaks, function(x){length(x)/length(eqtl.snp.gr)})

random.overlap.peaks <- lapply(peak.markers.union, function(x){plyranges::join_overlap_inner(rand.gr, x)})
random.overlap.peaks.prop <- sapply(random.overlap.peaks, function(x){length(x)/length(rand.gr)})

overlap.peaks.df <- data.frame(overlap.peaks.prop = c(eqtl.overlap.peaks.prop, random.overlap.peaks.prop))
overlap.peaks.df$celltypes <- factor(rep(names(peak.markers.union),2), levels = c("union",celltype_ideal_order))
overlap.peaks.df$snpType <- rep(c("eQTLs","Random"), each = length(peak.markers.union))

pdf('manuscript_figures/figure4/Fig4B_celltype_overlap.pdf', width = 10, height=8)

ggplot(overlap.peaks.df, aes(x = celltypes, y = overlap.peaks.prop, fill = snpType)) + 
    geom_bar(stat='identity', position='dodge') + ggClean(rotate_axis = T) +
    scale_fill_manual(values=c("black","grey")) + 
    xlab('Cell Types') + ylab('% eQTLs')

dev.off()

# 4F eGene expression

srna <- readRDS('seurat/Heart_RNA_Processed_Combined_NoAtrium.rds')
srna.exp <- Seurat::AverageExpression(srna)
srna.exp.mat <- as.matrix(log2(srna.exp$RNA + 1))
srna.exp.mat.z <- sweep(srna.exp.mat - rowMeans(srna.exp.mat), MARGIN = 1, STATS = matrixStats::rowSds(srna.exp.mat), FUN = '/')

finemap.res <- readr::read_tsv('eQTL_enrich/gtex_finemapping/GTEx_v8_finemapping_DAPG/Heart_LV_Finemapping_CS95.txt', col_names = F)
colnames(finemap.res) <- c("tissue","gene","cluster_id","cluster_pip","variant_id","variant_pip")
finemap.res <- finemap.res[finemap.res$variant_id %in% eqtl.snp.gr$SNP,] %>% 
    dplyr::select(variant_id, gene) %>% rename(gene_id = gene, SNP = variant_id) %>% mutate(gene_id = sub('\\..*', '', gene_id))

ensembl.to.symbol <- AnnotationDbi::select(x = org.Hs.eg.db::org.Hs.eg.db, keys = finemap.res$gene_id, columns = "SYMBOL", keytype = "ENSEMBL") %>% 
    rename(gene_id = ENSEMBL) %>%
    group_by(gene_id) %>% slice(1)

finemap.res <- left_join(finemap.res, ensembl.to.symbol)

peak.markers.union <- peak.markers
eqtl.overlap.peaks <- lapply(peak.markers.union, function(x){plyranges::join_overlap_inner(eqtl.snp.gr, x)})
eqtl.overlap.peaks$union <- NULL
ns <- rep(names(eqtl.overlap.peaks), lengths(eqtl.overlap.peaks))
eqtl.overlap.peaks <- unlist(GRangesList(eqtl.overlap.peaks))
eqtl.overlap.peaks$cellType <- factor(ns, levels = celltype_ideal_order)
eqtl.overlap.peaks.df <- eqtl.overlap.peaks %>% as_tibble() %>% select(SNP, cellType)

finemap.res <- finemap.res %>% left_join(., eqtl.overlap.peaks.df, on='SNP')

finemap.res.inRNA <- finemap.res[finemap.res$SYMBOL %in% rownames(srna.exp.mat),] %>% arrange(cellType) %>% filter(!is.na(cellType))
srna.egene.mat <- srna.exp.mat.z[finemap.res.inRNA$SYMBOL, ]

mat.list <- list()
for(ct in celltype_ideal_order){
    curr.mat <- srna.egene.mat[finemap.res.inRNA$cellType == ct,]
    mat.list[[ct]] <- curr.mat[order(curr.mat[, ct], decreasing = T),]
}
srna.egene.mat.sorted <- Reduce(rbind, mat.list)
srna.egene.mat.sorted <- srna.egene.mat.sorted[,celltype_ideal_order]
rownames(srna.egene.mat.sorted) <- make.unique(rownames(srna.egene.mat.sorted))

genes.to.show <- c("ESRRB","MAST4","SOX9","FGD4","CAMK1D","COL18A1","EVA1B","DLK1","IL32.1","CDK11B","MICB.1","STX3")

pdf('manuscript_figures/figure4/Fig4_celltype_specific_eGene_exp.pdf', width=10, height=6)
p <- Heatmap(t(srna.egene.mat.sorted), show_row_names = T, show_column_names = F,
        cluster_rows = F, cluster_columns = F, 
        bottom_annotation  = HeatmapAnnotation(which = "column", 
                                            cluster = finemap.res.inRNA$cellType,
                                            col = list(cluster = palette[celltype_ideal_order])),
        top_annotation = HeatmapAnnotation(which = "column", 
                                            genestoshow = anno_mark(at = match(genes.to.show, rownames(srna.egene.mat.sorted)), 
                                                                    labels = sub('\\..*','',genes.to.show), side = 'top')),
        col = circlize::colorRamp2(c(-3, 0, 3), c("lightblue","white","firebrick")),
        heatmap_legend_param = list(direction = "horizontal"))
draw(p, heatmap_legend_side = "bottom")
dev.off()

# shared vs specific eGene expression
srna <- readRDS('seurat/Heart_RNA_Processed_Combined_NoAtrium.rds')
srna.exp <- Seurat::AverageExpression(srna)
srna.exp.mat <- as.matrix(log2(srna.exp$RNA + 1))
srna.exp.mat.z <- sweep(srna.exp.mat - rowMeans(srna.exp.mat), MARGIN = 1, STATS = matrixStats::rowSds(srna.exp.mat), FUN = '/')

finemap.res <- readr::read_tsv('eQTL_enrich/gtex_finemapping/GTEx_v8_finemapping_DAPG/Heart_LV_Finemapping_CS95.txt', col_names = F)
colnames(finemap.res) <- c("tissue","gene","cluster_id","cluster_pip","variant_id","variant_pip")
finemap.res <- finemap.res[finemap.res$variant_id %in% eqtl.snp.gr$SNP,] %>% 
    dplyr::select(variant_id, gene) %>% rename(gene_id = gene, SNP = variant_id) %>% mutate(gene_id = sub('\\..*', '', gene_id))

ensembl.to.symbol <- AnnotationDbi::select(x = org.Hs.eg.db::org.Hs.eg.db, keys = finemap.res$gene_id, columns = "SYMBOL", keytype = "ENSEMBL") %>% 
    rename(gene_id = ENSEMBL) %>%
    group_by(gene_id) %>% slice(1)

finemap.res <- left_join(finemap.res, ensembl.to.symbol)

peak.markers.union <- peak.markers
peak.markers.union$Shared <- peak.set[!(peak.set$peakID %in% GRToString(unlist(peak.markers))),]
eqtl.overlap.peaks <- lapply(peak.markers.union, function(x){plyranges::join_overlap_inner(eqtl.snp.gr, x)})

ns <- rep(names(eqtl.overlap.peaks), lengths(eqtl.overlap.peaks))
eqtl.overlap.peaks <- unlist(GRangesList(eqtl.overlap.peaks))
eqtl.overlap.peaks$cellType <- factor(ns, levels = c(celltype_ideal_order,"Shared"))
eqtl.overlap.peaks.df <- eqtl.overlap.peaks %>% as_tibble() %>% dplyr::select(SNP, cellType)

finemap.res <- finemap.res %>% left_join(., eqtl.overlap.peaks.df, on='SNP')

finemap.res.inRNA <- finemap.res[finemap.res$SYMBOL %in% rownames(srna.exp.mat),] %>% arrange(cellType) %>% filter(!is.na(cellType))
srna.egene.mat <- srna.exp.mat.z[finemap.res.inRNA$SYMBOL, ]

mat.list <- list()
for(ct in celltype_ideal_order){
    curr.mat <- srna.egene.mat[finemap.res.inRNA$cellType == ct,]
    mat.list[[ct]] <- curr.mat[order(curr.mat[, ct], decreasing = T),]
}
mat.list[["Shared"]] <- srna.egene.mat[finemap.res.inRNA$cellType=="Shared",]
srna.egene.mat.sorted <- Reduce(rbind, mat.list)
srna.egene.mat.sorted <- srna.egene.mat.sorted[,celltype_ideal_order]
rownames(srna.egene.mat.sorted) <- make.unique(rownames(srna.egene.mat.sorted))

Heatmap(srna.egene.mat.sorted, 
        cluster_rows = T, 
        cluster_columns = F, 
        show_row_names = F, 
        left_annotation = rowAnnotation(type = finemap.res.inRNA$cellType),
        col = circlize::colorRamp2(c(-3, 0, 3), c("lightblue","white","firebrick")))


# prepare eQTL co-accessibility
enhancer.coacc <- readRDS('ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coaccess_enhancers_corr_0.5_maxDist_1Mb_hg38.gr.rds')
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
eqtl.gene.coacc <- eqtl.gene.coacc[!is.na(eqtl.gene.coacc$coacc_gene_name),]
eqtl.gene.coacc <- eqtl.gene.coacc[eqtl.gene.coacc$eGene == eqtl.gene.coacc$coacc_gene_name,]
eqtl.gene.coacc.tbl <- eqtl.gene.coacc %>% as_tibble() %>% group_by(snp, coacc_gene_name) %>% arrange(-correlation) %>% slice(1)


# 4E

bw.files <- list.files(path = 'ArchR/ArchR_heart_latest_noAtrium/GroupBigWigs/CellTypes/', pattern = '*.bw$', full.names = T)
peaks <- satac@peakSet

bw.list <- lapply(bw.files, function(x){rtracklayer::import(x)})
names(bw.list) <- sub('\\.', ' ', sub('-.*', '', basename(bw.files)))
bw.list <- lapply(celltype_ideal_order, function(x){ bw.list[[x]] })
names(bw.list) <- celltype_ideal_order

atac.tracks <- lapply(names(bw.list), function(x){
    DataTrack(range = bw.list[[x]], type = 'h', col = palette[x], name = x, showAxis=T, ylim=c(0,0.06))
})

union.calls <- DataTrack(range = peaks, type = 'h', genome = "hg38", name = "Union", col = "grey", showAxis=FALSE, ylim = c(0,1))

interest.eqtl <- eqtl.snp.gr[eqtl.snp.gr$eGene=="ESRRB",]
curr.locus.gr <- GRanges(seqnames = seqnames(interest.eqtl), ranges = IRanges(start = start(interest.eqtl)-10000, end = end(interest.eqtl)+20000))

promoter.gr <- StringToGR(eqtl.gene.coacc$promoterID)
snp.promoter <- GenomicInteractions::GenomicInteractions(anchor1 = eqtl.gene.coacc, anchor2 = promoter.gr)
snp.promoter$counts <- round(100*eqtl.gene.coacc$correlation)
snp.promoter <- snp.promoter[snp.promoter$anchor1.coacc_gene_name == "ESRRB",]

coacc.track <- InteractionTrack(snp.promoter, name = "eQTL_Coaccess")
dpars <- list(col.interactions="red", 
              col.anchors.fill ="blue", 
              col.anchors.line = "black", 
              interaction.dimension = "height",  
              interaction.measure = "counts", 
              plot.trans = FALSE, 
              plot.outside = FALSE,  
              col.outside="lightblue",  
              anchor.height = 0.1)
displayPars(coacc.track) <- dpars

gene.track <- knownGeneObject(curr.locus.gr = curr.locus.gr, genome = "hg38")

pdf('manuscript_figures/figure4/ESRRB_eQTL_trackplot.pdf', width=6, height=5)
plotTracks(c(coacc.track, union.calls, atac.tracks, gene.track),
           chromosome = as.character(seqnames(curr.locus.gr)), 
           transcriptAnnotation = "symbol", 
           collapseTranscripts= 'longest', 
           from = start(curr.locus.gr), 
           to = end(curr.locus.gr), sizes = c(1, 0.2, rep(0.5, 8), 1), 
           panel.only = T)
dev.off()


