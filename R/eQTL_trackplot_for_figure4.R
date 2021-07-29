library(ArchR)
library(Gviz)
library(GenomicInteractions)
library(ComplexHeatmap)
setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

satac <- suppressMessages(loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/'))
peak.markers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')
peak.set <- satac@peakSet
peak.set$peakID <- GRToString(peak.set)

# prepare eQTL co-accessibility
enhancer.coacc <- readRDS('ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coacc_ENHANCERS_AllCellTypes_overlapCutoff50_k200_corr_cut_-1_maxDist_1Mb_hg38.gr.rds')
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
eqtl.gene.coacc.tbl <- eqtl.gene.coacc %>% as_tibble() %>% group_by(snp, coacc_gene_name) %>% arrange(-correlation) %>% dplyr::slice(1)


# do the plotting

bw.files <- list.files(path = 'ArchR/ArchR_heart_latest_noAtrium/GroupBigWigs/CellTypes/', pattern = '*.bw$', full.names = T)
bw.files <- bw.files[c(1,2,3)]
peaks <- satac@peakSet

bw.list <- lapply(bw.files, function(x){rtracklayer::import(x)})
names(bw.list) <- sub('\\.', ' ', sub('-.*', '', basename(bw.files)))

atac.tracks <- lapply(names(bw.list), function(x){
    DataTrack(range = bw.list[[x]], type = 'h', col = palette[x], name = x, showAxis=T, ylim=c(0,0.06))
})

par.path <- '/project2/gca/Heart_Atlas/Nuc_seq/SP-HE-HE200915RNA-175RNA/output/'
celltypes.bams <- c("Cardiomyocyte.bam","Endothelial.bam","Fibroblast.bam")
celltypes <- sub('\\.bam','',celltypes.bams)
bam.files <- paste0(par.path, celltypes.bams)
names(bam.files) <- celltypes
#nreads <- c(11188127, 20603239, 16849499)
nreads <- c(16669819, 16374895, 16374895)
names(nreads) <- celltypes

rna.tracks <- lapply(names(bam.files), function(x){
    DataTrack(range = bam.files[[x]], type = 'h', col = palette[x], name = x, showAxis=T, ylim=c(0,5e-6),transformation = function(y){y/nreads[[x]]})
})

union.calls <- DataTrack(range = peaks, type = 'h', genome = "hg38", name = "Union", col = "grey", showAxis=FALSE, ylim = c(0,1))

interest.eqtl <- eqtl.snp.gr[eqtl.snp.gr$eGene=="SOX9",]
curr.locus.gr <- GRanges(seqnames = seqnames(interest.eqtl), 
                         ranges = IRanges(start = start(interest.eqtl)-2000, 
                                          end = end(interest.eqtl)+20000))

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
plotTracks(c(coacc.track, union.calls, atac.tracks, rna.tracks, gene.track),
           chromosome = as.character(seqnames(curr.locus.gr)), 
           transcriptAnnotation = "symbol", 
           collapseTranscripts= "longest", 
           from = start(curr.locus.gr), 
           to = end(curr.locus.gr), sizes = c(1, 0.2, rep(0.5, 6), 2), 
           panel.only = F)
dev.off()
