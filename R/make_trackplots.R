### Track plotting 
library(GenomicRanges) # required for all Gviz input
library(tidyverse) # manipulating tibbles
library(rtracklayer) # loading bigwigs/bed files
library(bigsnpr) # loading genotype data from 1000Genomes for LD calculation
library(Gviz) # main plotting
library(GenomicInteractions) # hic plots

setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')
celltype_ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")
palette <- readRDS('notebooks/palette.rds')

# genotype data and annotations hg19
big.snp <- bigsnpr::snp_attach('/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds')
gene.annots <- readRDS('genomic_annotations/hg19_gtf_genomic_annots.gr.rds')

# finemapping results
prior_res <- readRDS('GWAS/finemapping/aFib_Finemapped.tble.rds')
finemap.gr <- GenomicRanges::GRanges(seqnames = prior_res$chr, 
                                     ranges = IRanges::IRanges(start = prior_res$pos, end = prior_res$pos),
                                     snp = prior_res$snp,
                                     susie_pip = prior_res$susie_pip,
                                     chr = prior_res$chr,
                                     pos = prior_res$pos,
                                     locus = prior_res$locus)                                                                                                                                                                                                                       
seqlevelsStyle(finemap.gr) <- "UCSC"
finemap.gr <- finemap.gr[finemap.gr$pip > 1e-5, ]

# cell-type ATAC-seq bigWigs
bw.files <- list.files(path = 'ArchR/ArchR_heart_latest_noAtrium/GroupBigWigs/CellTypes/', pattern = 'Hg19.*', full.names = T)
bw.files <- bw.files[c(1,2,3)]

bw.list <- lapply(bw.files, function(x){rtracklayer::import(x)})
names(bw.list) <- sub('Hg19_', '' , sub('\\.', ' ', sub('-.*', '', basename(bw.files))))

# make an ATAC track for each cell type
atac.tracks <- lapply(names(bw.list), function(x){
    DataTrack(range = bw.list[[x]], type = 'h', col = palette[x], name = x, showAxis=F, ylim = c(0,0.8))
})

# rna tracks
par.path <- '/project2/gca/Heart_Atlas/Nuc_seq/SP-HE-HE200915RNA-397RNA/output/'
celltypes.bams <- c("Cardiomyocyte.bam","Endothelial.bam","Fibroblast.bam")
celltypes <- sub('\\.bam','',celltypes.bams)
bam.files <- paste0(par.path, celltypes.bams)
names(bam.files) <- celltypes
nreads <- c(11188127, 20603239, 16849499)
names(nreads) <- celltypes

rna.tracks <- lapply(names(bam.files), function(x){
    DataTrack(range = bam.files[[x]], 
              type = 'h',
              col = palette[x], 
              name = x, 
              showAxis=T, 
              transform = function(y){y/nreads[[x]]})
})


# h3k27ac tracl
h3k27 <- rtracklayer::import('ENCODE/H3k27ac_gwas_hg19/hg19_mapped/H3K27ac_heart_concat.bed')
seqlevelsStyle(h3k27) <- "UCSC"
seqlevels(h3k27, pruning.mode = "coarse") <- paste0("chr",1:22)
h3k27$score <- 1

fetal <- rtracklayer::import('ENCODE/H3k27ac_gwas_hg19/FetalHeart_E083-DNase_hg19_cleaned_narrowPeak.bed.gz')
seqlevelsStyle(fetal) <- "UCSC"
seqlevels(fetal, pruning.mode = "coarse") <- paste0("chr",1:22)
fetal$score <- 1

# make bed region tracks with Gviz
h3k27.track <- DataTrack(range = h3k27, type = "h", genome = "hg19", name = "H3K27ac", col = "navy", showAxis=FALSE, ylim = c(0,1))
fetal.track <- DataTrack(range = fetal, type="h", geome="hg19", name = "FetalDHS", col="gray", showAxis=FALSE, ylim=c(0,1))

# prepare HiC/coacc track data
pcHic <- readRDS('HiC/iPSC_CM_pcHiC_protein_Hg19.gr.rds') %>% as_tibble()

enhancer.pcHiC.gr <- GRanges(seqnames = pcHic$seqnames, ranges = IRanges(start = pcHic$start, end = pcHic$end), score = pcHic$score, gene = pcHic$gene_name)
promoter.pcHiC.gr <- GRanges(seqnames = pcHic$promoter_chr, ranges = IRanges(start = pcHic$promoter_start, end =pcHic$promoter_end), score = pcHic$score)
pchic.obj <- GenomicInteractions::GenomicInteractions(anchor1 = enhancer.pcHiC.gr, anchor2 = promoter.pcHiC.gr)
pchic.obj$counts <- round(pchic.obj$anchor1.score)

enhancer.coaccess <- readRDS('ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coacc_ENHANCERS_AllCellTypes_overlapCutoff50_k200_corr_cut_-1_maxDist_1Mb_hg19.gr.rds') %>% as_tibble() %>% 
    dplyr::select(enhancer_chr, enhancer_start, enhancer_end, correlation, idx, gene_name) 
promoter.coacc <- readRDS('ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coacc_PROMOTERS_AllCellTypes_overlapCutoff50_k200_corr_cut_-1_maxDist_1Mb_hg19.gr.rds') %>% as_tibble() %>% 
    dplyr::select(promoter_chr, promoter_start, promoter_end, idx)
coacc.df <- enhancer.coaccess %>% inner_join(., promoter.coacc, on = 'idx')

enhancer.coacc.df <- GRanges(seqnames = coacc.df$enhancer_chr, ranges = IRanges(start = coacc.df$enhancer_start, end = coacc.df$enhancer_end), corr = coacc.df$correlation, gene = coacc.df$gene_name)
promoter.coacc.df <- GRanges(seqnames = coacc.df$promoter_chr, ranges = IRanges(start = coacc.df$promoter_start, end = coacc.df$promoter_end))

coacc.obj <- GenomicInteractions::GenomicInteractions(anchor1 = enhancer.coacc.df, anchor2 = promoter.coacc.df)
coacc.obj <- coacc.obj[coacc.obj$anchor1.corr > 0.5,]
coacc.obj$counts <- round(100*(coacc.obj$anchor1.corr))

# actual plotting
final.mat <- readRDS('GWAS/finemapping/aFib_Finemapped_GeneMapped_ActivePromoter_07242021.gr.rds')
high.conf.snp.df <- final.mat %>% dplyr::filter(pip > 0.2) %>% group_by(snp) %>% arrange(-gene_pip) %>% dplyr::slice(1) 
gene.gr <- gene.annots$genes[match(high.conf.snp.df$gene_name, gene.annots$genes$gene_name),]
gene.gr.tss <- GRanges(seqnames = seqnames(gene.gr), ranges = IRanges(start = ifelse(strand(gene.gr)=="+", start(gene.gr), end(gene.gr))), symbol = gene.gr$gene_name)
high.conf.snp.df$gene.start.tss <- start(gene.gr.tss)

#genes.of.interest <- c("TBX5","ETV1","PLN","PITX2","HCN4","NKX2-5","PRRX1","TTN","SCN10A","FGF5","HAND2","CAMK2B","NFKB2","GATA4")
genes.of.interest <- 'CAMK2D'
locus.snp.dist <- high.conf.snp.df %>% 
    filter(gene_name %in% genes.of.interest) %>% 
    group_by(locus) %>% 
    arrange(-pip) %>% 
    dplyr::slice(1) %>% 
    mutate(distToTSS = pos-gene.start.tss) %>% 
    dplyr::select(locus, pos, distToTSS, gene_name)

for(l in 1:nrow(locus.snp.dist)){
    
    LOCUS = locus.snp.dist$locus[l]
    distToTSS = locus.snp.dist$distToTSS[l] 
    pip.df <- prior_res[prior_res$locus == LOCUS,] 
    snp.p <- locus.snp.dist$pos[l]
    
    ext <- 24000
    if(distToTSS < 0){ #snp is upstream
        ss <- snp.p - ext
        ee <- snp.p + abs(distToTSS) + ext
    } else{ # snp is downstream
        ss <- snp.p - abs(distToTSS) - ext
        ee <- snp.p + ext
    }
    
    # This one-item GRanges determines the boundary of the region we will visualize
    curr.locus.gr <- GRanges(seqnames = paste0("chr",pip.df$chr[1]), IRanges(start =  ss, 
                                                                             end = ee ))
    # prepare P-value track
    ## Add LD information
    pip.df$isAnnotate <- pip.df$susie_pip > 0.1
    top.snp <- pip.df$bigSNP_index[which.max(pip.df$pval)]
    top.snp.G <- big.snp$genotypes[,top.snp]
    G.mat <- big.snp$genotypes[,pip.df$bigSNP_index]
    r2.vals <- as.vector(cor(top.snp.G, G.mat))^2
    r2.brackets <- cut(r2.vals, breaks = c(0,0.1, 0.25, 0.75, 0.9, 1), labels = c("0-0.1","0.1-0.25","0.25-0.75","0.75-0.9","0.9-1"))
    pip.df$r2 <- r2.brackets
    
    ## actually make the track with Gviz
    pval.df <- pip.df[,c("chr","pos", "pval","r2")] %>% mutate(start = pos, end = pos) %>% dplyr::select(-pos) %>% pivot_wider(names_from = r2, values_from = "pval") 
    pval.df.gr <- makeGRangesFromDataFrame(pval.df, keep.extra.columns = T)
    seqlevelsStyle(pval.df.gr) <- "UCSC"
    pval.track <- DataTrack(range = pval.df.gr,  genome = "hg19", groups = names(mcols(pval.df.gr)), col = c("black","blue","green","orange","red"), name = "-log10 pvalue")
    
    # axis track - so we know where we are in the genome 
    axisTrack <- GenomeAxisTrack()

    #PIP TRACK
    pip.track <- DataTrack(data = pip.df$susie_pip, chromosome = pip.df$chr, start = pip.df$pos, end = pip.df$pos, genome = "hg19", name = "PIP")
    
    # gene track
    gene.track.track <- knownGeneObject(curr.locus.gr = curr.locus.gr, genome = "hg19")
    
    # hic and/or co-accessibility track, focused on the current region and only links with current gene 
    pchic.obj.filt <- pchic.obj[which(pchic.obj$anchor1.gene == locus.snp.dist$gene_name[l]),]
    hic.track <- InteractionTrack(pchic.obj.filt, name = "pcHiC")
    
    coacc.obj.filt <- coacc.obj[which(coacc.obj$anchor1.gene == locus.snp.dist$gene_name[l]),]
    coacc.track <- InteractionTrack(coacc.obj.filt, name = "Coaccess")
    
    dpars <- list(col.interactions="red", 
                  col.anchors.fill ="blue", 
                  col.anchors.line = "black", 
                  interaction.dimension = "height",  
                  interaction.measure = "counts", 
                  plot.trans = FALSE, 
                  plot.outside = FALSE,  
                  col.outside="lightblue",  
                  anchor.height = 0.1,
                  interaction.dimension.transform = "log10")
    displayPars(hic.track) <- dpars
    displayPars(coacc.track) <- dpars
    
    # put all tracks into a list
    list.of.tracks <- c(pval.track, pip.track, atac.tracks, h3k27.track, gene.track.track, hic.track, axisTrack)
    
    # highlight a particular SNP
    ht1 <- HighlightTrack(trackList = list.of.tracks, 
                          start = c(snp.p-500), width = 1000, 
                          chromosome = as.character(seqnames(curr.locus.gr)), col = 'pink')
    
    # plot
    pdf(paste0('manuscript_figures/figure5/',locus.snp.dist$gene_name[l],'_2_gviz_tracks_yaxis.pdf'), width=12, height=8)
    plotTracks(ht1,
               chromosome = as.character(seqnames(curr.locus.gr)), 
               transcriptAnnotation = "symbol", 
               collapseTranscripts= 'longest', 
               from = start(curr.locus.gr), 
               to = end(curr.locus.gr),
               sizes = c(1, 0.4, rep(0.2, 3), 0.1, 0.5, 1, 0.4),
               panel.only = F)
    dev.off()
    
}


