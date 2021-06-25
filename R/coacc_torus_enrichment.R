setwd('/project2/gca/aselewa/heart_atlas_project/')
library(tidyverse)
source('R/analysis_utils.R')

lv.sumstats <- vroom::vroom('eQTL_enrich/summary_statistics/Heart_LV_eQTL_sumstats.txt.gz')

gene.annots <- readRDS('genomic_annotations/hg38_gtf_genomic_annots.gr.rds')
gene.annots <- gene.annots$genes

satac <- suppressMessages(ArchR::loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/'))
peak.set <- satac@peakSet
peak.set$peakID <- GRToString(peak.set)

coacc.enhancers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coaccess_enhancers_corr_cut_0_maxDist_1Mb_hg38.gr.rds')
coacc.enhancers.high <- coacc.enhancers[coacc.enhancers$correlation > 0.5,]
coacc.enhancers.low <- coacc.enhancers[coacc.enhancers$correlation < 0.1]

coacc.genes <- unique(coacc.enhancers.high$coacc_gene_name)
coacc.enhancers.low <- coacc.enhancers.low[coacc.enhancers.low$coacc_gene_name %in% coacc.genes,]

gene.ids <- gene.annots[gene.annots$gene_name %in% coacc.genes,]$gene_id

gene.ids.noDot <- sub(pattern = '\..*', replacement = '', x = gene.ids)
gene.ids.noDot <- sub(pattern = '\\..*', replacement = '', x = gene.ids)
gtex.gene.ids.noDot <- sub('\\..*', '', lv.sumstats$gene_id)

lv.sumstats.sub <- lv.sumstats[gtex.gene.ids.noDot %in% gene.ids.noDot,]
lv.sumstats.sub %>% write_tsv(path = 'eQTL_enrich/summary_statistics/Heart_LV_eQTL_CoaccGenes_sumstats.txt.gz')

# prepare snp map and gene map

snp.map <- readRDS('eQTL_enrich/metadata/hg38_SNP_map.gr.rds')

snp.map2 <- snp.map[snp.map$snp %in% lv.sumstats.sub$variant_id,]
snp.map2 %>% as_tibble()  %>% select(snp,seqnames, start) %>% write_tsv('eQTL_enrich/metadata/hg38_snp_map_Coacc.txt.gz', col_names = FALSE)

gene.map <- vroom::vroom('eQTL_enrich/metadata/hg38_gene_tss_map.txt.gz', col_names = F)
gene.map2 <- gene.map[gene.map$X1 %in% lv.sumstats.sub$gene_id,]
gene.map2 %>% write_tsv('eQTL_enrich/metadata/hg38_gene_tss_map_Coacc_genes.txt.gz', col_names=F)

# prepare annotations

coacc.enhancers.list <- list("0-20kb" = coacc.enhancers.high[coacc.enhancers.high$distToCoaccGene < 20000,],
                             "20kb-100kb" = coacc.enhancers.high[coacc.enhancers.high$distToCoaccGene >= 20000 & coacc.enhancers.high$distToCoaccGene < 100000,],
                             ">100kb" = coacc.enhancers.high[coacc.enhancers.high$distToCoaccGene >= 100000])

sapply(names(coacc.enhancers.list),
       function(x){
           snp.map.coacc <- IRanges::subsetByOverlaps(snp.map2, coacc.enhancers.list[[x]])
           snpsIn <- unique(snp.map.coacc$snp)
           snp.map2$coacc <- 1*(snp.map2$snp %in% snpsIn)
           
           snp.map.df <- snp.map2 %>% as_tibble() %>% select(snp, coacc) %>% rename(SNP = snp)
           vroom::vroom_write(snp.map.df, path = paste0('eQTL_enrich/annotations/Coaccess_',x,'_overlapped_annot.txt.gz'), delim = '\t')
       })



control.peak.list <- list("0-20kb" = coacc.enhancers.low[coacc.enhancers.low$distToCoaccGene < 20000,],
                             "20kb-100kb" = coacc.enhancers.low[coacc.enhancers.low$distToCoaccGene >= 20000 & coacc.enhancers.low$distToCoaccGene < 100000,],
                             ">100kb" = coacc.enhancers.low[coacc.enhancers.low$distToCoaccGene >= 100000])

sapply(names(control.peaks.list),
       function(x){
           snp.map.coacc <- IRanges::subsetByOverlaps(snp.map2, control.peaks.list[[x]])
           snpsIn <- unique(snp.map.coacc$snp)
           snp.map2$control <- 1*(snp.map2$snp %in% snpsIn)
           
           snp.map.df <- snp.map2 %>% as_tibble() %>% select(snp, control) %>% rename(SNP = snp)
           vroom::vroom_write(snp.map.df, path = paste0('eQTL_enrich/annotations/Control_',x,'_overlapped_annot.txt.gz'), delim = '\t')
       })

