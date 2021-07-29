
library(coloc)
setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

finemap.res <- read_csv('GWAS/aFib_Finemapped_minPIP_0.2_06252021.csv')
gwas <- readRDS('GWAS/finemapping/aFib_Finemapped.tble.rds')
eqtl <- readRDS('misc/V8_Signif_eQTLs_lifted_hg19.rds') %>% as_tibble()

ldblocks <- finemappeR::LD_Blocks
ldblocks.gr <- GRanges(seqnames = ldblocks$X1, ranges = IRanges(start = ldblocks$X2, end = ldblocks$X3), num = ldblocks$X4)
varids <- paste0(gwas$chr,'_',gwas$pos,'_b37')
gwas$varids <- varids

loci <- unique(finemap.res[!is.na(finemap.res$eQTL_Symbols),]$locus)
coloc.results <- list()
for(l in loci){
    print(l)
    
    curr.gwas <- gwas[gwas$locus==l,]
    curr.varids <- curr.gwas$varids
    
    causal.snps <- finemap.res$SNP[finemap.res$locus==l]
    causal.snp.ids <- curr.varids[curr.gwas$snp %in% causal.snps]
    curr.genes <- finemap.res$eQTL_Symbols[finemap.res$locus==l]
    
    curr.eqtl <- eqtl[eqtl$variant_ids %in% curr.varids,]

    curr.gwas <- curr.gwas[match(curr.eqtl$variant_ids, curr.gwas$varids),]
    
    curr.eqtl$variant_ids <- make.unique(curr.eqtl$variant_ids)
    curr.gwas$varids <- curr.eqtl$variant_ids

    gwas.obj <- list(beta = curr.gwas$beta, varbeta = curr.gwas$se^2, snp = curr.gwas$varids, position = curr.gwas$pos, type="cc", sdY=1)
    eqtl.obj <- list(beta = curr.eqtl$slope, varbeta = curr.eqtl$slope_se^2, snp = curr.gwas$varids, position = curr.gwas$pos, type="quant", sdY=1)
    
    coloc.results[[`l`]] <- coloc.abf(dataset1=gwas.obj, dataset2=eqtl.obj) 
    #coloc.results[as.character(l)] <- coloc.susie(dataset1 = gwas.obj, dataset2 = eqtl.obj, susie.args = list(prior_weights = curr.gwas$torus_pip, L=1))
}

coloc.results[sapply(coloc.results, is.null)] <- NULL
saveRDS(coloc.results, file = 'misc/coloc_analysis.rds')

h0.probability <- sapply(coloc.results, function(x){x$summary['PP.H0.abf']*100})
h3.probability <- sapply(coloc.results, function(x){x$summary['PP.H3.abf']*100})
h4.probability <- sapply(coloc.results, function(x){x$summary['PP.H4.abf']*100})

hist(h4.probability, breaks = 30, main = 'PP both traits share a single causal variant', xlab='Probability', ylab='Number of Loci')

coloc.results.high <- coloc.results[h4.probability > 50]
shared.snps <- list()
for(i in 1:length(coloc.results.high)){
  
  curr.locus.df <- coloc.results.high[[i]]$results %>% 
    arrange(-SNP.PP.H4) %>% 
    filter(SNP.PP.H4 > 0.1) %>%
    dplyr::select(snp, SNP.PP.H4)
  
  shared.snps[[i]] <- curr.locus.df
  
}


