
setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

gwas <- readRDS('GWAS/finemapping/aFib_Finemapped.tble.rds')
eqtl <- vroom::vroom('eQTL_enrich/summary_statistics/Heart_LV_eQTL_b37.txt.gz')
bigsnp <- bigsnpr::snp_attach('/project2/xinhe/1kg/bigsnpr/EUR.1kg.rds')

ldblocks <- finemappeR::LD_Blocks
ldblocks.gr <- GRanges(seqnames = ldblocks$X1, ranges = IRanges(start = ldblocks$X2, end = ldblocks$X3), num = ldblocks$X4)
varids <- paste0(gwas$chr,'_',gwas$pos,'_',gwas$a0,'_',gwas$a1,'_b37')
gwas$varid <- varids

loci <- unique(gwas$locus)
coloc.results <- list()
for(l in loci){
    print(l)
    
    curr.gwas <- gwas[gwas$locus==l,]
    curr.varids <- curr.gwas$varids
    
    curr.eqtl <- eqtl[eqtl$variant_id %in% curr.varids,]
    
    if(nrow(curr.eqtl) > 5){
        curr.eqtl <- curr.eqtl %>% 
            mutate(zscore = slope/slope_se) %>% 
            group_by(variant_id) %>% 
            arrange(-abs(zscore)) %>% 
            slice(1)
        
        
        curr.gwas <- curr.gwas[match(curr.eqtl$variant_id, curr.gwas$varids),]
        names <- paste0("s",1:ncol(Rgwas))
        
        X <- bigsnp$genotypes[ , curr.gwas$bigSNP_index]
        X <- scale(X, center = T, scale = T)
        zhat <- curr.gwas$zscore
        Rgwas <- cov2cor((crossprod(X) + tcrossprod(zhat))/nrow(X))
        colnames(Rgwas) <- curr.gwas$varids
        rownames(Rgwas) <- curr.gwas$varids
        
        zhat <- curr.eqtl$zscore
        Reqtl <- cov2cor((crossprod(X) + tcrossprod(zhat))/nrow(X))
        colnames(Reqtl) <- curr.gwas$varids
        rownames(Reqtl) <- curr.gwas$varids
        
        gwas.obj <- list(beta = curr.gwas$beta, varbeta = curr.gwas$se^2, snp = curr.gwas$varids, position = curr.gwas$pos, type="cc", sdY=1, LD=Rgwas)
        eqtl.obj <- list(beta = curr.eqtl$slope, varbeta = curr.eqtl$slope_se^2, snp = curr.gwas$varids, position = curr.gwas$pos, type="quant", sdY=1, LD=Reqtl)
        
        #coloc.results[[`l`]] <- coloc.abf(dataset1=gwas.obj, dataset2=eqtl.obj) 
        coloc.results[as.character(l)] <- coloc.susie(dataset1 = gwas.obj, dataset2 = eqtl.obj, susie.args = list(prior_weights = curr.gwas$torus_pip, L=1))
    }
    
}

coloc.results[sapply(coloc.results, is.null)] <- NULL
saveRDS(coloc.results, file = 'misc/coloc_analysis.rds')

h0.probability <- sapply(coloc.results, function(x){x$summary['PP.H0.abf']*100})
h3.probability <- sapply(coloc.results, function(x){x$summary['PP.H3.abf']*100})
h4.probability <- sapply(coloc.results, function(x){x$summary['PP.H4.abf']*100})

par(mfrow=c(1,2))
hist(100-h0.probability, breaks = 30, main = 'PP atleast one trait has genetic assoc.', xlab='Probability', ylab='Number of Loci')
hist(h4.probability, breaks = 30, main = 'PP both traits share a single causal variant', xlab='Probability', ylab='Number of Loci')
