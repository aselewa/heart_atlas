library(finemappeR)

setwd('/project2/gca/aselewa/heart_atlas_project/')

gwas.sumstats <- list.files(path = 'GWAS/summary_statistics', pattern = '*.df.rds', recursive = T, full.names = T)
trait <- basename(dirname(gwas.sumstats))
names(gwas.sumstats) <- trait

rds.list <- sapply(
  trait,
  function(x){
    readRDS(gwas.sumstats[x])
  }
)

annotations <- list.files(path = 'GWAS/annotations_hg19', pattern = '*.bed', full.names = T)

enrich.res <- list()

for(z in seq_along(rds.list)){
  enrich.res[[trait[z]]] <- list()
  for(w in seq_along(annotations)){
    
    cleaned.gwas.annots <- annotator(rds.list[[z]], annotations = annotations[w])
  
    readr::write_tsv(x = cleaned.gwas.annots[,-c(1:6,8:12)], path = 'torus_annotations.txt.gz', col_names = T)
    readr::write_tsv(x = rds.list[[z]][,c('snp','locus','zscore')], path = 'torus_zscore.txt.gz', col_names = T)
    
    torus.result <- RunTorus(torus_annot_file = 'torus_annotations.txt.gz', torus_zscore_file = 'torus_zscore.txt.gz') 
    
    enrich.res[[trait[z]]][[basename(annotations[w])]] <- torus.result$enrich
  }
}

saveRDS(enrich.res, 'GWAS/Torus_Enrichment_Results_Univariate.df.rds')
