library(finemappeR)

bigSNP <- bigsnpr::snp_attach(rdsfile = '/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds')
gwas <- RunCleaner(sumstats = 'GWAS/summary_statistics/', 
                   ColsToKeep = c('CHR','POS_GRCh37','Effect_A2','StdErr','A1','A2','MarkerName','Pvalue'), 
                   bigSNP = bigSNP)
saveRDS('GWAS/summary_statistics/nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.df.rds')

gwas <- readRDS('GWAS/summary_statistics/nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.df.rds')

annotations <- list.files(path = '../GWAS/annotations_hg19/', pattern = '*.bed', full.names = T)
enrich.list <- list()

for(i in 1:length(annotations)){
  print(annotations[i])
  
  system(paste0('mkdir -p ../GWAS/annotation_bed; cp ', annotations[i] ,' ../GWAS/annotation_bed/'))
  
  torus.files <- PrepareTorusFiles(gwas, bed_annotations = '../GWAS/annotation_bed/')
  torus.result <- RunTorus(torus_annot_file = torus.files[[1]], torus_zscore_file = torus.files[[2]])
  
  system('rm -rf ../GWAS/annotation_bed')
  system(paste0('rm ', torus.files[[1]]))
  system(paste0('rm ', torus.files[[2]]))
  
  enrich.list[[i]] <- torus.result$enrich
}

enrich.df <- Reduce(rbind, enrich.list)
enrich.df <- enrich.df[enrich.df$term!="Intercept",]
saveRDS(enrich.df, file='../GWAS/Torus_Enrichment_Results_Univariate.df.rds')
