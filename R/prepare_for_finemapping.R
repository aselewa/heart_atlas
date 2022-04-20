library(finemappeR)
library(liftOver)
library(ComplexHeatmap)
source('R/analysis_utils.R')
setwd('/project2/gca/aselewa/heart_atlas_project/')

# Run torus on joint model
rds <- readRDS('GWAS/summary_statistics/aFib/ebi-a-GCST006414_aFib.df.rds')
annotations <- list.files(path = 'GWAS/annotations_for_finemapping_hg19/', pattern = '*.bed', full.names = T)

cleaned.gwas.annots <- annotator(rds, annotations = annotations)

annot_columns = colnames(cleaned.gwas.annots)[endsWith(x=colnames(cleaned.gwas.annots), suffix='bed_d')]
cleaned.gwas.annots = cleaned.gwas.annots[,c("snp",annot_columns)]

# get SNP prior
readr::write_tsv(x = cleaned.gwas.annots, path = 'torus_annotations.txt.gz', col_names = T)
readr::write_tsv(x = rds[,c('snp','locus','zscore')], path = 'torus_zscore.txt.gz', col_names = T)
torus.result <- RunTorus(torus_annot_file = 'torus_annotations.txt.gz', torus_zscore_file = 'torus_zscore.txt.gz') 

saveRDS(torus.result, file = 'GWAS/Torus_Enrichment_Results_Joint.rds')
snp_torus_pip <- torus.result$snp_pip

# get locus FDR
torus.result.fdr <- RunTorusFDR(torus_annot_file = 'torus_annotations.txt.gz', torus_zscore_file = 'torus_zscore.txt.gz') 

system('rm torus_annotations.txt.gz')
system('rm torus_zscore.txt.gz')

# add Torus PIPs to normalized GWAS and split by locus
sumstats_finemap <- PrepareSusieData(sumstats,torus_pip = snp_torus_pip, torus_fdr = torus.result.fdr, fdr_thresh = 0.1)

prnt('Will finemap this many loci:')
print(length(unique(sumstats_finemap$locus)))

# save each loci to its own RDS
loci <- unique(sumstats_finemap$locus)
for(l in loci){
  curr.sumtats <- sumstats_finemap[sumstats_finemap$locus == l,]
  saveRDS(curr.sumtats, file=paste0('GWAS/finemapping/locus_files_finemapping/aFib_SignifLoci_',l,'.df.rds'))
}

# run finemapping on sbatch using R/run_finemapping.R script on each locus saved above

