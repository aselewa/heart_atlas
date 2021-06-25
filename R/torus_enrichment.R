library(finemappeR)
library(liftOver)
library(ComplexHeatmap)
source('R/analysis_utils.R')
setwd('/project2/gca/aselewa/heart_atlas_project/')
celltype_ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")

#gwas.sumstats <- list.files(path = 'GWAS/summary_statistics', pattern = '*.df.rds', recursive = T, full.names = T)
trait <- basename(dirname(gwas.sumstats))
names(gwas.sumstats) <- trait

rds.list <- lapply(
  trait,
  function(x){
    readRDS(gwas.sumstats[x])
  }
)

# prepare annotations
markers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')
path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch <- import.chain(path)
markers.hg19 <- lapply(markers, function(x){unlist(liftOver(x, ch))})
for(i in 1:length(markers.hg19)){
    seqlevelsStyle(markers.hg19[[i]]) <- "NCBI"
}
system('mkdir -p GWAS/bed_annotations_hg19')
lapply(names(markers.hg19), function(x){rtracklayer::export(markers.hg19[[x]], 
                                                            format = 'bed', 
                                                            con = paste0('GWAS/bed_annotations_hg19/',x,'_narrowPeaks.bed'))})
for(i in 1:length(markers.hg19)){
    seqlevelsStyle(markers.hg19[[i]]) <- "UCSC"
}
system('mkdir -p GWAS/bed_annotations_hg19')
lapply(names(markers.hg19), function(x){rtracklayer::export(markers.hg19[[x]], 
                                                            format = 'bed', 
                                                            con = paste0('GWAS/bed_annotations_hg19/UCSC_',x,'_narrowPeaks.bed'))})

annotations <- list.files(path = 'GWAS/bed_annotations_hg19/', pattern = '*.bed', full.names = T)

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

saveRDS(enrich.res, 'GWAS/Torus_CellType_Enrichment_Results_Univariate_MORE.df.rds')
enrich.res <- readRDS('GWAS/Torus_CellType_Enrichment_Results_Univariate_MORE.df.rds')

pval_from_ci <- function(mean, upper, ci){
    
    nsamp <- length(mean)
    pval.res <- rep(0, nsamp)
    for(i in 1:nsamp){
        alph <- (1-ci)/2
        zval <- qnorm(p = 1-alph)
        se <- (upper[i]-mean[i])/zval
        
        pval.res[i] <- 1 - pnorm(q = mean[i] / se)   
    }
    pval.res
}

res <- lapply(enrich.res, function(x){ Reduce(x = x, f = rbind)})
res <- lapply(res, function(x){x[x$term != "Intercept",]})
for(i in 1:length(res)){
    res[[i]]$pvalue <- pval_from_ci(mean = res[[i]]$estimate, upper = res[[i]]$high, ci = 0.95)
}
estimates <- as.data.frame(sapply(res, function(x){x["estimate"]}))
pvalues <- as.data.frame(sapply(res, function(x){x["pvalue"]}))
fdr <- matrix(p.adjust(unlist(pvalues), method = 'BH'), nrow = nrow(pvalues))

rnames <- basename(annotations)
names.order <- c("aFib", "PR_Interval","heart_rate","heart_failure",
                 "CAD","DiastolicBP","asthma","BMI","Height")
#celltype_ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")
celltype_ideal_order <- c("Cardiomyocyte","Pericyte","Endothelial","Fibroblast")


row.names(estimates) <- sub('_narrowPeaks.bed','',rnames)
colnames(estimates) <- names(enrich.res)
estimates <- estimates[celltype_ideal_order,names.order]
estimates <- t(estimates)

row.names(fdr) <- sub('_narrowPeaks.bed','',rnames)
colnames(fdr) <- names(enrich.res)
fdr <- fdr[celltype_ideal_order,names.order]
fdr <- t(fdr)

star.mat <- matrix('ns', nrow = nrow(fdr), ncol = ncol(fdr))
star.mat[fdr < 0.05] <- '*'
star.mat[fdr < 0.0001] <- '***'
rownames(star.mat) <- rownames(fdr)
colnames(star.mat) <- colnames(fdr)
    
mat.to.viz <- estimates/log(2)
mat.to.viz[mat.to.viz < 0] <- 0


### plotting

lgd_list <- list()

col_fun <- c("lightblue", "orange", "firebrick")
names(col_fun) <- c("ns", '*', '***')

lgd_list[["fdr"]] <- Legend(title = "fdr (binned)",
                            labels = c("ns", '*', '***'),
                            legend_gp = gpar(fill = col_fun))

tic_vec <- c(0, 2, 4)
lgd_list[["log2_enrich"]] <- Legend(title = "log2_enrich",
                                    labels = tic_vec,
                                    # labels_gp = gpar(fontsize = 14),
                                    grid_height = unit(6, "mm"),
                                    grid_width = unit(6, "mm"),
                                    graphics = list(
                                        function(x, y, w, h) grid.circle(x, y, r = (tic_vec[1]/10 + 0.2) * unit(2.5, "mm"),
                                                                         gp = gpar(fill = "black")),
                                        function(x, y, w, h) grid.circle(x, y, r = (tic_vec[2]/10 + 0.2) * unit(2.5, "mm"),
                                                                         gp = gpar(fill = "black")),
                                        function(x, y, w, h) grid.circle(x, y, r = (tic_vec[3]/10 + 0.2) * unit(2.5, "mm"),
                                                                         gp = gpar(fill = "black"))
                                    ))

map1 <- Heatmap(star.mat,
                name = "Association Effect Size",
                col = col_fun,
                rect_gp = gpar(type = "none"),
                cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width = width, height = height, 
                              gp = gpar(col = NA, fill = NA))
                    grid.circle(x = x, y = y,
                                r = (mat.to.viz[i, j]/10 + 0.2) * unit(2.5, "mm"),
                                gp = gpar(fill = col_fun[star.mat[i, j]], col = NA))
                },
                border_gp = gpar(col = "black"),
                row_title = "Trait",
                column_title = "Cell Type",
                cluster_rows = F, cluster_columns = F,
                show_heatmap_legend = F,
                row_names_gp = gpar(fontsize = 10.5),
                column_names_rot = 45,
                column_names_side = "top", 
                use_raster = T)

pdf('manuscript_figures/figure5/Fig5A_torus_enrichment_clean.pdf', width=5, height=4)
draw(map1, annotation_legend_list = lgd_list)
dev.off()

