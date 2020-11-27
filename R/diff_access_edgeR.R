library(ArchR)
library(edgeR)

setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

satac <- loadArchRProject('ArchR/ArchR_project_03231/')
peak.info <- getPeakSet(satac)
peak.names <- GRToString(peak.info)

donor.chamber.celltype <- paste0(satac$individual,'_',satac$regions,'_',satac$CellTypes)
satac$agg <- donor.chamber.celltype
agg.insert.mat <- getGroupSE(ArchRProj = satac, useMatrix = "PeakMatrix", groupBy = "agg", divideN = F)
agg.insert.mat <- agg.insert.mat@assays@data$PeakMatrix
row.names(agg.insert.mat) <- peak.names

cell.types <- unique(satac$CellTypes)
gr.list <- list()

for(c in cell.types){
  print(c)
  condition <- factor(1*(grepl(colnames(agg.insert.mat), pattern = c)))
  edgeR_obj <- DGEList(counts = agg.insert.mat, group = condition)
  edgeR_obj <- calcNormFactors(edgeR_obj)
  design_obj <- model.matrix(~ condition) 
  
  edgeR_obj <- estimateDisp(edgeR_obj, design = design_obj)
  fit <- glmFit(edgeR_obj, design = design_obj)
  lrt <- glmLRT(fit, coef = 2)
  toptag <- topTags(lrt, n = Inf, adjust.method = "BH")
  toptagtable <- toptag$table
  
  signif.gr <- StringToGR(rownames(toptagtable))
  signif.gr$FDR <- toptagtable$FDR
  signif.gr$Log2FC <- toptagtable$logFC *log2(exp(1))
  signif.gr$Pvalue <- toptagtable$PValue
  gr.list[[c]] <- signif.gr
}

saveRDS(gr.list, 'ArchR/ArchR_project_03231/PeakCalls/edgeR_DA_peaks.gr.rds')





