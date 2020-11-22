library(Seurat)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0){
  RNA_DIR <- args[1]
} else{
  RNA_DIR <- '/project2/gca/Heart_Atlas/Nuc_seq/'
  setwd('/project2/gca/aselewa/heart_atlas_project/')
}

# GLOBAL PARAMETERS
# GLOBAL PARAMETERS
RNA_SAMPLES <- c("MW200804RA","MW200804RB","MW200804RC","MW200804RD","SP-HE-MW200928E1RNA-411RNA","SP-HE-MW200928E1RNA-413RNA",
                  "SP-HE-HE200915RNA-175RNA","SP-HE-HE200915RNA-359RNA","SP-HE-HE200915RNA-360RNA","SP-HE-HE200915RNA-397RNA","SP-HE-MW200928E1RNA-396RNA","SP-HE-MW200928E1RNA-398RNA",
                  "SP-HE-MW200928E2RNA-175RNA", "SP-HE-MW200928E2RNA-366RNA","SP-HE-MW200928E2RNA-367RNA","SP-HE-MW200928E2RNA-407RNA","SP-HE-MW200928E1RNA-406RNA","SP-HE-MW200928E1RNA-408RNA")

RNA_INDIVIDUALS <- c(rep("02207",6),rep("02336",6),rep("03231",6))

RNA_REGIONS <- rep(c("Septum","Right Atrium","Right Ventricle","Left Ventricle","Left Atrium","Apex"), 3)

process_RNA <- function(){
  
  print("Loading RNA samples...")
  
  # Load RNA data
  obj.list <- list()
  sat <- c()
  ncells <- c()
  for(i in 1:length(RNA_SAMPLES)){
    s <- RNA_SAMPLES[i]
    ind <- RNA_INDIVIDUALS[i]
    region <- RNA_REGIONS[i]
    
    curr_dirs <- list.dirs(paste0(RNA_DIR,s))
    genefull_dir <- curr_dirs[grepl(pattern = "/GeneFull$", x = curr_dirs)]
    genefull_raw_dir <- curr_dirs[grepl(pattern = "GeneFull/raw", x = curr_dirs)]
    
    meta <- readr::read_csv(paste0(genefull_dir,"/","Summary.csv"), col_names = F)
    sat <- c(sat, 100*meta$X2[3])
    ncells <- c(ncells, meta$X2[10])

    curr <- Read10X(genefull_raw_dir)
    curr_obj <- CreateSeuratObject(curr, min.features = 10, min.cells = 10, assay = "RNA")
    curr_obj$individual <- ind
    curr_obj$dataset <- s
    curr_obj$region <- region
    obj.list[[s]] <- curr_obj
  }
  
  first <- obj.list[[1]]
  obj.list <- obj.list[-1]
  heart_rna <- merge(x = first, y = obj.list)
  heart_rna$percent.mito <- PercentageFeatureSet(heart_rna, pattern = '^MT-')
  
  saveRDS(heart_rna, file = "seurat/Heart_Seurat_RNA_all_samples_raw.RDS")
  
  meta.df <- data.frame(saturation=sat, ncells=ncells, sample=RNA_SAMPLES, individual=RNA_INDIVIDUALS, regions=RNA_REGIONS)
  saveRDS(object = meta.df, file = 'seurat/RNA_seq_metadata.rds')
  
}

# if running from command line/Makefile
if(length(args) > 0){
  process_RNA()
} 
