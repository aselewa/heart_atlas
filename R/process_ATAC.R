library(ArchR)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0){
  ATAC_DIR <- args[1]
} else{
  ATAC_DIR <- '/project2/gca/Heart_Atlas/ATAC_seq/'
  setwd('/project2/gca/aselewa/heart_atlas_project/')
}

# GLOBAL PARAMETERS
ATAC_SAMPLES <- c(paste0('MW200804A',c('A','B','C','D')),
                 paste0("SP-HE-HE200915ATAC-",c('175ATAC','359ATAC','360ATAC','397ATAC')),
                 paste0("SP-HE-MW200928E2ATAC-",c('175ATAC','366ATAC','367ATAC','407ATAC')),
                 paste0("SP-HE-MW200928E1ATAC-",c('396ATAC','398ATAC','406ATAC','408ATAC','411ATAC','413ATAC')))

ATAC_INDIVIDUALS <- c(rep("02207",4),rep("02336",4),rep("03231",4),rep("02336",2),rep("03231",2),rep("02207",2))

ATAC_REGIONS <- c("Septum","Right Atrium","Right Ventricle","Left Ventricle",
                 "Septum","Right Atrium","Right Ventricle","Left Ventricle",
                 "Septum","Right Atrium","Right Ventricle","Left Ventricle",
                 rep(c("Left Atrium","Apex"), 3))

make_meta_df <- function(){
  sat <- c()
  ncells <- c()
  for(s in ATAC_SAMPLES){
    curr.summary <- paste0(ATAC_DIR,s,'/outs/summary.json')
    meta <- fromJSON(txt = curr.summary)
    curr_sat <-  100*meta$bulk_estimated_saturation
    if(length(curr_sat)==0){
      curr_sat <- NA
    }
    sat <- c(sat, curr_sat)
    ncells <- c(ncells, meta$cells_detected)
  }
  meta.df <- data.frame(saturation=sat, ncells=ncells, sample=ATAC_SAMPLES, individual=RNA_INDIVIDUALS, regions=RNA_REGIONS)
  saveRDS(object = meta.df, file = 'ArchR_project_heart/ATAC_metadata_df.rds')
}

make_ATAC_arrows <- function(){

  inputFiles <- sapply(ATAC_SAMPLES, function(x){paste0(ATAC_DIR,x,'/outs/fragments.tsv.gz')})
  names(inputFiles) <- ATAC_SAMPLES

  addArchRThreads(threads = 5)
  addArchRGenome("hg38")
  system('mkdir ArchR_ArrowFiles')
  
  setwd('ArchR_ArrowFiles/')
  ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = names(inputFiles),
    filterTSS = 6,
    filterFrags = 3000,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
  )
  setwd('../')
}

main <- function() {
  
  make_meta_df()
  make_ATAC_arrows()
  
  individuals <- unique(ATAC_INDIVIDUALS)
  for(i in 1:length(individuals)){
    curr_samples <- ATAC_SAMPLES[ATAC_INDIVIDUALS == individuals[i]]
  
    projHeart <- ArchRProject(
      ArrowFiles = paste0("ArchR_ArrowFiles/",curr_samples, ".arrow"),
      outputDirectory = paste0("ArchR_project_",individuals[i]),
      copyArrows = TRUE 
    )  
    projHeart$individual <- individuals[i]
    projHeart$regions <- plyr::mapvalues(x = projHeart$Sample, from = curr_samples, to = ATAC_REGIONS[ATAC_INDIVIDUALS == individuals[i]])
    
    saveArchRProject(projHeart)
    
  }
  
}

if(length(args)>0){
  main()
}
