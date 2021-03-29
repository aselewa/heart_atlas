library(ArchR)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0){
  ATAC_DIR <- args[1]
} else{
  ATAC_DIR <- '/project2/gca/Heart_Atlas/ATAC_seq/'
  setwd('/project2/gca/aselewa/heart_atlas_project/ArchR/')
}

# GLOBAL PARAMETERS
ATAC_SAMPLES <- c("MW200804AA","Deep_MW200804AC","Deep_MW200804AD","SP-HE-MW200928E1ATAC-413ATAC",
                 "SP-HE-HE200915ATAC-175ATAC","Deep_SP-HE200915ATAC-360ATAC","Deep_SP-HE200915ATAC-397ATAC","SP-HE-MW200928E1ATAC-398ATAC",
                 "SP-HE-MW200928E2ATAC-175ATAC","Deep_SP-MW200928E2ATAC-367ATAC","Deep_SP-MW200928E2ATAC-407ATAC","SP-HE-MW200928E1ATAC-408ATAC")

ATAC_INDIVIDUALS <- c(rep("02207",4),rep("02336",4),rep("03231",4))

ATAC_REGIONS <- rep(c("Septum","Right Ventricle","Left Ventricle","Apex"), 3)

make_ATAC_arrows <- function(){

  inputFiles <- sapply(ATAC_SAMPLES, function(x){paste0(ATAC_DIR,x,'/outs/fragments.tsv.gz')})
  names(inputFiles) <- ATAC_SAMPLES

  addArchRThreads(threads = 5)
  addArchRGenome("hg38")
  system('mkdir ArchR_ArrowFiles_5kmin_6TSSmin')
  
  setwd('ArchR_ArrowFiles_5kmin_6TSSmin/')
  ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = names(inputFiles),
    filterTSS = 6,
    filterFrags = 5000
  )
  
  arrows <- list.files(path = '.', pattern = '*.arrow', full.names = F)
  x <- addDoubletScores(input = arrows, k = 5)
  
  setwd('../')
  
}

main <- function() {

  make_ATAC_arrows()
  
  individuals <- unique(ATAC_INDIVIDUALS)
  for(i in 1:length(individuals)){
    curr_samples <- ATAC_SAMPLES[ATAC_INDIVIDUALS == individuals[i]]
    projDir <- paste0("ArchR/ArchR_project_",individuals[i],'_3')
    
    projHeart <- ArchRProject(
       ArrowFiles = paste0("ArchR/ArchR_ArrowFiles/",curr_samples, ".arrow"),
       outputDirectory = projDir,
       copyArrows = TRUE 
    )
    projHeart$individual <- individuals[i]
    projHeart$regions <- plyr::mapvalues(x = projHeart$Sample, from = curr_samples, to = ATAC_REGIONS[ATAC_INDIVIDUALS == individuals[i]])
    
    saveArchRProject(projHeart)
    
    frag.sizes <- getFragmentSizes(projHeart)
    saveRDS(frag.sizes, paste0(projDir, '/fragment_sizes.rds'))
    
  }
}

if(length(args)>0){
  main()
}


