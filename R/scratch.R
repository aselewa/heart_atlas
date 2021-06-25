library(tidyverse)
setwd('/project2/gca/aselewa/heart_atlas_project/')
# processing ENCODE batch downloads
fname.map <- metad %>% filter(`File format` %in% "bed narrowPeak") %>% select(`File accession`, `Biosample term name`) %>% write_tsv('ENCODE/H3K27ac_extensive/file_name_biosample.tsv')
fname.map2 <- fname.map %>% group_by(`Biosample term name`) %>% mutate(bedfiles = paste0(`File accession`,'.bed.gz')) %>% summarise(og_fnames = paste0(bedfiles, collapse = ' '))

for(i in 1:nrow(fname.map2)){ 
    fname=fname.map2$og_fnames[i]
    new_fname=gsub(' ', '_', paste0(fname.map2$`Biosample term name`[i],'.bed.gz'), fixed = T)
    cmd = paste0('cat ',fname,' > ',new_fname)
    system(cmd)
}

