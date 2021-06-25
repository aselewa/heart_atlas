library(tidyverse)

gene.annots <- readRDS('genomic_annotations/hg38_gtf_genomic_annots.gr.rds')

# CM HiC data
cm.hic1 <- read_tsv('HiC/CM+iPSC_Hi-C/capt-CHiCAGO_interactions-CM-1.ibed') #hg19
cm.hic1$sample <- 1
cm.hic2 <- read_tsv('HiC/CM+iPSC_Hi-C/capt-CHiCAGO_interactions-CM-2.ibed') #hg19
cm.hic2$sample <- 2
cm.hic3 <- read_tsv('HiC/CM+iPSC_Hi-C/capt-CHiCAGO_interactions-CM-3.ibed') #hg19
cm.hic3$sample <- 3

cm.hic <- bind_rows(list(cm.hic1, cm.hic2, cm.hic3))

bait.gene.names <- sub(pattern = "\\*.*", replacement = "" , cm.hic$bait_name)
cm.hic.enhancer.gr <- GRanges(seqnames = cm.hic$otherEnd_chr, IRanges(start = cm.hic$otherEnd_start, end = cm.hic$otherEnd_end)) 
cm.hic.enhancer.gr$promoter_chr <- cm.hic$bait_chr
cm.hic.enhancer.gr$promoter_start <- cm.hic$bait_start
cm.hic.enhancer.gr$promoter_end <- cm.hic$bait_end
cm.hic.enhancer.gr$gene_name <- bait.gene.names
cm.hic.enhancer.gr$score <- cm.hic$score
cm.hic.enhancer.gr$sample <- cm.hic$sample
    
cm.hic.enhancer.protein.gr <- cm.hic.enhancer.gr[cm.hic.enhancer.gr$gene_name %in% gene.annots$genes$gene_name,] # restrict to protein coding genes

saveRDS(cm.hic.enhancer.protein.gr, 'HiC/iPSC_CM_pcHiC_protein_Hg19.gr.rds')
saveRDS(cm.hic.enhancer.gr, 'HiC/iPSC_CM_pcHiC_all_Hg19.gr.rds')

# ipsc HiC data
ipsc.hic1 <- read_tsv('HiC/CM+iPSC_Hi-C/capt-CHiCAGO_interactions-iPSC-1.ibed') #hg19
ipsc.hic1$sample <- 1
ipsc.hic2 <- read_tsv('HiC/CM+iPSC_Hi-C/capt-CHiCAGO_interactions-iPSC-2.ibed') #hg19
ipsc.hic2$sample <- 2
ipsc.hic3 <- read_tsv('HiC/CM+iPSC_Hi-C/capt-CHiCAGO_interactions-iPSC-2.ibed') #hg19
ipsc.hic3$sample <- 3

ipsc.hic <- bind_rows(list(ipsc.hic1, ipsc.hic2, ipsc.hic3))

bait.gene.names <- sub(pattern = "\\*.*", replacement = "" , ipsc.hic$bait_name)
ipsc.hic.enhancer.gr <- GRanges(seqnames = ipsc.hic$otherEnd_chr, IRanges(start = ipsc.hic$otherEnd_start, end = ipsc.hic$otherEnd_end)) 
ipsc.hic.enhancer.gr$promoter_chr <- ipsc.hic$bait_chr
ipsc.hic.enhancer.gr$promoter_start <- ipsc.hic$bait_start
ipsc.hic.enhancer.gr$promoter_end <- ipsc.hic$bait_end
ipsc.hic.enhancer.gr$gene_name <- bait.gene.names
ipsc.hic.enhancer.gr$score <- ipsc.hic$score
ipsc.hic.enhancer.gr$sample <- ipsc.hic$sample

ipsc.hic.enhancer.protein.gr <- ipsc.hic.enhancer.gr[ipsc.hic.enhancer.gr$gene_name %in% gene.annots$genes$gene_name,] # restrict to protein coding genes

saveRDS(ipsc.hic.enhancer.protein.gr, 'HiC/iPSC_pcHiC_protein_Hg19.gr.rds')
saveRDS(ipsc.hic.enhancer.gr, 'HiC/iPSC_pcHiC_all_Hg19.gr.rds')


# 
ipsc.hic.enhancer.gr <- readRDS('HiC/iPSC_pcHiC_all_Hg19.gr.rds')
cm.hic.enhancer.gr <-  readRDS('HiC/iPSC_CM_pcHiC_all_Hg19.gr.rds')

ipsc.hic.promoter.gr <- GRanges(seqnames = ipsc.hic.enhancer.gr$promoter_chr, ranges = IRanges(start = ipsc.hic.enhancer.gr$promoter_start, end = ipsc.hic.enhancer.gr$promoter_end),
                          sample = ipsc.hic.enhancer.gr$sample)
cm.hic.promoter.gr <- GRanges(seqnames = cm.hic.enhancer.gr$promoter_chr, ranges = IRanges(start = cm.hic.enhancer.gr$promoter_start, end = cm.hic.enhancer.gr$promoter_end),
                          sample = cm.hic.enhancer.gr$sample)


enhancer.hits <- GenomicRanges::findOverlaps(query = cm.hic.enhancer.gr, subject = ipsc.hic.enhancer.gr, maxgap = 1000)
promoter.hits <- GenomicRanges::findOverlaps(query = cm.hic.promoter.gr, subject = ipsc.hic.promoter.gr, maxgap = 1000)

ehits.str <- paste0(queryHits(enhancer.hits),'-',subjectHits(enhancer.hits))
phits.str <- paste0(queryHits(promoter.hits),'-',subjectHits(promoter.hits))

ehitsIn <- queryHits(enhancer.hits)
ehitsIn <- ehitsIn[ehits.str %in% phits.str]
ehitsIn <- unique(ehitsIn)

cm.hic.enhancer.gr <- cm.hic.enhancer.gr[-ehitsIn,]

saveRDS(cm.hic.enhancer.gr, 'HiC/CM_Only_pcHiC_all_Hg19.gr.rds')

