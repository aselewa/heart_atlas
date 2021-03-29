library(ArchR)
library(tidyverse)
library(liftOver)
source('R/analysis_utils.R')
setwd('/project2/gca/aselewa/heart_atlas_project/')

satac <- ArchR::loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/')
coacc <- ArchR::getCoAccessibility(satac, corCutOff = 0.001, resolution = 1, returnLoops = F) %>% as_tibble()
gene.annots <- ArchR::geneAnnoHg38$genes
geneStart <- gene.annots %>% 
    as_tibble() %>% 
    mutate(geneStart = ifelse(strand(gene.annots)=="+",start(gene.annots),end(gene.annots))) %>% 
    select(geneStart, symbol) %>% 
    filter(!is.na(symbol)) %>%
    rename(coacc_gene_name = symbol)

peak.granges <- ArchR::getPeakSet(satac)
coacc$queryType <- peak.granges[coacc$queryHits,]$peakType
coacc$subjectType <- peak.granges[coacc$subjectHits,]$peakType

coacc <- coacc[coacc$queryHits < coacc$subjectHits,]


# keep links that involve a promoter, make queryHits be non-promoter peaks.
coacc <- coacc[xor(coacc$queryType == "Promoter", coacc$subjectType == "Promoter"),]
needSwap <- which(coacc$queryType == "Promoter")
temp <- coacc[needSwap, "subjectHits"] 
coacc[needSwap, "subjectHits"] <- coacc[needSwap, "queryHits"]
coacc[needSwap, "queryHits"] <- temp

other.ranges <- peak.granges[coacc$queryHits,]
other.ranges$correlation <- coacc$correlation
promoter.ranges <- peak.granges[coacc$subjectHits,]

isNaGene <- is.na(promoter.ranges$nearestGene)
other.ranges <- other.ranges[!isNaGene]
promoter.ranges <- promoter.ranges[!isNaGene]

other.ranges$distToPromoterPair <- abs(start(other.ranges) - start(promoter.ranges))
other.ranges$coacc_gene_name <- promoter.ranges$nearestGene
other.ranges$distToCoaccGene <- other.ranges %>% as_tibble() %>% 
    left_join(., geneStart, on="coacc_gene_name") %>% 
    mutate(distToCoaccGene = abs((start+end)/2 - geneStart)) %>%
    .$distToCoaccGene
other.ranges$peakID <- GRToString(other.ranges)
other.ranges$promoterID <- GRToString(promoter.ranges)

cor.vals <- seq(0, 1, 0.1)
links.per.promoter <- rep(0, length(cor.vals))
links.per.distal <- rep(0, length(cor.vals))

for(i in 1:length(cor.vals)){
    curr <- other.ranges[other.ranges$correlation >= cor.vals[i],] %>% as_tibble()
    links.per.promoter[i] <- median( curr %>% group_by(promoterID) %>% summarise(nLinks = n()) %>% .$nLinks )
    links.per.distal[i] <- median( curr %>% group_by(peakID) %>% summarise(nLinks = n()) %>% .$nLinks )
}

corr.df <- data.frame(links.per = c(links.per.promoter, links.per.distal),
                      corr.cutoff = c(cor.vals, cor.vals),
                      type = c(rep("Distal sites per promoter", length(links.per.promoter)),
                               rep("Promoter sites per distal", length(links.per.distal))))

pdf('manuscript_figures/Coaccess_Median_sites_per_corr_cutoff.pdf', width=10, height=6)
ggplot(corr.df, aes(x =corr.cutoff, y = links.per, color = type)) + geom_point() + geom_line() + ggClean() + xlab('Correlation cutoff') + ylab('Median')
dev.off()

other.ranges <- other.ranges[other.ranges$correlation > 0.5,]
promoter.ranges <- promoter.ranges[promoter.ranges$correlation > 0.5,]

saveRDS(object = other.ranges, file = 'ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coaccess_enhancers_corr_0.5_maxDist_1Mb_hg38.gr.rds')
saveRDS(object = promoter.ranges, file = 'ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coaccess_promoters_corr_0.5_maxDist_1Mb_hg38.gr.rds')

# lift over to hg19 for GWAS stuff

other.ranges.clean <- GRanges(seqnames = seqnames(other.ranges), ranges = IRanges(start = start(other.ranges), end = end(other.ranges)))
other.ranges.clean$idx <- 1:length(other.ranges.clean)
other.ranges.clean$correlation <- other.ranges$correlation

promoter.ranges.clean <- GRanges(seqnames = seqnames(promoter.ranges), ranges = IRanges(start = start(promoter.ranges), end = end(promoter.ranges)))
promoter.ranges.clean$idx <- 1:length(promoter.ranges.clean)

other.ranges.clean$gene_name <- promoter.ranges$nearestGene

path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch <- import.chain(path)

other.ranges.clean.hg19 <- unlist(liftOver(other.ranges.clean, ch))
other.ranges.clean.hg19$enhancer_chr <- as.character(seqnames(other.ranges.clean.hg19))
other.ranges.clean.hg19$enhancer_start <- start(other.ranges.clean.hg19)
other.ranges.clean.hg19$enhancer_end <- end(other.ranges.clean.hg19)
other.ranges.clean.hg19 <- other.ranges.clean.hg19[!duplicated(other.ranges.clean.hg19$idx),]

promoter.ranges.clean.hg19 <- unlist(liftOver(promoter.ranges.clean, ch))
promoter.ranges.clean.hg19$promoter_chr <- as.character(seqnames(promoter.ranges.clean.hg19))
promoter.ranges.clean.hg19$promoter_start <- start(promoter.ranges.clean.hg19)
promoter.ranges.clean.hg19$promoter_end <- end(promoter.ranges.clean.hg19)
promoter.ranges.clean.hg19 <- promoter.ranges.clean.hg19[!duplicated(promoter.ranges.clean.hg19$idx),]

saveRDS(object = other.ranges.clean.hg19, file = 'ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coaccess_enhancers_hg19.gr.rds')
saveRDS(object = promoter.ranges.clean.hg19, file = 'ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coaccess_promoters_hg19.gr.rds')



