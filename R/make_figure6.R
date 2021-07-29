setwd('/project2/gca/aselewa/heart_atlas_project/')
source('../R/analysis_utils.R')
library(ggpubr)
library(rstatix)

finemap.res <- suppressMessages(readr::read_tsv('eQTL_enrich/gtex_finemapping/GTEx_v8_finemapping_DAPG/Heart_LV_Finemapping_CS95.txt', col_names = F))
colnames(finemap.res) <- c("tissue","gene","cluster_id","cluster_pip","variant_id","variant_pip")


high_pip_eqtls <- finemap.res[finemap.res$variant_pip>0.8,] %>% 
  group_by(variant_id) %>% 
  arrange(-variant_pip) %>% 
  slice(1) %>% 
  select(variant_id, variant_pip, gene)

eqtl.gr <- snpIDtoGR(high_pip_eqtls$variant_id)
eqtl.gr$SNP <- high_pip_eqtls$variant_id
eqtl.gr$gene <- high_pip_eqtls$gene
eqtl.gr$finemap.pip <- high_pip_eqtls$variant_pip


# load matched LD SNPs

snp_match <- vroom::vroom('matched_SNPs/eQTL_top_pip_hg19_5batches/SNPmatch_10_all.txt', col_names = F)
snp_match$X1 <- NULL
snp_match$X2 <- NULL

random.snp.gr.list <- list()
path <- '/project2/gca/software/liftOver_Chains/hg19ToHg38.over.chain'
ch <- import.chain(path)

for(i in 1:ncol(snp_match)){
  curr <- snp_match %>% pull(i)
  chroms <- sapply(strsplit(curr, split = ':'), function(x){x[1]})
  pos <- as.integer(sapply(strsplit(curr, split = ':'), function(x){x[2]}))
  curr.gr <- GRanges(seqnames = chroms, ranges = IRanges(start = pos), SNP = curr)
  seqlevelsStyle(curr.gr) <- "UCSC"
  random.snp.gr.list[[i]] <- unlist(liftOver(curr.gr, ch))
}

# Get tissue activity of each eQTL from https://zenodo.org/record/3727189 (GTEx v8 paper)
eqtl.lfsr <- suppressMessages(vroom::vroom('eQTL_enrich/broadinstitute-gtex-v8-a014b43/data/Fig6C_all_top.z_lfsr.sig.pruned.txt.gz'))

same.snps <- intersect(eqtl.gr$SNP, eqtl.lfsr$variant)
heart.eqtl.lfsr <- eqtl.lfsr[eqtl.lfsr$variant %in% same.snps,]
eqtl.gr <- eqtl.gr[eqtl.gr$SNP %in% same.snps]
saveRDS(eqtl.gr, 'eQTL_enrich/gtex_finemapping/HighPIP_eQTLs_w_LFSR.gr.rds')

#add tissue activity
eqtl.gr$ntissues.active <- heart.eqtl.ntissues$mean_tissues[match(eqtl.gr$SNP, heart.eqtl.ntissues$eqtl)]


#################### 6A -- tissue sharing #################### 

# all eqtls
tissues.active <- rowSums(eqtl.lfsr[,3:ncol(eqtl.lfsr)] < 0.01, na.rm = T) 
all.eqtl.ntissues <- data.frame(eqtl=eqtl.lfsr$variant, ntissues=tissues.active) %>% group_by(eqtl) %>% summarise(mean_tissues = mean(ntissues))
all.eqtl.ntissues <- all.eqtl.ntissues[all.eqtl.ntissues$mean_tissues > 0,]

# heart eqtls
tissues.active <- rowSums(heart.eqtl.lfsr[,3:ncol(heart.eqtl.lfsr)] < 0.01, na.rm = T) 
heart.eqtl.ntissues <- data.frame(eqtl=heart.eqtl.lfsr$variant, ntissues=tissues.active) %>% group_by(eqtl) %>% summarise(mean_tissues = mean(ntissues))
heart.eqtl.ntissues <- heart.eqtl.ntissues[heart.eqtl.ntissues$mean_tissues > 0,]


betterHist <- function(X){
  bks <- seq(0, 50, 5)
  labs <- sapply(1:(length(bks)-1), function(x){paste0(bks[x]+1,'-',bks[x+1])})
  bin.count <- table(cut(X$mean_tissues, breaks = bks, labels = labs))
  bin.count <- bin.count/sum(bin.count)
  bin.count.df <- data.frame(count=as.numeric(bin.count), breaks=factor(names(bin.count), levels = labs))
}

all.eqtl.bin <- betterHist(all.eqtl.ntissues)
heart.eqtl.bin <- betterHist(heart.eqtl.ntissues)

result <- rbind(all.eqtl.bin, heart.eqtl.bin)
result$eQTLs <- c(rep(c("All","Heart"), each=nrow(all.eqtl.bin)))

ggplot(result, aes(x=breaks, y=count, fill=eQTLs)) + 
  geom_bar(stat="identity", position = "dodge", width=0.8) + 
  ggClean() + xlab('Tissues with LFSR < 0.01') + 
  ylab('Proportion of eQTLs') + 
  ggtitle('Tissue Sharing of Heart eQTLs') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_fill_manual(values = c("#D6604D","#4393C3"))


#################### 6B eQTL overlap with OCRs #################### 

## prepare annotations
# remove exons
satac <- suppressMessages(ArchR::loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/', showLogo = F))
union.set <- ArchR::getPeakSet(satac)
union.set$peakID <- GRToString(union.set)
union.exons <- union.set[union.set$peakType == "Exonic", ]
union.nonExon <- union.set[union.set$peakType != "Exonic",]

annots <- readRDS('genomic_annotations/hg38_gtf_genomic_annots.gr.rds')
exon.annots <- annots$exons
utr.annots <- annots$UTRs
exon.annots <- removeOverlaps(X = exon.annots, to.remove = utr.annots)
intron.annots <- annots$introns
intron.nonpeak <- removeOverlaps(intron.annots, to.remove = union.set)

# get cell-type specific markers and shared category
markers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')
markers <- lapply(markers, function(x){removeOverlaps(X = x, to.remove = union.exons)})

marker.ids <- GRToString(unlist(GRangesList(markers)))
marker.ids.count <- table(marker.ids)
shared.markers <- names(marker.ids.count)[marker.ids.count > 1 ]
ct.specific.markers <- names(marker.ids.count)[marker.ids.count == 1]

markers <- lapply(markers, function(x){x[GRToString(x) %in% ct.specific.markers,]})
union.set.shared <- union.nonExon[!(union.nonExon$peakID %in% ct.specific.markers),]

# cluster the shared set into multiple cell types
mean.acc.obj <- getGroupSE(satac, useMatrix = 'PeakMatrix', groupBy = 'CellTypes', divideN = T, scaleTo = 10^5)
mean.acc.mat <- mean.acc.obj@assays@data$PeakMatrix
rdata <- rowData(mean.acc.obj)
rownames(mean.acc.mat) <- paste0(rdata$seqnames,':',rdata$start,'-',rdata$end)
mean.acc.mat <- mean.acc.mat[union.set.shared$peakID,]

for(i in 1:ncol(mean.acc.mat)){
  curr.thresh <- as.numeric(quantile(mean.acc.mat[,i],0.75))
  mean.acc.mat[,i] <- 1*(mean.acc.mat[,i] > curr.thresh)
}

peak.count <- rowSums(mean.acc.mat)

# do the overlapping
markers$Neuronal <- NULL
markers$`Smooth Muscle` <- NULL
peaks.list <- markers
 
peaks.list[["2-3"]] <- StringToGR(rownames(mean.acc.mat)[peak.count >= 2 & peak.count <= 3])
peaks.list[["4+"]] <- StringToGR(rownames(mean.acc.mat)[peak.count >= 4])

peaks.list[["Exons"]] <- exon.annots
peaks.list[["UTR"]] <- utr.annots
peaks.list[["Intron"]] <- intron.nonpeak

peaks.eqtls <- join_overlap_list(gr.list = peaks.list, X = eqtl.gr)
peaks.random <- lapply(random.snp.gr.list, function(x){
  join_overlap_list(gr.list = peaks.list, X = x)
})

random.dist.df.list <- list()
for(i in 1:length(peaks.random)){
  curr <- peaks.random[[i]]
  snpsIn <- unique(unlist(sapply(curr, function(x){x$SNP})))
  curr$Unassigned <- random.snp.gr.list[[i]][!(random.snp.gr.list[[i]]$SNP %in% snpsIn),]
  
  curr.df <- as.data.frame(sapply(curr, FUN = function(x){length(x)/length(random.snp.gr.list[[i]])}))
  colnames(curr.df) <- c("freq")
  random.dist.df.list[[i]] <- curr.df
}

random.dist.avg <- data.frame(freq = rowMeans(Reduce(cbind, random.dist.df.list)))
random.dist.avg$category <- rownames(random.dist.avg)

snpsIn <- unique(
  unlist(
    sapply(peaks.eqtls, function(x){x$SNP})
  )
)
peaks.eqtls$Unassigned <- eqtl.gr[!(eqtl.gr$SNP %in% snpsIn),]

## do the plotting

plot.levels <- c(names(peaks.list),"Unassigned")
peak.dist.df <- as.data.frame(sapply(peaks.eqtls, FUN = function(x){length(x)/length(eqtl.gr)}))
colnames(peak.dist.df) <- c("freq")
peak.dist.df$category <- rownames(peak.dist.df)

peak.set.dist.df <- Reduce(rbind, list(random.dist.avg, peak.dist.df))
peak.set.dist.df$SNPs <- c(rep("Random SNPs", nrow(random.dist.avg)),
                           rep("eQTLs", nrow(peak.dist.df)))

plot.levels <- c(names(peaks.gr.list),"Unassigned")
peak.set.dist.df$category <- factor(peak.set.dist.df$category, 
                                    levels = plot.levels)

n_cat <- nrow(peak.dist.df)
enrich <- peak.dist.df$freq / random.dist.avg$freq
enrich.df <- data.frame(category= peak.dist.df$category,
                        group1 = rep("eQTLs", n_cat),
                        group2 = rep("Random SNPs", n_cat),
                        p = round(enrich,2),
                        y.position = sapply(1:n_cat, function(x)max(peak.dist.df$freq[x], random.dist.avg$freq[x]))) %>% 
  as_tibble() %>% 
  mutate(group1 = unfactor(group1), group2 = unfactor(group2))

enrich.df$category <- factor(enrich.df$category, levels = plot.levels)
enrich.df <- enrich.df %>% add_x_position(x = "category", dodge = 0.75)
enrich.df$p[enrich.df$p == 0] <- "ns"

pdf('manuscript_figures/eQTL_overlap_with_OCRs.pdf', width=8, height=6)
ggbarplot(peak.set.dist.df, x = "category", y = "freq", fill = "SNPs", position = position_dodge(.75), palette = "Paired") + 
  stat_pvalue_manual(enrich.df, tip.length = 0.01, step.increase = 0, bracket.nudge.y = 0.01) + coord_cartesian(ylim = c(0, 0.4)) +
  ggClean(rotate_axis = T)
dev.off()

#################### 6B Tissue sharing of each category #################### 

peaks.tissue.sharing <- sapply(peaks.eqtls, function(x){x$ntissues.active})
peaks.tissue.sharing <- setNames(unlist(peaks.tissue.sharing, use.names=F),rep(names(peaks.tissue.sharing), lengths(peaks.tissue.sharing)))
peaks.tissue.sharing.df <- data.frame(category=names(peaks.tissue.sharing), ntissues=peaks.tissue.sharing)
peaks.tissue.sharing.df$category <- factor(peaks.tissue.sharing.df$category, levels=plot.levels)
peaks.tissue.sharing.df <- peaks.tissue.sharing.df[!is.na(peaks.tissue.sharing.df$ntissues),]

peaks.tissue.sharing.df$category <- factor(peaks.tissue.sharing.df$category, levels = plot.levels)

pdf('manuscript_figures/eQTL_Tissue_Sharing_Violin.pdf', width=20, height=3)
ggplot(peaks.tissue.sharing.df, aes(x=category, y=ntissues)) + 
  geom_violin(fill="lightblue") + 
  geom_boxplot(width=0.1) + 
  ggClean() + LegendOff() + ylab("Tissues Active (LFSR < 0.01)") + xlab("") +
  scale_y_continuous(breaks=seq(0,50,10))
dev.off()

#################### 6C Comparison in other tissues #################### 

dapg.eqtl <- readr::read_tsv(file = 'eQTL_enrich/gtex_finemapping/GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.txt.gz')

peaks.eqtls.atac <- list(peaks.eqtls$Cardiomyocyte, 
                         peaks.eqtls$Endothelial,
                         peaks.eqtls$Fibroblast,
                         c(peaks.eqtls$Lymphoid, peaks.eqtls$Myeloid))

names(peaks.eqtls.atac) <- c("Cardiomyocyte","Endothelial","Fibroblast","Immune")
tissues.to.keep <- c("Whole_Blood", "Lung","Heart_Atrial_Appendage","Testis","Brain_Hippocampus")
#tissues.to.keep <- unique(dapg.eqtl$tissue_id)
dapg.eqtl.filt <- dapg.eqtl %>% 
  filter(variant_pip > 0.8) %>%
  filter(tissue_id %in% tissues.to.keep) %>% 
  mutate(eqtl_identifier = paste0(variant_id, '__', gene_id))

res <- lapply(peaks.eqtls.atac, function(x){
  eqtlsIn <- x %>% as_tibble() %>% mutate(eqtl_identifier = paste0(SNP, '__', gene)) %>% .$eqtl_identifier;
  dapg.eqtl.filt %>% group_by(tissue_id) %>% summarise(non_cm_overlap = sum(eqtlsIn %in% eqtl_identifier)/length(eqtlsIn))
})
res.mat <- Reduce(rbind, res)
res.mat$type <- rep(names(res), each=length(tissues.to.keep))
res.mat.cm <- res.mat[res.mat$type=="Cardiomyocyte",]
res.mat$tissue_id <- factor(res.mat$tissue_id, levels  = c("Heart_Atrial_Appendage", "Lung", "Whole_Blood","Testis","Brain_Hippocampus"))

pdf('manuscript_figures/eQTL_overlap_other_tissues.pdf', width=8, height=6)
ggplot(res.mat, aes(x= tissue_id, y = non_cm_overlap, fill=type)) + geom_bar(stat='identity', position='dodge') + ggClean(rotate_axis = F) +
  xlab('Proportion of eQTLs') + ylab('Tissue') +   
  scale_fill_manual(values = c("#b22222","#8DD3C7","#BEBADA", "#FB8072")) + scale_y_continuous(breaks = seq(0, 0.55, 0.1)) +
  coord_flip()
dev.off()

#################### 6D power simulation #################### 



