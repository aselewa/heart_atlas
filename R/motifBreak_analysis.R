library(motifbreakR)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(MotifDb)
library(tidyverse)

setwd('/project2/gca/aselewa/heart_atlas_project/')

finemapped.res <- readr::read_csv('GWAS/aFib_Finemapped_minPIP_0.2_03252021.csv')
#finemapped.res <- finemapped.res %>% filter(`Chromatin status` == 'CM Specific ATAC' | `Chromatin status` == 'CM Shared ATAC')

variants <- snps.from.rsid(rsid = finemapped.res$SNP,
                           dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37::SNPlocs.Hsapiens.dbSNP142.GRCh37,
                           search.genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)

drops <- c("rs3731746:T","rs60632610:A","rs60632610:G")
variants <- variants[!(names(variants) %in% drops),]

chroms <- sapply(strsplit(curr, split = ':'), function(x){x[1]})
pos <- as.integer(sapply(strsplit(curr, split = ':'), function(x){x[2]}))
finemapped.res <- finemapped.res %>% left_join(., variants %>% as_tibble %>% dplyr::select(SNP_id, REF, ALT) %>% rename(SNP = SNP_id))

#motifs.of.interest <- c("TBX5", "GATA4", "NKX2-5", "PITX2", "MEF2A","MEF2C", "TEAD1", "PRRX1", "HAND2", "ETV1", "NR2F2")
motif.cor.res <- read_tsv('manuscript_figures/figure2/All_TFMotif_Expression_Correlations_CisBP.tsv')
cardiac.tfs <- motif.cor.res$gene_name[motif.cor.res$correlation > 0.4]
cardiac.tfs <- c(cardiac.tfs, "NKX2-5","PRRX1")

res <- query(MotifDb, andStrings = c("hsapiens"), orStrings = cardiac.tfs)
res <- res[!duplicated(mcols(res)$geneSymbol)]
res <- res[mcols(res)$geneSymbol %in% cardiac.tfs,]

hand2.pwm <- read.table('misc/HAND2_motif.txt', sep = '')
colnames(hand2.pwm) <- seq(1:10)
rownames(hand2.pwm) <- c("A","C","G","T")
hand2.pwm <- t(t(hand2.pwm)/colSums(hand2.pwm))
res$HAND2 <- hand2.pwm
mcols(res)["HAND2",] <- "HAND2"

motif.breaks <- motifbreakR(snpList = variants, 
                            show.neutral = T,
                            filterp = TRUE,
                            pwmList = res, 
                            threshold = 1e-4,
                            method = "ic",
                            bkg = c(A=0.25, C=0.25, G=0.25, T=0.25), 
                            BPPARAM = BiocParallel::MulticoreParam(workers = 16))

saveRDS(motif.breaks, 'misc/cardiac_tf_motifbreakR_analysis_CM.gr.rds')

motif.breaks.tbl <- motif.breaks %>% as_tibble() %>% distinct(SNP_id, geneSymbol, .keep_all = T) %>% dplyr::select(SNP_id, geneSymbol, effect) %>%
    mutate(tf.info = paste0(geneSymbol," - ", effect)) %>% dplyr::select(-geneSymbol, -effect) %>% dplyr::rename(SNP = SNP_id) %>%
    group_by(SNP) %>% summarise(tf.info = paste0(tf.info, collapse=','))

finemapped.res <- finemapped.res %>% left_join(., motif.breaks.tbl)
finemapped.res[is.na(finemapped.res)] <- ""

finemapped.res %>% readr::write_csv('GWAS/aFib_Finemapped_minPIP_0.2_04062021_CardiacTFs.csv')


res <- query(MotifDb, andStrings = c("hsapiens"))
hand2.pwm <- read.table('misc/HAND2_motif.txt', sep = '')
colnames(hand2.pwm) <- seq(1:10)
rownames(hand2.pwm) <- c("A","C","G","T")
hand2.pwm <- t(t(hand2.pwm)/colSums(hand2.pwm))
res$HAND2 <- hand2.pwm
mcols(res)["HAND2",] <- "HAND2"
res <- res[!duplicated(mcols(res)$geneSymbol)]
res <- res[-contains('::', ignore.case = T, mcols(res)$geneSymbol)]
res <- res[!is.na(mcols(res)$geneSymbol)]

motif.breaks.all <- motifbreakR(snpList = variants, 
                            filterp = TRUE,
                            pwmList = res, 
                            threshold = 1e-4,
                            method = "ic",
                            bkg = c(A=0.25, C=0.25, G=0.25, T=0.25), 
                            BPPARAM = BiocParallel::MulticoreParam(workers = 16))

saveRDS(motif.breaks.all, 'misc/motifbreakR_analysis_all_TFs.gr.rds')

motif.breaks.tbl <- motif.breaks.all %>% as_tibble() %>% 
    rename(SNP = SNP_id) %>% 
    mutate(snp.score = abs(scoreRef-scoreAlt)) %>%
    group_by(SNP) %>% arrange(-snp.score) %>% dplyr::slice(1:5) %>%
    summarise(top5_TFs = paste0(geneSymbol, collapse=',')) %>% select(SNP, top5_TFs)

finemapped.res <- finemapped.res %>% left_join(., motif.breaks.tbl)
finemapped.res %>% readr::write_csv('GWAS/aFib_Finemapped_minPIP_0.2_04062021_TFs.csv')


for(i in 1:nrow(res)){
    if((res$a0[i] == res$ALT[i]) & (res$a1[i] == res$REF[i])){
        tt <- res$a1[i]
        res$a1[i] <- res$a0[i]
        res$a0[i] <- tt
        res$beta[i] <- -res$beta[i]
    }
}
