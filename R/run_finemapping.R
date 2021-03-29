library(finemappeR)

args <- commandArgs(trailingOnly = T)
EUR_LD_1KG <- '/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds'
bigSNP <- bigsnpr::snp_attach(EUR_LD_1KG)

sumstats <- readRDS(args[1])
z <- unique(sumstats$locus)

susie.res <- run.susie(sumstats = sumstats, bigSNP = bigSNP, ldchunk = z, L = 1, prior = F)
susie.list <- list()
susie.list[[as.character(z)]] <- susie.res

sumstats.fm <- merge_susie_sumstats(susie.list, sumstats)

saveRDS(sumstats.fm, paste0(args[2],'_',z,'.df.rds'))