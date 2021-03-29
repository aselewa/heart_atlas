library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

#gtf.file <- '/project2/gca/software/annotations/gencode.v19.annotation.gtf.gz'
gtf.file <- '/project2/xinhe/shared_data/gencode/gencode.v35.GRCh38.gtf.gz'

my.gtf <- rtracklayer::import(con = gtf.file, format = 'gtf')
seqlevels(my.gtf, pruning.mode = "coarse") <- paste0("chr",1:22)

my.gtf.protein <- my.gtf[which(my.gtf$gene_type=="protein_coding"),]

my.genes <- my.gtf.protein[my.gtf.protein$type=="gene",]

canonical.transcripts <- my.gtf.protein %>% 
    as_tibble() %>% 
    dplyr::filter(type == "transcript") %>% 
    group_by(gene_id) %>% 
    mutate(transLength = abs(end - start)) %>% 
    arrange(-transLength) %>% 
    dplyr::slice(1)

canoncial.transcripts.str <- canonical.transcripts %>% .$transcript_id

my.gtf.protein.canonical <- my.gtf.protein[my.gtf.protein$transcript_id %in% canoncial.transcripts.str,] # keep only canonical transcripts
my.exons <- my.gtf.protein.canonical[my.gtf.protein.canonical$type=="exon",]
my.UTR <- my.gtf.protein.canonical[my.gtf.protein.canonical$type=="UTR",]
my.introns <- GenomicRanges::setdiff(my.genes, my.exons)


my.genes.plus <- my.genes[strand(my.genes)=="+",]
my.genes.neg <- my.genes[strand(my.genes)=="-",]
my.promoters.plus <- GRanges(seqnames(my.genes.plus), ranges = IRanges(start = start(my.genes.plus) - 2000, end = start(my.genes.plus)), strand = strand(my.genes.plus))
my.promoters.neg <- GRanges(seqnames(my.genes.neg), ranges = IRanges(start = end(my.genes.neg), end = end(my.genes.neg) + 2000), strand = strand(my.genes.neg))
my.promoters <- append(my.promoters.plus, my.promoters.neg)
my.promoters$gene_name <- c(my.genes.plus$gene_name, my.genes.neg$gene_name)

exon.chr <- c(as.character(seqnames(my.exons)), as.character(seqnames(my.exons)))
exon.pos <- c(start(my.exons), end(my.exons))
exon.gene.name <- c(my.exons$gene_name, my.exons$gene_name)
my.splice.junc <- GRanges(seqnames = exon.chr, ranges = IRanges(start = exon.pos-100, end = exon.pos+100), gene_name = exon.gene.name)

annots <- list(
    genes = my.genes,
    exons = my.exons,
    introns = my.introns,
    UTRs = my.UTR,
    promoters = my.promoters,
    splice_junctions = my.splice.junc
)

saveRDS(annots, 'hg38_gtf_genomic_annots.gr.rds')
