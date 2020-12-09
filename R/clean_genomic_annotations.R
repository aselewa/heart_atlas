# load protein coding annotations
genomic.annots <- rtracklayer::import('/project2/xinhe/alan/hg38_annotations/gencode.v29.annotation.gtf.gz', format = 'gtf')
genomic.annots <- genomic.annots[genomic.annots$source=="HAVANA",]
genomic.annots <- genomic.annots[genomic.annots$gene_type=="protein_coding",]

gene.annots <- genomic.annots[genomic.annots$type=="gene",]

# select canonical transcript per gene
exon.annots <- genomic.annots[genomic.annots$type=="exon",]
exon.annots$widths <- GenomicRanges::width(exon.annots)
gene.trans.size <- as.data.frame(exon.annots@elementMetadata) %>% group_by(gene_id, transcript_id) %>% summarise(size=sum(widths)) 
largest.trans.id <- gene.trans.size %>% group_by(gene_id) %>% arrange(-size) %>% slice(1) %>% pull(transcript_id)
genomic.annots <- genomic.annots[genomic.annots$transcript_id %in% largest.trans.id, ]

exon.annots <- genomic.annots[genomic.annots$type=="exon",]
utr.annots <- genomic.annots[genomic.annots$type=="UTR",]

# get intronic annotations
intron.annots <- GenomicRanges::setdiff(gene.annots, exon.annots)
intron.annots$type <- "intron"

# get promoter annotations
promoter.annots <- gene.annots
start(promoter.annots[strand(promoter.annots)=="+",]) <- start(promoter.annots[strand(promoter.annots)=="+",]) - 5000
end(promoter.annots[strand(promoter.annots)=="+",]) <- start(promoter.annots[strand(promoter.annots)=="+",]) + 1000
end(promoter.annots[strand(promoter.annots)=="-",]) <- end(promoter.annots[strand(promoter.annots)=="-",]) + 5000
start(promoter.annots[strand(promoter.annots)=="-",]) <- end(promoter.annots[strand(promoter.annots)=="-",]) - 1000
promoter.annots$type <- "promoter"

# get exons without UTR annotations
exon.annots <- GenomicRanges::setdiff(exon.annots, utr.annots)
exon.annots$type <- "exon"

disjoint.annots <- GenomicRanges::bindROWS(exon.annots, list(utr.annots, intron.annots, promoter.annots))
disjoint.annots <- disjoint.annots[,"type"]
saveRDS(disjoint.annots, file = '../eQTL_enrich/annotations/hg38_disjoint_annotations.gr.rds')
