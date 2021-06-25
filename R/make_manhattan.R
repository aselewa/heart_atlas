setwd('/project2/gca/aselewa/heart_atlas_project/')

library(tidyverse)

# gwas <- readRDS('GWAS/summary_statistics/aFib/ebi-a-GCST006414_aFib.df.rds')
# gwas <- gwas[gwas$pval>1,]
# gwas$pval[gwas$pval > 50] <- 50
# gwas$pos_kb <- gwas$pos/1e3

finemap.genes <- readRDS('GWAS/finemapping/aFib_Finemapped_GeneMapped.tble.rds')
genomic.annots <- readRDS('genomic_annotations/hg19_gtf_genomic_annots.gr.rds')
gene.cords <- genomic.annots$genes
gene.cords$chr <- as.integer(sub("chr", "", as.character(seqnames(gene.cords))))
gene.cords$start <- start(gene.cords) / 1e3

gene.cord.df <- gene.cords@elementMetadata %>% as_tibble() %>% select(gene_name, gene_id, chr, start)
finemap.genes.pip <- finemap.genes %>% dplyr::select(gene_name, gene_pip) %>% distinct()
gene.pip.summary <- left_join(finemap.genes.pip, gene.cord.df, on="gene_name")
gene.pip.summary$is_highlight <- gene.pip.summary$gene_pip > 0.5

null.genes.gr <- gene.cords[!(gene.cords$gene_name %in% gene.pip.summary$gene_name),]
null.gene.pip.summary <- data.frame(gene_name = null.genes.gr$gene_name, 
                                    gene_pip = 0, 
                                    gene_id = null.genes.gr$gene_id, 
                                    chr = as.integer(seqnames(null.genes.gr)),
                                    start = start(null.genes.gr) / 1e3,
                                    is_highlight = F, stringsAsFactors = F) %>% as_tibble()

full.gene.pip.summary <- bind_rows(gene.pip.summary, null.gene.pip.summary)

don <- full.gene.pip.summary %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(start)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(full.gene.pip.summary, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate( BPcum=start+tot)

axisdf <- don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

pdf(file = 'manuscript_figures/Fig5D_gene_manhattan.pdf', width = 12, height=6)

ggplot(don, aes(x=BPcum, y=gene_pip)) +
  
  # Show all points
  ggrastr::geom_point_rast( aes(color=as.factor(chr)), size=2) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(0,1.25)) +     # remove space between plot area and x axis
  
  # Add highlighted points
  ggrastr::geom_point_rast(data=subset(don, is_highlight==T), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  ggrepel::geom_label_repel( data=subset(don, is_highlight==T), aes(label=gene_name), size=4, min.segment.length = 0) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    text = element_text(size = 18),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  xlab("Chromosome") + 
  ylab("PIP")

dev.off()


