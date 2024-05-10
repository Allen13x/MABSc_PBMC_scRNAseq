library(clusterProfiler)
library(DESeq2)

load(#res, res_MAC, res_MABs, 
  file = 'GSE205161_res.rda')
load(#dds, dds_MAC, dds_MABs, 
  file = 'GSE205161_dds.rda')
sig_gene1 <- readRDS(file = 'sig_gene1.rds')

Bioplanet_geneSets <- read.csv('Bioplanet_pathways.csv')
TERM2GENE <- Bioplanet_geneSets %>% dplyr::select(PATHWAY_NAME, GENE_SYMBOL)

res_MABs_filt <- res_MABs[apply(counts(dds_MABs)>10,1,all),]
gene_list <- res_MABs_filt$log2FoldChange
names(gene_list) <- rownames(res_MABs_filt)
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_res <- GSEA(gene_list, seed = 0,
                 TERM2GENE = TERM2GENE, 
                 eps = 0, pvalueCutoff=2)

gsea_sign <- GSEA(gene_list, seed = 0,
                  TERM2GENE = data.frame(TERM='signature', 
                                         GENE=sig_gene1), 
                  eps = 0, pvalueCutoff=2, minGSSize = 0)

