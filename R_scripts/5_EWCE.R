library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(DESeq2)

# follow author's pipeline for gene selection
#---------------
counts <- read.csv(gzfile('GSE205161_20220525-geo_raw_counts.csv.gz'))
rownames(counts) <- counts$ensemble_id_version
counts$ensemble_id_version <- NULL
# convert to gene symbols
gsymb <- mapIds(EnsDb.Hsapiens.v86, sub('\\..*$', '', 
                                        rownames(counts)), 'GENENAME', 'GENEID')
# remove not converted genes and assign new names (symbols)
good_conversion <- !is.na(gsymb) & !duplicated(gsymb)
counts <- counts[good_conversion,]
rownames(counts) <- gsymb[good_conversion]
# get gene description
gene_descr <- mapIds(org.Hs.eg.db, rownames(counts), 'GENENAME', 'SYMBOL')
# get chromosomes
gene_chrom <- mapIds(EnsDb.Hsapiens.v86, rownames(counts), 'SEQNAME', 'SYMBOL')
#-----------
mtdata <- read.csv('doi_10.5061_dryad.np5hqbzx2__v6/NTM_clinical.csv')
rownames(mtdata) <- mtdata$CFB_study_id
mtdata$CFB_study_id <- NULL
mtdata <- mtdata[colnames(counts),]
numeric_cols <- c('age_years','age_sample','height_cms', 
                  'base_mean_FEV1', 'base_FEV1', 'bmi_baseline')
date_cols <- grep('date', colnames(mtdata), value = T)
for(colname in setdiff(colnames(mtdata), c(numeric_cols, date_cols)))
  mtdata[[colname]] <- factor(mtdata[[colname]])
# filter genes
keep_hemo <- !grepl('hemoglobin subunit', gene_descr)
keep_sex <- !grepl('^X|Y', gene_chrom)
keep_counts <- rowSums(counts) < 7000000
table(keep <- keep_hemo & keep_sex & keep_counts)
counts <- counts[keep, ]
#------------
idx <- ! (rownames(mtdata) %in% c('CFB2006'))
counts <- counts[,idx]
mtdata <- mtdata[idx,]

#-- all
dds <- DESeqDataSetFromMatrix(counts, colData = mtdata, design = ~ ntm_disease)
dds <- DESeq(dds)
# vsd <- vst(dds)
#-- MAC
subset_idx <- mtdata$mycobacteria_species == 'MAC'
counts_MAC <- counts[,subset_idx]
mtdata_MAC <- mtdata[subset_idx,]
dds_MAC <- DESeqDataSetFromMatrix(counts_MAC, colData = mtdata_MAC, design = ~ ntm_disease)
dds_MAC <- DESeq(dds_MAC)
# vsd_MAC <- vst(dds_MAC)
#-- MABs
subset_idx <- mtdata$mycobacteria_species == 'MABs'
counts_MABs <- counts[,subset_idx]
mtdata_MABs <- mtdata[subset_idx,]
dds_MABs <- DESeqDataSetFromMatrix(counts_MABs, colData = mtdata_MABs, design = ~ ntm_disease)
dds_MABs <- DESeq(dds_MABs)
# vsd_MABs <- vst(dds_MABs)

res <- results(dds)
res_MAC <- results(dds_MAC)
res_MABs <- results(dds_MABs)

save(res, res_MAC, res_MABs, file = 'GSE205161_res.rda')
save(dds, dds_MAC, dds_MABs, file = 'GSE205161_dds.rda')
# 
library(EWCE)
library(Seurat)

scdata <- readRDS('scdata.RDS')

lung_mrna <- list()
lung_mrna$exp <- scdata@assays$RNA@counts
lung_mrna$exp_scT_normed <- EWCE::sct_normalize(lung_mrna$exp)

saveRDS(lung_mrna, file = 'ewce_lung_mrna.rds')

#-----------

library(EWCE)
lung_mrna <- readRDS(file = 'ewce_lung_mrna.rds')
scdata_meta <-
  readRDS('scdata_meta.RDS')
scdata_meta <- scdata_meta[colnames(lung_mrna$exp_scT_normed),]

lung_mrna$annot <- list()
lung_mrna$annot$level1class <- lung_mrna$annot$level2class <- scdata_meta$celltype

exp_Lung_DROPPED <- EWCE::drop_uninformative_genes(
  exp = lung_mrna$exp_scT_normed, 
  input_species = "human",
  output_species = "human",
  level2annot = lung_mrna$annot$level2class) 

annotLevels <- list(level1class=lung_mrna$annot$level1class,
                    level2class=lung_mrna$annot$level2class)

fNames_LungOnly <- EWCE::generate_celltype_data(
  exp = exp_Lung_DROPPED,
  annotLevels = annotLevels,
  groupName = "ewce1_6_Human_PBMCs", 
  savePath = '')

#----------

library(EWCE)
library(ggplot2)

load('ctd_ewce1_6_Human_PBMCs.rda')
load(#res, res_MAC, res_MABs, 
  file = 'GSE205161_res.rda')

ewceRes <- function(res, ntop, reps, seed=1) {
  set.seed(seed)
  hitsDN <- head(rownames(res[order(res$stat),]), ntop)
  hitsUP <- head(rownames(res[order(res$stat, decreasing = T),]), ntop)
  #------ enrich
  mouse.hits <- c(hitsUP, hitsDN)

  full_resultsUP = bootstrap_enrichment_test(sct_data=ctd,hits=hitsUP,#bg=mouse.bg,
                                             reps=reps,annotLevel=1, sctSpecies='human', 
                                             genelistSpecies='human')
  full_resultsDN = bootstrap_enrichment_test(sct_data=ctd,hits=hitsDN,#bg=mouse.bg,
                                             reps=reps,annotLevel=1, sctSpecies='human',
                                             genelistSpecies='human')
  return(list(up=full_resultsUP, dn=full_resultsDN, hitsDN=hitsDN, hitsUP=hitsUP))
}
ntop <- 250 # 250 per il paper (con 1000 vengono fuori le cd16 per MAC)
reps <- 10000
ewce <- ewceRes(res, ntop, reps)
ewce_MAC <- ewceRes(res_MAC, ntop, reps)
ewce_MABs <- ewceRes(res_MABs, ntop, reps)

# UP-reg 

ewce_UP_res <- rbind(
  # cbind(ewce$up$results, mycobacteria='all'),
  cbind(ewce_MAC$up$results, mycobacteria='MAC'),
  cbind(ewce_MABs$up$results, mycobacteria='MABSC')
)
ewce_UP_res$p[ewce_UP_res$p==0] <- 1/reps/10
ewce_UP_res$log10.p <- log10(ewce_UP_res$p)
ewce_UP_res$sd_from_mean[ewce_UP_res$sd_from_mean<0] <- 0
# (ewce_up <- ggplot(ewce_UP_res, aes(mycobacteria, CellType, size=sd_from_mean, color=log10.p)) + 
#     geom_point() + scale_color_gradient(high = 'white', low = 'red')+ theme_bw() + 
#     scale_size_area()) + ggtitle('EWCE for UP-regulated genes')

# DOWN-reg

ewce_DN_res <- rbind(
  # cbind(ewce$dn$results, mycobacteria='all'),
  cbind(ewce_MAC$dn$results, mycobacteria='MAC'),
  cbind(ewce_MABs$dn$results, mycobacteria='MABSC')
)
ewce_DN_res$p[ewce_DN_res$p==0] <- 1/reps/10
ewce_DN_res$log10.p <- log10(ewce_DN_res$p)
ewce_DN_res$sd_from_mean[ewce_DN_res$sd_from_mean<0] <- 0
# (ewce_dn <- ggplot(ewce_DN_res, aes(mycobacteria, CellType, size=sd_from_mean, color=log10.p)) + 
#     geom_point() + scale_color_gradient(high = 'white', low = 'blue') + theme_bw() + scale_size_area()) + 
#   ggtitle('EWCE for DOWN-regulated genes')

# UP + DOWN

ewce_UP_rev <- ewce_UP_res
ewce_UP_rev$log10.p <- -ewce_UP_rev$log10.p
ewce_res_all <- rbind(
  cbind(ewce_DN_res, sign='down'), 
  cbind(ewce_UP_rev, sign='up')
)
ewce_res_all$sign <- factor(ewce_res_all$sign, 
                            levels=c('up','down'))
ewce_res_all$CellType <- factor(ewce_res_all$CellType,
                                levels = sort(unique(ewce_res_all$CellType), 
                                              decreasing = T))
ewce_res_all$mycobacteria <- factor(ewce_res_all$mycobacteria)

levels(ewce_res_all$mycobacteria) <- 
  c('MABSC-PD\nvs\nMABSC\n(5 vs 6)','MAC-PD\nvs\n7 MAC\n(7 vs 15)')

(ewce_all <- ggplot(ewce_res_all, aes(mycobacteria, CellType, size=sd_from_mean, color=log10.p)) + facet_wrap( ~ sign ) + 
    geom_point() + scale_color_gradient2(high = 'red', mid = 'white', low = 'blue') + theme_bw() + scale_size_area() + 
  xlab('') + ylab(''))
#ggtitle('EWCE for DOWN-regulated genes')

ggsave('ewce.pdf', ewce_all, 
       height = 4, width = 6) 
