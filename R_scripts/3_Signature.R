
# Import Libraries --------------------------------------------------------



library(Seurat)
library(patchwork)
library(RColorBrewer)
library(Matrix)
library(readxl)
library(pheatmap)
library(fgsea)
library(ggpubr)
library(enrichR)
library(rstatix)
library(ggvenn)
library(monocle3)
library(future)
library(SeuratWrappers)
library(clusterProfiler)
library(enrichplot)  
library(org.Hs.eg.db)
library(topGO)
library(EnhancedVolcano)
library(hdWGCNA)
library(tidyverse)





# Define signature --------------------------------------------------------

## Import DEGS
read_delim('write_R/scDeg.csv',delim=';')->Deg_list


## Extract Signature for Cluster 1

Deg_list%>%
  filter(Comparison%in%c('MABSC-PD,vs,MABSC','MABSC-PD,vs,CF')) %>% 
  filter(Cluster==1) %>% 
  filter(!str_detect(Genes,'^MT|^RP|^AC|^X|^Y|^AL')) %>%
  mutate(fold=avg_log2FC>0) %>% 
  filter(p_val_adj<0.05,abs(avg_log2FC)>0.25) %>% 
  dplyr::group_by(Genes,fold) %>% 
  filter(n()==2,fold) %>% 
  ungroup() %>% dplyr::select(Genes) %>% distinct() %>% pull(Genes)->sig_gene1

saveRDS(sig_gene1, file = 'sig_gene1.rds')


# Subset Cluster 1 --------------------------------------------------------


## Generate subset of the cluster 1

scdata<-readRDS('write_R/scdata.rds')

scdata_s<-subset(scdata,seurat_clusters==1)

DefaultAssay(scdata_s)<-'integrated'

## Process the subset
FindVariableFeatures(scdata_s)->scdata_s
scdata_s<-ScaleData(scdata_s,verbose=T)

scdata_s <- Seurat::RunPCA(scdata_s, npcs = 30, verbose = FALSE,)
scdata_s <- RunUMAP(scdata_s, reduction = "pca", dims = 1:30)
scdata_s <- FindNeighbors(scdata_s, reduction = "pca", dims = 1:30)
scdata_s <- FindClusters(scdata_s, resolution = 0.1)

## Optional: visualize

wrap_plots(DimPlot(scdata_s,label=T,label.box = T)+ggtitle(''),
           DimPlot(scdata_s,group.by = 'Groups_f',
                   cols=rev(c('#b11600','#e69f11','green4')))+
             ggtitle(''))

ggsave('image/sub_UMAP_CL1.png',height=5,width=7)


## Create signature module


DefaultAssay(scdata_s)<-'RNA'


scdata_s<-AddModuleScore(scdata_s,features = list(c(sig_gene1)),name='c1_sig')


## Visualize signature at group level in cluster 1

lapply(names(table(scdata_s$Groups_f)),function(x){
  FeaturePlot(subset(scdata_s,Group_f==x),
              #split.by = 'sample_f',
              #title=x,
              cols = viridis::viridis(n=1000,option='turbo'),
              #keep.scale = NULL,
              features='c1_sig1')+ggtitle('')})->plot_list

wrap_plot(plot_list,ncol=3)

ggsave('image/Signature_Cluster1.png',dpi=300)



## Visualize signature at group level in whole sample


DefaultAssay(scdata)<-'RNA'


scdata<-AddModuleScore(scdata,features = list(c(sig_gene1)),name='c1_sig')

lapply(names(table(scdata$Groups_f)),function(x){
  FeaturePlot(subset(scdata,Groups_f==x),
              #split.by = 'sample_f',
              #title=x,
              cols = viridis::viridis(n=1000,option='turbo'),
              #keep.scale = NULL,
              features='c1_sig1')+ggtitle('')})->plot_list

wrap_plot(plot_list,ncol=3)

ggsave('image/Signature_WholeMap.png',dpi=300)


