
# Import Libraries --------------------------------------------------------
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(Matrix)
library(readxl)
library(pheatmap)
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




# Import datasets ---------------------------------------------------------


#### Change the input folder in the respective foder

INPUT='PATH_TO_COUNTMATRIXES/'


#### Create a vector for indexing the samples

sample<-c(1,3,4,5,6,9,12,13,14,15)

datasc<-list()


for (q in sample) {
  
##Import raw countmatrix
  i=paste('adata_',i,sep='')
  raw=read_delim(paste(INPUT,q,'.counts.tsv.gz',sep=''),delim='\t') %>% 
    spread(cell,count) %>% 
    column_to_rownames('gene') %>% 
    replace(is.na(.),0)
##CreateSeuratObject  
  datasc[[i]]<-CreateSeuratObject(raw,min.cells = 5,min.features = 200)
  
##ADD MT percentage  
  datasc[[i]][['percent_mt']]<- PercentageFeatureSet(datasc[[i]],pattern="^MT-")
  
## Filter Out cells with less than 200 feature or more than 8000 and mt percentage higher than 25%
  datasc[[i]]<-subset(datasc[[i]],subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent_mt < 25)
  
## Normalizazion  
  
  
  datasc[[i]] <- NormalizeData(datasc[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
## Find Variabel Features  
  
  datasc[[i]] <- FindVariableFeatures(datasc[[i]], selection.method = "vst", nfeatures = 2000)
## Scale Data
  all.genes <- rownames(datasc[[i]])
  datasc[[i]] <- ScaleData(datasc[[i]], features = all.genes)
  
## Dimensional reduction
  datasc[[i]] <- RunPCA(datasc[[i]], features = VariableFeatures(object = datasc[[i]]))
  
## Perform Clustering and UMAP
  
  
  datasc[[i]] <- FindNeighbors(datasc[[i]], dims = 1:10)
  datasc[[i]] <- FindClusters(datasc[[i]], resolution = 0.5)
  
  
  datasc[[i]] <- RunUMAP(datasc[[i]], dims = 1:10)
  
  
  
## Optional: Save Seurat Objects  
  #saveRDS(datasc[[i]],paste('write_R/',i,'.rds',sep=''))
  
}




# Integration -------------------------------------------------------------



## Select Integration Features

features <- SelectIntegrationFeatures(object.list = datasc)

## Find Anchors

anchors<- FindIntegrationAnchors(object.list=datasc,anchor.features=features)

## Create Integrated Data

scdata<-IntegrateData(anchorset = anchors)

## Optional Save
saveRDS(scdata,'write_R/scdata.rds')

