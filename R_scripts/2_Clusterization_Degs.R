
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



# PCA + Clusterization ----------------------------------------------------

scdata<-readRDS('write_R/scdata.RDS')

## Set defaultassay to integrated
DefaultAssay(scdata) <- "integrated"

## ScaleData
scdata <- ScaleData(scdata, verbose = FALSE)


## PCA
scdata <- Seurat::RunPCA(scdata, npcs = 30, verbose = FALSE)
## UMAP
scdata <- RunUMAP(scdata, reduction = "pca", dims = 1:30)

#Clustering
scdata <- FindNeighbors(scdata, reduction = "pca", dims = 1:30)
scdata <- FindClusters(scdata, resolution = 0.1)



## Optional Visualize
DimPlot(scdata, reduction = "umap", label = TRUE, repel = TRUE)



# Metadata ----------------------------------------------------------------



## Add samples metadata

scdata@meta.data %>% 
  rownames_to_column('CellID') %>% 
  mutate(sample=orig.ident,
         Groups=case_when(sample%in%c(1,12,13)~'MABSC-PD',
                          sampel%in%c(3,4,9)~'MABSC',
                          sample%in%c(5,6,14,15)~'CF'),
         sample_f=factor(as.character(sample),levels = as.character(c(5,6,14,15,3,4,9,1,12,13))),
         Groups_f=factor(Groups,levels('CF','MABSC','MABSC-PD'))) %>% 
  column_to_rownames('CellID')->scdata@meta.data




# Find groups-independent conserved markers -------------------------------


DefaultAssay(scdata)<-'RNA'


sc_markers_list<-list()
for (i in 0:13) {
  idx=paste('cluster_',i,sep='')
  sc_markers_list[[idx]]<-FindConservedMarkers(scdata,ident.1=i,grouping.var = "Groups",verbose=T)
}


## Check markers 

bind_rows(sc_markers_list,.id='cluster') %>% 
  rownames_to_column('gene') %>% 
  mutate(gene=str_remove(gene,"\\.\\.\\..*"),
         cluster=(str_remove(cluster,'cluster_'))) %>% tibble() %>% 
  filter(!str_detect(gene,'^IB|^RP|^LINC')) %>% 
  filter(if_all(ends_with("avg_log2FC"),is.positive)) %>% 
  group_by(cluster) %>% 
  mutate(PAr=rank(-`PA-CS_avg_log2FC`),
         CSr=rank(-`NTM-CS_avg_log2FC`),
         CDr=rank(-`NTM-CD_avg_log2FC`),
         Rank=(PAr+CSr+CDr)/3,
         Rank=rank(Rank),
         cluster=as.numeric(cluster))->clustermarkers  
## Optional: save markers  
write_excel_csv('write_R/markers.csv',delim=';')


## Filter out non-classificable cells

scdata<-subset(scdata,seurat_clusters<10)

## Annotate celltype
scdata$celltype<-plyr::mapvalues(scdata$seurat_clusters,
                                 from = c(1,2,8,4,6,7,3,5,10,9),
                                 to=c("CD4+ Naive T","CD14+ Monocyte","CD4+ Effector","NK","CD16+ Monocyte","B Cell ","CD8+ Effector ","CD8+ Naive T","DC","Platelet"))




# Images EXPORT -----------------------------------------------------------

## Umap by groups




DimPlot(scdata,split.by = 'Groups_f',group.by = 'Groups_f',cols=c('green4','#e69f11','#b11600'))+
  ggtitle('')

ggsave('images/UMAP_groups.png',dpi=300,height = 5,width = 10)


## Umaps celltype

palcluster<-c("#0075DC","#9DCC00","#A588FF","#2BCE48","#005C31","#00998F","#808080","#C20088","#FF5005","#FF0010")



DimPlot(scdata,group.by = c('celltype'), reduction = "umap", label = TRUE, 
        repel = TRUE,label.box = T,cols=palcluster)

ggsave('images/UMAP_celltype.png',dpi=300,height=6,width=9)





# Differentially expression analysis --------------------------------------




## MABSC+ vs CF




Idents(scdata)<-'cellcluster'

scdata_s<-subset(scdata,sample!='2')

scdata_s$NTM<-ifelse(str_detect(scdata_s$Groups,'MABSC'),'MABSC','noMABSC')
scdata_s$cellcluster.cond<-paste(scdata_s$cellcluster,scdata_s$NTM,sep='_')
Idents(scdata_s)<-'cellcluster.cond'



Deg_list<-list()

for (i in 0:9) {
  #for (j in 1:dim(Conditions)[1]){
  a=paste(i,'MABSC',sep='_')
  b=paste(i,'noMABSC',sep='_')
  idx=paste(i,'MABSC+,vs,CF',sep='_')
  Deg_list[[idx]]<-FindMarkers(scdata_s,ident.1=a,ident.2=b,verbose=T,logfc.threshold = 0)%>% 
    as.data.frame() %>% rownames_to_column('Genes') %>% tibble()
  #}
}
remove(scdata_s)

## MABSC vs CF




scdata_s<-subset(scdata,Groups!='MABSC-PD')

scdata_s$NTM<-ifelse(str_detect(scdata_s$Groups,'MABSC'),'MABSC','noMABSC')
scdata_s$cellcluster.cond<-paste(scdata_s$cellcluster,scdata_s$NTM,sep='_')
Idents(scdata_s)<-'cellcluster.cond'





for (i in 0:9) {
  #for (j in 1:dim(Conditions)[1]){
  a=paste(i,'MABSC',sep='_')
  b=paste(i,'noMABSC',sep='_')
  idx=paste(i,'MABSC,vs,CF',sep='_')
  Deg_list[[idx]]<-FindMarkers(scdata_s,ident.1=a,ident.2=b,verbose=T,logfc.threshold = 0) %>% 
    as.data.frame() %>% rownames_to_column('Genes') %>% tibble()
  #}
}






## MABSC-PD vs CF+MABSC




scdata_s<-subset(scdata,sample!='2')

scdata_s$NTM<-ifelse(str_detect(scdata_s$Groups,'MABSC-PD'),'PD','CS')
scdata_s$cellcluster.cond<-paste(scdata_s$cellcluster,scdata_s$NTM,sep='_')
Idents(scdata_s)<-'cellcluster.cond'





for (i in 0:9) {
  #for (j in 1:dim(Conditions)[1]){
  a=paste(i,'PD',sep='_')
  b=paste(i,'CS',sep='_')
  idx=paste(i,'MABSC-PD,vs,CF + MABSC',sep='_')
  Deg_list[[idx]]<-FindMarkers(scdata_s,ident.1=a,ident.2=b,verbose=T,logfc.threshold = 0) %>% 
    as.data.frame() %>% rownames_to_column('Genes') %>% tibble()
  #}
}





## MABSC-PD vs CF




scdata_s<-subset(scdata,sample!='2'&Groups!='MABSC')

scdata_s$NTM<-ifelse(str_detect(scdata_s$Groups,'MABSC-PD'),'MABSC-PD','CS')
scdata_s$cellcluster.cond<-paste(scdata_s$cellcluster,scdata_s$NTM,sep='_')
Idents(scdata_s)<-'cellcluster.cond'





for (i in 0:9) {
  #for (j in 1:dim(Conditions)[1]){
  a=paste(i,'MABSC-PD',sep='_')
  b=paste(i,'CS',sep='_')
  idx=paste(i,'MABSC-PD,vs,CF',sep='_')
  Deg_list[[idx]]<-FindMarkers(scdata_s,ident.1=a,ident.2=b,verbose=T,logfc.threshold = 0) %>% 
    as.data.frame() %>% rownames_to_column('Genes') %>% tibble()
  #}
}







## MABSC-PD vs MABSC




scdata_s<-subset(scdata,sample!='2'&Groups!='CF')

scdata_s$NTM<-ifelse(str_detect(scdata_s$Groups,'PD'),'CD','CS')
scdata_s$cellcluster.cond<-paste(scdata_s$cellcluster,scdata_s$NTM,sep='_')
Idents(scdata_s)<-'cellcluster.cond'





for (i in 0:9) {
  #for (j in 1:dim(Conditions)[1]){
  a=paste(i,'CD',sep='_')
  b=paste(i,'CS',sep='_')
  idx=paste(i,'MABSC-PD,vs,MABSC',sep='_')
  Deg_list[[idx]]<-FindMarkers(scdata_s,ident.1=a,ident.2=b,verbose=T,logfc.threshold = 0) %>% 
    as.data.frame() %>% rownames_to_column('Genes') %>% tibble()
  #}
}

remove(scdata_s)


## Output Degslist
bind_rows(Deg_list,.id='aa') %>% 
  separate(aa,c('Cluster','Comparison'),sep='_') %>% 
  write_delim('write_R/scDeg.csv',delim=';')










