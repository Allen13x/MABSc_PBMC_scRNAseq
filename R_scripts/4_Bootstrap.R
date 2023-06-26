
# Import libraries --------------------------------------------------------


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


# Prepare bootstrap -------------------------------------------------------

## Import data
scdata<-readRDS('write_R/scdata.rds')




Idents(scdata)<-'cellcluster'


final_sub<-subset(scdata,sample!='2')


comparison_listA<-list()
comparison_listB<-list()


comparison_listA[['MABSC+_vs_CF']]<-c('MABSC-PD','MABSC')
comparison_listA[['MABSC_vs_CF-CS']]<-c('MABSC')
comparison_listA[['MABSC-PD_vs_CF_+_MABSC']]<-c('MABSC-PD')
comparison_listA[['MABSC-PD_vs_CF']]<-c('MABSC-PD')
comparison_listA[['MABSC-PD_vs_MABSC']]<-c('MABSC-PD')


comparison_listB[['CF-NTM+_vs_CF-CS']]<-c('CF')
comparison_listB[['CF-MABSC_vs_CF-CS']]<-c('CF')
comparison_listB[['CF-MABSC-PD_vs_CF-CS_+_CF-MABSC']]<-c('CF','MABSC')
comparison_listB[['CF-MABSC-PD_vs_CF-CS']]<-c('CF')
comparison_listB[['CF-MABSC-PD_vs_CF-MABSC']]<-c('MABSC')





boot_min_deg<-list()
for (i in names(comparison_listA)){
  
  #i='CF-NTM-CD_vs_CF-CS_+_CF-NTM-CS'
  
  
  final_sub$comparison<-paste(
    case_when(final_sub$Groups%in%comparison_listA[[i]]~'A',
              final_sub$Groups%in%comparison_listB[[i]]~'B'),final_sub$cellcluster,sep='_')
  
  
  
  table(final_sub$sample,final_sub$comparison)
  
  
  
  for (j in 0:9) {
    
    final_subs<-subset(final_sub,comparison==paste('A',j,sep='_')|comparison==paste('B',j,sep='_'))
    
    Idents(final_subs)<-'comparison'
    
    table(final_subs$sample,final_subs$comparison) %>% as.data.frame() %>% 
      filter(Freq>0) %>% 
      summarise(m=min(Freq)) %>% pull(m)->m
    
    table(final_subs$sample,final_subs$comparison) %>% as.data.frame() %>% 
      filter(Freq>0) %>% pull(Var1) ->sub_sample
    for (z in 1:20){ 
      
      paste(i,j,z,sep='__')->b_index
      old_sub_cell=sub_cell 
      print(b_index)
      #    x <- sample(1:1000,1)
      set.seed(z)
      final_subs@meta.data %>% 
        rownames_to_column('cell') %>% 
        filter(sample%in%sub_sample) %>% 
        group_by(sample) %>% 
        slice_sample(n=m) %>%ungroup() %>%  pull(cell)->sub_cell
      
      
      if(setequal(old_sub_cell, sub_cell)){
        
        break}
      
      
      boot_min_deg[[b_index]]<-FindMarkers(final_subs[,sub_cell],
                                           ident.1=paste('A',j,sep='_'),ident.2=paste('B',j,sep='_'),
                                           verbose=T)
      
      
      
      
    }
    
  }
  
  
  
}




## Visualize 


bind_rows(boot_min_deg,.id='boot') %>% rownames_to_column('genes') %>% 
  separate(boot,c('comparison','cellcluster','boot'),sep='__') %>% 
  filter(p_val_adj<0.05) %>% 
  group_by(comparison,cellcluster,boot) %>% 
  summarise(n=n()) %>% 
  ggplot(aes(x=cellcluster,y=n,fill=cellcluster)) +
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values=palcluster,labels=data.frame(scdata$cellcluster,scdata$celltype) %>% tibble() %>% 
                      distinct() %>% arrange(scdata.cellcluster) %>% pull(scdata.celltype))+
  facet_wrap(~comparison,scales = 'free')

ggsave('image/subsampling_DEGS.png',dpi=300,height=7,width=10)