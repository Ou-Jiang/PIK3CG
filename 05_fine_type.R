library(Seurat)
library(tidyverse)
library(patchwork)

sce_interest<-readRDS('fine_type.rds')

fine_type<-c('CD8+Tm',
'CD8+Tex','Treg','M2c','CD8+Tm',
'M1','M2a','M2a','CD8+Tm',
'CD8+Tex','CD8+Tm','M1','M2a',
'M1','CD8+Tm','M1','M2a',
'CD8+Teff','CD8+Tn')
Idents(sce_interest)<-'RNA_snn_res.0.5'
names(fine_type)<-levels(sce_interest)
sce_interest<-RenameIdents(sce_interest,fine_type)
sce_interest@meta.data$fine_type<-Idents(sce_interest)

DefaultAssay(sce_interest)<-'RNA'

p1<-DimPlot(sce_interest,reduction='umap',label=TRUE,repel=TRUE)
ggsave('umap_fine.jpg')
ggsave('umap_fine.pdf')

marker_fine<-readRDS('marker_fine.rds')
p<-DotPlot(sce_interest,features=marker_fine,group.by='fine_type')+RotatedAxis()
ggsave('dotplot_qc_fine.jpg',width=12,height=8)
ggsave('dotplot_qc_fine.pdf',width=12,height=8)


