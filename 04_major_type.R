library(Seurat)
library(tidyverse)
library(patchwork)

scRNA1<-readRDS('before_anno.rds')
major_type<-c('macrophage',
'CD8T','B cell','CD8T','epithelial',
'monocyte','macrophage','epithelial','monocyte',
'mast cell','endothelial/fibroblast','epithelial','neutrophil',
'epithelial','macrophage','monocyte','macrophage',
'B cell','CD8T','CD4T','epithelial','NK','Treg','CD4T')
Idents(scRNA1)<-'integrated_snn_res.0.5'
names(major_type)<-levels(scRNA1)
scRNA1<-RenameIdents(scRNA1,major_type)
scRNA1@meta.data$major_type<-Idents(scRNA1)

DefaultAssay(scRNA1)<-'RNA'

p1<-DimPlot(scRNA1,reduction='umap',label=TRUE,repel=TRUE)
ggsave('umap_major.jpg')
ggsave('umap_major.pdf')

p2<-DimPlot(scRNA1,reduction='umap',label=TRUE,repel=TRUE,group.by='integrated_snn_res.0.5')
p<-p1+p2
ggsave('umap_combine_major_integ.jpg',width=14)
ggsave('umap_combine_major_integ.pdf',width=14)

marker_main<-readRDS('marker_main.rds')
p<-DotPlot(scRNA1,features=marker_main,group.by='major_type')+RotatedAxis()
ggsave('dotplot_qc.jpg',width=12,height=8)
ggsave('dotplot_qc.pdf',width=12,height=8)

