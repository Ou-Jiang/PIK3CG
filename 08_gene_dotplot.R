library(Seurat)
library(tidyverse)
library(patchwork)

scRNA1<-readRDS('major_type.rds')
DefaultAssay(scRNA1)<-'RNA'

gene<-c('PIK3CA','PIK3CB','PIK3CD','PIK3CG')
p<-DotPlot(scRNA1,features=gene)+RotatedAxis()
ggsave('gene_dotplot_major.jpg',width=5,height=7)
ggsave('gene_dotplot_major.pdf',width=5,height=7)

scRNA1<-readRDS('fine_type.rds')
DefaultAssay(scRNA1)<-'RNA'

gene<-c('PIK3CA','PIK3CB','PIK3CD','PIK3CG')
p<-DotPlot(scRNA1,features=gene)+RotatedAxis()
ggsave('gene_dotplot_fine.jpg',width=5,height=7)
ggsave('gene_dotplot_fine.pdf',width=5,height=7)

scRNA1<-readRDS('major_type.rds')
DefaultAssay(scRNA1)<-'RNA'
Idents(scRNA1)<-'merged_type'

gene<-c('PIK3CA','PIK3CB','PIK3CD','PIK3CG')
p<-DotPlot(scRNA1,features=gene)+RotatedAxis()
ggsave('gene_dotplot_merged.jpg',width=5,height=7)
ggsave('gene_dotplot_merged.pdf',width=5,height=7)
