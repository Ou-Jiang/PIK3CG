library(Seurat)
library(tidyverse)
library(patchwork)
scRNA1<-readRDS('major_type.rds')
DefaultAssay(scRNA1)<-'RNA'
gene<-c('PIK3CA','PIK3CB','PIK3CD','PIK3CG')
p<-VlnPlot(scRNA1,features = gene,stack=T,pt.size=0,flip = T,add.noise = T)+
theme(axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
      legend.position = 'none')
ggsave('gene_vlnplot_major.jpg',height=15,width=10)
ggsave('gene_vlnplot_major.pdf',height=15,width=10)

Idents(scRNA1)<-'merged_type'
gene<-c('PIK3CA','PIK3CB','PIK3CD','PIK3CG')
p<-VlnPlot(scRNA1,features = gene,stack=T,pt.size=0,flip = T,add.noise = T)+
theme(axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
      legend.position = 'none')
ggsave('gene_vlnplot_merged.jpg',height=15,width=10)
ggsave('gene_vlnplot_merged.pdf',height=15,width=10)

scRNA1<-readRDS('fine_type.rds')
DefaultAssay(scRNA1)<-'RNA'
gene<-c('PIK3CA','PIK3CB','PIK3CD','PIK3CG')
p<-VlnPlot(scRNA1,features = gene,stack=T,pt.size=0,flip = T,add.noise = T)+
theme(axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
      legend.position = 'none')
ggsave('gene_vlnplot_fine.jpg',height=15,width=10)
ggsave('gene_vlnplot_fine.pdf',height=15,width=10)
