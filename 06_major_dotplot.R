library(Seurat)
library(tidyverse)

sce<-readRDS('major_type.rds')
marker_main<-readRDS('marker_main.rds')
levels(sce)<-c('Treg','CD4T','CD8T','B cell','NK','mast cell','neutrophil',
'monocyte','macrophage','epithelial','endothelial/fibroblast')
p<-DotPlot(sce,features=marker_main)+RotatedAxis()
ggsave('dotplot_qc.jpg',width=10,height=7)
ggsave('dotplot_qc.pdf',width=10,height=7)
