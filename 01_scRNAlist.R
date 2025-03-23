library(Seurat)
library(data.table)
data<-fread('GSE207422_NSCLC_scRNAseq_UMI_matrix.txt')
data<-as.data.frame(data)
rownames(data)<-data[,1]
data<-data[,-1]
sce<-CreateSeuratObject(data)

library(readxl)
meta<-read_excel('GSE207422_NSCLC_scRNAseq_metadata.xlsx')
meta<-meta[1:15,]
meta<-as.data.frame(meta)

library(tidyverse)
orig<-str_split(colnames(sce),'_',simplify=T)
orig<-paste(orig[,1],orig[,2],sep='_')
sce$orig.ident<-orig

orig<-data.frame(Sample=orig)
orig_add<-left_join(orig,meta)
rownames(orig_add)<-colnames(sce)
orig_add<-orig_add[,-1]

sce<-AddMetaData(sce,orig_add)
scRNAlist<-SplitObject(sce,split.by='orig.ident')
scRNAlist<-saveRDS(scRNAlist,'scRNAlist.rds')
