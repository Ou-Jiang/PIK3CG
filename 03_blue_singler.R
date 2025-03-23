library(Seurat)
scRNA1<-readRDS('before_anno.rds')
DefaultAssay(scRNA1)<-'RNA'
Idents(scRNA1)<-'integrated_snn_res.0.5'

library(celldex)
blue<-BlueprintEncodeData()
pbmc_for_singler<-GetAssayData(scRNA1,layer='data')

library(SingleR)
options(matrixStats.useNames.NA = "deprecated")
pbmc.blue<-SingleR(test=pbmc_for_singler,ref=blue,labels=blue$label.main)
scRNA1$blueprint_main<-pbmc.blue$labels

saveRDS(scRNA1,'before_anno.rds')
savehistory('blue_singler.R')
