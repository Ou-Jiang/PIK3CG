library(Seurat)
library(tidyverse)
library(patchwork)
#######library(future)
#######plan('multisession',workers=10)
#######options(future.globals.maxSize = 20000 * 1024^2)

scRNAlist<-readRDS('scRNAlist.rds')
scRNAlist<-scRNAlist[-c(1,7,8)]
#integrate seurat objects
for (i in 1:length(scRNAlist)) {
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst")
}

#########plan('sequential')
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)
scRNA3 <- IntegrateData(anchorset = scRNA.anchors)

#########plan('multisession',workers=4)
##########options(future.globals.maxSize = 64000 * 1024^2)

#caculate mitochondria content
scRNA <- scRNA3
DefaultAssay(scRNA) <- "RNA"
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

#quality control
violin <-VlnPlot(scRNA, group.by = "orig.ident",
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave('qc_vlnplot.jpg')
ggsave('qc_vlnplot.pdf')

violin <-VlnPlot(scRNA, group.by = "orig.ident",
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0)
ggsave('pt0_qc_vlnplot.jpg')
ggsave('pt0_qc_vlnplot.pdf')

scRNA <- subset(scRNA, subset =  nFeature_RNA < 6500 & percent.mt < 25)

#dimensionality reduction and clustering
scRNA1<-scRNA
DefaultAssay(scRNA1)<-'integrated'
scRNA1 <- ScaleData(scRNA1, features = VariableFeatures(scRNA1))
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1))
plot2 <- ElbowPlot(scRNA1, ndims=50, reduction="pca") 
ggsave('elbowplot.jpg')
ggsave('elbowplot.pdf')

pc.num=1:50
scRNA1 <- FindNeighbors(scRNA1, dims = pc.num) 
scRNA1 <- FindClusters(scRNA1, resolution = 0.5)
scRNA1 <- RunUMAP(scRNA1, dims = pc.num)
plot3 = DimPlot(scRNA1, reduction = "umap", label=T) 
plot4 = DimPlot(scRNA1, reduction = "umap", group.by='orig.ident')
plotc <- plot3+plot4
ggsave('umap_before_anno_0.5.jpg',width=14,height=7)
ggsave('umap_before_anno_0.5.pdf',width=14,height=7)

scRNA1 <- FindClusters(scRNA1, resolution = 1.0)
scRNA1 <- RunUMAP(scRNA1, dims = pc.num)
plot3 = DimPlot(scRNA1, reduction = "umap", label=T) 
plot4 = DimPlot(scRNA1, reduction = "umap", group.by='orig.ident')
plotc <- plot3+plot4
ggsave('umap_before_anno_1.0.jpg',width=14,height=7)
ggsave('umap_before_anno_1.0.pdf',width=14,height=7)

saveRDS(scRNA1,'before_anno.rds')
