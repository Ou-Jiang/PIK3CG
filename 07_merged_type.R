library(Seurat)
sce_sub<-readRDS('fine_type.rds')
sce_total<-readRDS('major_type.rds')
library(tidyverse)
total_meta<-sce_total@meta.data
sub_meta<-sce_sub@meta.data
sub_meta$cellid<-rownames(sub_meta)
total_meta$cellid<-rownames(total_meta)
merged_meta<-merge(total_meta,sub_meta,by='cellid',all.x=TRUE)

# 提取需要的列，添加细胞编号列
sub_meta$cellid <- rownames(sub_meta)
total_meta$cellid <- rownames(total_meta)
# 从sub_meta提取细胞亚类
sub_fine_type <- sub_meta[, c("cellid", "fine_type")]
colnames(sub_fine_type) <- c("cellid", "Cell_Type")  # 统一列名
# 从total_meta提取细胞大类
total_major_type <- total_meta[, c("cellid", "major_type")]
colnames(total_major_type) <- c("cellid", "Cell_Type")  # 统一列名
# 合并两个数据框，优先sub_meta中的信息
merged_meta <- merge(total_major_type, sub_fine_type, by = "cellid", all.x = TRUE)
# 创建新列，优先使用细胞亚类（从sub_meta中）
merged_meta$Final_Cell_Type <- ifelse(!is.na(merged_meta$Cell_Type.y),
                                      merged_meta$Cell_Type.y,
                                      merged_meta$Cell_Type.x)
# 清理结果，保留最终所需的列
merged_meta <- merged_meta[, c("cellid", "Final_Cell_Type")]
library(dplyr)
# 提取细胞亚类和细胞大类，并统一列名
sub_fine_type <- sub_meta %>%
  select(cellid = cellid, Cell_Type = fine_type)
total_major_type <- total_meta %>%
  select(cellid = cellid, Cell_Type = major_type)
# 确保相关列是字符类型
sub_meta$fine_type <- as.character(sub_meta$fine_type)
total_meta$major_type <- as.character(total_meta$major_type)
# 使用merge方法
sub_meta$cellid <- rownames(sub_meta)
total_meta$cellid <- rownames(total_meta)
sub_fine_type <- sub_meta[, c("cellid", "fine_type")]
colnames(sub_fine_type) <- c("cellid", "Cell_Type")
total_major_type <- total_meta[, c("cellid", "major_type")]
colnames(total_major_type) <- c("cellid", "Cell_Type")
merged_meta <- merge(total_major_type, sub_fine_type, by = "cellid", all.x = TRUE)
merged_meta$Final_Cell_Type <- ifelse(!is.na(merged_meta$Cell_Type.y),
                                      merged_meta$Cell_Type.y,
                                      merged_meta$Cell_Type.x)
# 确保最终列为字符类型
merged_meta$Final_Cell_Type <- as.character(merged_meta$Final_Cell_Type)
merged_meta <- merged_meta[, c("cellid", "Final_Cell_Type")]
rownames(merged_meta)<-merged_meta$cellid
sce_merged<-AddMetaData(sce_total,metadata=merged_meta)

colnames(sce_merged@meta.data)[22]<-'merged_type'
sce_merged$cellid<-NULL
#saveRDS(sce_merged,'major_type.rds')

Idents(sce_merged)<-'merged_type'
DimPlot(sce_merged,label=T,repel=T)
p<-DimPlot(sce_merged,label=T,repel=T)
ggsave('umap_merged.pdf')
ggsave('umap_merged.jpg')
