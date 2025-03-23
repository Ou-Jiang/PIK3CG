library(Seurat)
library(tidyverse)
# 假设sce是已分好细胞群体的Seurat对象
sce<-readRDS('major_type.rds')
# Step 1: 提取目标基因表达数据
gene <- "PIK3CG"
if (gene %in% rownames(sce)) {
  sce$PIK3CG_expression <- FetchData(sce, vars = gene)[, 1]
  
  # Step 2: 根据表达值进行分组
  # 设定一个阈值，例如表达值大于1认为高表达
 # expression_threshold <- 1
  sce$PIK3CG_group <- ifelse(sce$PIK3CG_expression==0, 
                                    "PIK3CG-","PIK3CG+")
  
  # Step 3: 合并现有群体信息，形成新的分群变量
  sce$merged_PIK3CG_group <- paste(sce$merged_type, sce$PIK3CG_group, sep = "_")
  
  # Step 4: 可视化或下游分析
  # 比如使用UMAP展示新的分组
  Idents(sce)<-'merged_PIK3CG_group'
  p<-DimPlot(sce, split.by = "PIK3CG_group", label = TRUE,repel=TRUE,label.size = 4) + ggtitle("PIK3CG Expression Groups")
  ggsave('umap_PIK3CG_group.jpg',width=20,height=9)
  ggsave('umap_PIK3CG_group.pdf',width=20,height=9)

} else {
  message("PIK3CG gene not found in the Seurat object.")
}

saveRDS(sce,'major_type.rds')
