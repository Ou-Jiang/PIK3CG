library(Seurat)
library(tidyverse)
# 假设 sce 已经包含 merged_PIK3CG_group 变量
sce<-readRDS('major_type.rds')
# Step 1: 提取所有细胞群体
all_clusters <- unique(sce$merged_type)

# Step 2: 创建存储差异分析结果的列表
diff_results <- list()

# Step 3: 批量执行差异分析
for (cluster in all_clusters) {
  # 提取当前群体的细胞
  cells_in_cluster <- WhichCells(sce, expression = merged_type == cluster)
  
  # 提取当前群体细胞的 Seurat 子集
  cluster_subset <- subset(sce, cells = cells_in_cluster)
  
  # 设置分组标识为 PIK3CG_group
  Idents(cluster_subset) <- "PIK3CG_group"
  
  # 确保群体中存在阳性和阴性两种细胞
  if (length(unique(Idents(cluster_subset))) == 2) {
    # 差异分析：PIK3CG_high vs PIK3CG_low
    markers <- FindMarkers(cluster_subset, ident.1 = "PIK3CG+", ident.2 = "PIK3CG-")
    
    # 将结果存储到列表中
    diff_results[[cluster]] <- markers
  } else {
    message(paste("Skipping cluster", cluster, "- does not contain both PIK3CG+ and PIK3CG+ cells."))
  }
}

# Step 4: 将结果写入文件（可选）
for (cluster in names(diff_results)) {
  write.csv(diff_results[[cluster]], if(cluster=="endothelial/fibroblast"){"DiffMarkers_endothelial_fibroblast.csv"}else{paste0("DiffMarkers_", cluster, ".csv")})
}







library(clusterProfiler)
library(org.Hs.eg.db) # 假设分析的是人类数据
library(ggplot2)

# Step 1: 加载差异分析结果
# 假设 diff_results 是前面差异分析生成的结果列表
# 创建存储GO富集分析结果的列表
up_results <- list()

# Step 2: 对每个群体执行GO富集分析
for (cluster in names(diff_results)) {
  # 提取差异分析结果
  markers <- diff_results[[cluster]]
  
  # 筛选显著差异基因（例如 P.Val < 0.05 ）
  up_genes <- rownames(markers[markers$p_val < 0.05 & markers$avg_log2FC > 0, ])
  
  # 确保有足够的显著基因
  if (length(up_genes) > 10) {
    # 将基因名称转换为ENTREZ ID
    entrez_genes <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    # GO富集分析
    ego <- enrichGO(
      gene = entrez_genes$ENTREZID,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "ALL", # Biological Process
      pAdjustMethod = "none",
      pvalueCutoff = 0.05,
      qvalueCutoff = 1,
      readable = TRUE
    )
    
    # 将结果存储
    up_results[[cluster]] <- ego
   

    # 可视化前10个GO条目
    if (!is.null(ego)) {
	p<-dotplot(ego, split="ONTOLOGY",showCategory=10,color="pvalue")+
  facet_grid(ONTOLOGY~., scale="free")+
ggtitle(paste("GO Enrichment for", paste("PIK3CG+",cluster)))
 ggsave(if(cluster=="endothelial/fibroblast"){"GOenrich_endothelial_fibroblast_PIK3CG+.jpg"}else{paste0("GOenrich_", paste0(cluster,"_PIK3CG+"), ".jpg")},width = 9,height = 14)
ggsave(if(cluster=="endothelial/fibroblast"){"GOenrich_endothelial_fibroblast_PIK3CG+.pdf"}else{paste0("GOenrich_", paste0(cluster,"_PIK3CG+"), ".pdf")},width = 9,height = 14)
    }
  } else {
    message(paste("Skipping cluster", cluster, "- not enough significant genes."))
  }
}

# Step 3: 保存GO分析结果
# 将每个群体的GO富集结果保存到CSV文件
for (cluster in names(up_results)) {
  if (!is.null(up_results[[cluster]])) {
    write.csv(as.data.frame(up_results[[cluster]]), if(cluster=="endothelial/fibroblast"){"GO_Results_endothelial_fibroblast_PIK3CG+.csv"}else{paste0("GO_Results_", paste0(cluster,"_PIK3CG+"), ".csv")})
  }
}

# 如果需要提取和展示所有群体的富集结果，可以整合到一个表中
#combined_results <- do.call(rbind, lapply(names(up_results), function(cluster) {
#  if (!is.null(up_results[[cluster]])) {
#    df <- as.data.frame(up_results[[cluster]])
#    df$Cluster <- cluster
#    return(df)
#  }
#}))
#write.csv(combined_results, "Combined_GO_Results.csv")


# step 1: 加载差异分析结果
# 假设 diff_results 是前面差异分析生成的结果列表
# 创建存储GO富集分析结果的列表
down_results <- list()

# Step 2: 对每个群体执行GO富集分析
for (cluster in names(diff_results)) {
  # 提取差异分析结果
  markers <- diff_results[[cluster]]

  # 筛选显著差异基因（例如 P.Val < 0.05 ）
  down_genes <- rownames(markers[markers$p_val < 0.05 & markers$avg_log2FC < 0, ])

  # 确保有足够的显著基因
  if (length(down_genes) > 10) {
    # 将基因名称转换为ENTREZ ID
    entrez_genes <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    # GO富集分析
    ego <- enrichGO(
      gene = entrez_genes$ENTREZID,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "ALL", # Biological Process
      pAdjustMethod = "none",
      pvalueCutoff = 0.05,
      qvalueCutoff = 1,
      readable = TRUE
    )

    # 将结果存储
    down_results[[cluster]] <- ego
    # 可视化前10个GO条目
    if (!is.null(ego)) {
        p<-dotplot(ego, split="ONTOLOGY",showCategory=10,color="pvalue")+
  facet_grid(ONTOLOGY~., scale="free")+
ggtitle(paste("GO Enrichment for", paste("PIK3CG-",cluster)))
 ggsave(if(cluster=="endothelial/fibroblast"){"GOenrich_endothelial_fibroblast_PIK3CG-.jpg"}else{paste0("GOenrich_", paste0(cluster,"_PIK3CG-"), ".jpg")},width = 9,height = 14)
ggsave(if(cluster=="endothelial/fibroblast"){"GOenrich_endothelial_fibroblast_PIK3CG-.pdf"}else{paste0("GOenrich_", paste0(cluster,"_PIK3CG-"), ".pdf")},width = 9,height = 14)
    }
  } else {
    message(paste("Skipping cluster", cluster, "- not enough significant genes."))
  }
}

# Step 3: 保存GO分析结果
# 将每个群体的GO富集结果保存到CSV文件
for (cluster in names(down_results)) {
  if (!is.null(down_results[[cluster]])) {
    write.csv(as.data.frame(down_results[[cluster]]), if(cluster=="endothelial/fibroblast"){"GO_Results_endothelial_fibroblast_PIK3CG-.csv"}else{paste0("GO_Results_", paste0(cluster,"_PIK3CG-"), ".csv")})
  }
}



