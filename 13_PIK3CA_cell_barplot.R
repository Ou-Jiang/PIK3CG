library(Seurat)
library(tidyverse)

sce<-readRDS('major_type.rds')
# 确保元数据包含必要的列
if (!all(c("merged_type", "RECIST") %in% colnames(sce@meta.data))) {
  stop("元数据中缺少 'merged_type' 或 'RECIST' 列")
}

# 定义目标基因
gene <- "PIK3CA"
if (!(gene %in% rownames(sce))) {
  stop("目标基因不存在于 Seurat 对象中")
}

# 初始化结果列表
results_list <- list()

# 按细胞类群分组
merged_types <- unique(sce$merged_type)

for (merged_type in merged_types) {
  message(paste("正在分析细胞亚群:", merged_type))
  
  # 筛选当前细胞类群的细胞
  subset_cells <- WhichCells(sce, expression = merged_type == !!merged_type)
  subset_obj <- subset(sce, cells = subset_cells)
  
  # 提取基因表达和元数据
  subset_obj$gene_expression <- FetchData(subset_obj, vars = gene)[, 1]
  subset_data <- subset_obj@meta.data
  subset_data$Cell_ID <- rownames(subset_data)
  
  # 检验组间差异
  expr_by_group <- split(subset_obj$gene_expression, subset_obj$RECIST)
  if (length(expr_by_group) == 2) {
    test_result <- t.test(expr_by_group[[1]], expr_by_group[[2]])
    p_value <- test_result$p.value
  } else {
    test_result <- kruskal.test(subset_obj$gene_expression ~ subset_obj$RECIST)
    p_value <- test_result$p.value
  }
  
  # 添加显著性列
  subset_data$Significance <- ifelse(p_value < 0.05, "Significant", "Not Significant")
  
  # 添加目标基因表达量
  subset_data$Gene_Expression <- subset_obj$gene_expression
  
  # 合并到结果列表
  results_list[[merged_type]] <- subset_data[, c("Cell_ID", "RECIST", "merged_type", "Gene_Expression", "Significance")]
}

# 合并所有亚群的数据
final_results <- do.call(rbind, results_list)

# 重命名列
colnames(final_results) <- c("Cell_ID", "Immune_Response", "Cell_Type", "Gene_Expression", "Significance")

# 保存为CSV文件
#write.csv(final_results, "ImmuneResponse_Analysis_Results.csv", row.names = FALSE)

# 输出结果数据框
head(final_results)



dat2<-final_results
#dat2$Significance<-NULL
colnames(dat2)<-c('sample','group','ImmuneCell','Score','change')

library(rstatix)
stat_res <- dat2 %>% 
  group_by(ImmuneCell) %>% 
  t_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH")
stat_res$p1 <- signif(stat_res$p, 3)
stat_res$p2 <- ifelse(stat_res$p1<0.05, ifelse(stat_res$p1<0.01,ifelse(stat_res$p1<0.001,ifelse(stat_res$p1<0.0001,"****","***"),"**"),"*"), "ns")

#绘图
library(ggpubr)
p <- ggbarplot(dat2, x = "ImmuneCell", y = "Score", color = 'group',fill = NULL,
               outlier.shape = NA,
               bxp.errorbar = T,
               add = "mean_se",position = position_dodge(0.8)) +
  scale_colour_manual(values = c(SD="orange",PR="navy"))+
  stat_pvalue_manual(stat_res,
                     x = "ImmuneCell",
                     y.position = 0.6,
                     size = 3.5,
                     color = "black",
                     family = "Times",
                     label = "p2",
                     parse = F) +
  labs(x = "", y = "Gene Expression",color="Group") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 20, face = "bold", family = "Times"),
        axis.title.x =element_blank(),
        axis.text = element_text(size = 14, face = "bold", family = "Times"),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,vjust = 1),
        legend.text = element_text(size = 16, family = "Times"),
        text = element_text(family = "Times"),
        plot.title = element_text(hjust = 0.5,size = 20,face = "bold",family = "Times"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=15))+
   labs(title="PIK3CA")
ggsave(filename = "PIK3CA_cell_barplot.pdf", height = 5, width = 12)
ggsave(filename = "PIK3CA_cell_barplot.jpg", height = 5, width = 12)

