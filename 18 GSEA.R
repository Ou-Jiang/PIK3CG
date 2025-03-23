library(org.Hs.eg.db)
library(stringr)
library(BiocGenerics)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(enrichplot)
library(future)
library(future.apply)
library(conflicted)


setwd("F:/r/GSE30784")    

dir.create("result")

exprSet <- read.csv("GSE30784.csv",row.names = 1)


batch_cor <- function(gene){
  y <- as.numeric(exprSet[gene,])
  rownames <- rownames(exprSet)
  do.call(rbind,future_lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(exprSet[x,]),y,type="spearman")
    data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}


dd <- batch_cor("PIK3CG")

write.csv(dd,"结果/PIK3CG.spearman-cor.csv")

gene <- dd$mRNAs


gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_df <- data.frame(cor=dd$cor,SYMBOL = dd$mRNAs)

gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

geneList <- gene_df$cor
names(geneList)=gene_df$SYMBOL

geneList=sort(geneList,decreasing = T)

dotplot_internal <- function(object, x = "GeneRatio", color = "pvalue",
                             showCategory=10, size=NULL, split = NULL,
                             font.size=12, title = "", orderBy="x", decreasing=TRUE) {
 
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    if (is.null(size))
      size <- "Count"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
    if (is.null(size))
      size <- "GeneRatio"
  } else if (is(x, "formula")) {
    x <- as.character(x)[2]
    if (is.null(size))
      size <- "Count"
  } else {
    if (is.null(size))
      size  <- "Count"
  }
  df <- fortify(object, showCategory = showCategory, split=split)
  
  if (orderBy !=  'x' && !orderBy %in% colnames(df)) {
    message('wrong orderBy parameter; set to default `orderBy = "x"`')
    orderBy <- "x"
  }
  
  if (orderBy == "x") {
    df <- dplyr::mutate(df, x = eval(parse(text=x)))
  }
  
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description, levels=rev(unique(df$Description[idx])))
  ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
    geom_point() +
    scale_color_continuous(low="darkorange3", high="darkseagreen4", name = color, guide=guide_colorbar(reverse=TRUE)) +
    ylab(NULL) + ggtitle(title) + DOSE::theme_dose(font.size) + scale_size(range=c(3, 8))
  
}

GO<-read.gmt("c5.go.v2023.2.Hs.symbols.gmt")

GO <-GSEA(geneList,TERM2GENE = GO)

library(ggplot2)
library(enrichplot)

write.csv(data.frame(GO), file = "ruslt.csv", row.names = FALSE)


p1=gseaplot2(GO,
          geneSetID = 4,
          title=GO$Description[4],
          color="red", 
          base_size = 20, 
          subplots = 1:3, 
          pvalue_table = T) 
p1
ggsave("result/p1.pdf",width = 12,height = 10)

p2=gseaplot2(GO,
          geneSetID = 1:3,#这里展示前3个通路
          subplots = 1:3,
          ES_geom='line'
)
p2
ggsave("result/p2.pdf",width = 12,height = 10)

num=5
p3=gseaplot2(GO, geneSetID = rownames(GO@result)[head(order(GO@result$enrichmentScore),num)])
p3
ggsave("result/p3.pdf",width = 12,height = 10)

p4=gseaplot2(GO, geneSetID = rownames(GO@result)[tail(order(GO@result$enrichmentScore),num)])
p4
ggsave("result/p4.pdf",width = 12,height = 10)

num=5
p5=gseaplot2(GO, geneSetID = rownames(GO@result)[c(head(order(GO@result$enrichmentScore),num),tail(order(GO@result$enrichmentScore),num))])
p5
ggsave("result/p5.pdf",width = 12,height = 12)


p6=dotplot(GO,color="pvalue")
p6
ggsave("result/p6.pdf",width = 10,height = 10)


p7=dotplot(GO,split=".sign")+facet_grid(~.sign)
p7
ggsave("result/p7.pdf",width = 10,height = 10)

p8=ridgeplot(GO,showCategory = 10,
          fill = 'p.adjust',
          label_format = 25,
          core_enrichment = T)
p8
ggsave("result/p8.pdf",width = 10,height = 10)

library(cowplot)
pp <- lapply(1:3, function(i) {
  anno <- GO[i, c("NES", "pvalue", "p.adjust")]
  lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
  
  gsearank(GO, i, GO[i, 2]) + xlab(NULL) +ylab(NULL) +
    annotate("text", 10000, GO[i, "enrichmentScore"] * .75, label = lab, hjust=0.5, vjust=0.5)
})
p9=plot_grid(plotlist=pp, ncol=1)
p9
ggsave("result/p9.pdf",width = 10,height = 10)
