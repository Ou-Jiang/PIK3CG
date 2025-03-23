library(Seurat)
sce<-readRDS('major_type.rds')

expr_df<-data.frame(type=paste(sce$orig.ident,sce$merged_type,sep='.'),
pik3cg=sce$PIK3CG_expression,pik3ca=sce$PIK3CA_expression,
pik3cb=sce$PIK3CB_expression,pik3cd=sce$PIK3CD_expression)

write.csv(expr_df,'patient.merged_type_expr.csv')
