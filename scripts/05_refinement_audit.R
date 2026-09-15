# 输入：聚类对象和再聚类细胞表；输出：逐亚簇 marker 检出比例与共表达证据。
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source("scripts/00_common.R", encoding = "UTF-8")
suppressPackageStartupMessages(library(Seurat))
stage_start("05_refinement_audit")
ensure_dirs("results/annotation")
o <- readRDS("data/processed/04_seurat_clustered.rds")
c <- read.csv("results/annotation/lymphocyte_refinement_cells.csv")
m <- gene_mapping(); x <- counts_of(o)
genes <- c("CD3D","CD3E","TRAC","KLRF1","FCGR3A","GNLY","EPCAM","KRT8","KRT19","TFF2","LYZ","C1QA","TYROBP","MS4A1","CD79A")
ids <- m$seurat_feature[match(genes,m$gene_symbol)]
z <- x[ids,c$cell_id] > 0; rownames(z) <- genes
out <- do.call(rbind,lapply(unique(c$sub_key),function(k) {
 s <- c$sub_key == k
 data.frame(sub_key=k,n_cells=sum(s), t(as.matrix(Matrix::rowMeans(z[,s,drop=FALSE]))),check.names=FALSE)
}))
names(out)[-(1:2)] <- genes
out$t_lineage_2_of_3 <- sapply(out$sub_key,function(k) mean(Matrix::colSums(z[c("CD3D","CD3E","TRAC"),c$sub_key==k,drop=FALSE]) >=2))
out$epithelial_2_of_3 <- sapply(out$sub_key,function(k) mean(Matrix::colSums(z[c("EPCAM","KRT8","KRT19"),c$sub_key==k,drop=FALSE]) >=2))
write.csv(out,"results/annotation/refinement_detection_audit.csv",row.names=FALSE)
saveRDS(out,"data/processed/05_refinement_audit.rds")
print(out)
stage_end()
