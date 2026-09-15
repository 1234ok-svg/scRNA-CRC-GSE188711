# 输入：聚类对象、marker 表；输出：marker 支持图、候选注释评分、经审阅的细胞注释对象。
# 默认先输出证据；--apply 读取人工审阅后的 cluster_annotations.csv。
# 注释不会依据样本左右侧；不把所有上皮细胞直接称为恶性细胞。
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source("scripts/00_common.R", encoding = "UTF-8")
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)})
stage_start("05_annotation")
ensure_dirs("results/annotation", "figures/annotation")
object <- readRDS("data/processed/04_seurat_clustered.rds")
m <- gene_mapping(); expr <- data_of(object)
panels <- list(
  T_cells = c("CD3D", "CD3E", "TRAC", "CD2", "IL7R"),
  NK_cells = c("NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1"),
  B_cells = c("MS4A1", "CD79A", "CD79B", "CD37", "CD22"),
  Plasma_cells = c("MZB1", "JCHAIN", "SDC1", "DERL3", "XBP1"),
  Myeloid = c("LYZ", "LST1", "TYROBP", "FCER1G", "AIF1", "C1QA", "CD68"),
  Dendritic = c("FCER1A", "CD1C", "CLEC10A", "CLEC9A", "LILRA4"),
  Mast_cells = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "CDH1"),
  Fibroblasts = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"),
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "CLDN5"),
  Pericytes = c("RGS5", "CSPG4", "MCAM", "PDGFRB", "ACTA2"))
panel_table <- do.call(rbind, lapply(names(panels), function(t) data.frame(cell_type = t, gene_symbol = panels[[t]])))
write.csv(panel_table, "results/annotation/canonical_marker_panel.csv", row.names = FALSE)
clusters <- sort(unique(as.character(object$cluster)))
group <- factor(object$cluster, levels = clusters)
design <- Matrix::sparseMatrix(i = seq_along(group), j = as.integer(group), x = 1,
                              dims = c(length(group), length(clusters)))
n <- as.numeric(table(group)); means <- expr %*% design %*% Matrix::Diagonal(x = 1/n)
all_symbols <- unique(unlist(panels)); ids <- m$seurat_feature[match(all_symbols, m$gene_symbol)]
ok <- !is.na(ids) & ids %in% rownames(expr); ids <- ids[ok]; all_symbols <- all_symbols[ok]
avg <- as.matrix(means[ids, , drop = FALSE]); rownames(avg) <- all_symbols; colnames(avg) <- clusters
scaled <- t(scale(t(avg))); scaled[!is.finite(scaled)] <- 0
scores <- sapply(panels, function(genes) colMeans(scaled[intersect(genes, rownames(scaled)), , drop = FALSE]))
write.csv(data.frame(cluster = clusters, scores), "results/annotation/cluster_panel_scores.csv", row.names = FALSE)
write.csv(data.frame(gene_symbol = rownames(avg), avg, check.names = FALSE), "results/annotation/cluster_marker_means.csv", row.names = FALSE)
markers <- readRDS("data/processed/04_cluster_markers.rds")
tops <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 12, with_ties = FALSE)
write.csv(tops, "results/annotation/top12_cluster_markers.csv", row.names = FALSE)
info <- lapply(clusters, function(cl) {
  sel <- object$cluster == cl; tab <- table(object$sample_id[sel]); sc <- scores[cl, ]; ord <- order(sc, decreasing = TRUE)
  data.frame(cluster = cl, n_cells = sum(sel), n_samples = sum(tab > 0),
             dominant_sample = names(which.max(tab)), dominant_sample_fraction = max(tab)/sum(tab),
             candidate = names(sc)[ord[1]], score_margin = sc[ord[1]] - sc[ord[2]],
             median_mt = median(object$percent.mt[sel]))
})
write.csv(do.call(rbind, info), "results/annotation/cluster_annotation_candidates.csv", row.names = FALSE)
# 用明确 marker 显示表达比例和平均表达；基因标签显示 symbol，对象内仍使用稳定 ID。
dot_ids <- ids; names(dot_ids) <- all_symbols
p <- DotPlot(object, features = unname(dot_ids), group.by = "cluster") +
  scale_x_discrete(labels = setNames(all_symbols, ids)) + RotatedAxis() +
  labs(title = "Canonical marker evidence | cluster annotation", x = NULL, y = "Cluster")
save_plot(p, "figures/annotation/05_marker_dotplot", 18, 9)
if ("--apply" %in% commandArgs(trailingOnly = TRUE)) {
  ann <- read.csv("results/annotation/cluster_annotations.csv")
  stopifnot(!anyDuplicated(ann$cluster), setequal(as.character(ann$cluster), clusters),
            all(nzchar(ann$cell_type)), all(nzchar(ann$evidence)))
  object$cell_type <- ann$cell_type[match(object$cluster, ann$cluster)]
  object$cell_subtype <- ann$cell_subtype[match(object$cluster, ann$cluster)]
  object$annotation_confidence <- ann$confidence[match(object$cluster, ann$cluster)]
  object$annotation_cluster <- as.character(object$cluster)
  # 对混合谱系的 2/10 簇使用独立再聚类的已审阅表；完整覆盖、细胞身份和父簇必须一致。
  cells <- read.csv("results/annotation/lymphocyte_refinement_cells.csv")
  refined <- read.csv("results/annotation/lymphocyte_subcluster_annotations.csv")
  stopifnot(!anyDuplicated(cells$cell_id), !anyDuplicated(refined$sub_key),
            setequal(cells$cell_id, colnames(object)[object$cluster %in% c("2", "10")]),
            setequal(cells$sub_key, refined$sub_key))
  idx <- match(cells$cell_id, colnames(object)); j <- match(cells$sub_key, refined$sub_key)
  stopifnot(all(as.character(object$cluster[idx]) == as.character(cells$parent_cluster)))
  object$cell_type[idx] <- refined$cell_type[j]
  object$cell_subtype[idx] <- refined$cell_subtype[j]
  object$annotation_confidence[idx] <- refined$confidence[j]
  object$annotation_cluster[idx] <- cells$sub_key
  stopifnot(!anyNA(object$cell_type))
  saveRDS(object, "data/processed/05_seurat_annotated.rds")
  write.csv(object[[]], "results/annotation/cell_annotations.csv")
  save_plot(DimPlot(object, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE, raster = FALSE, pt.size = .1),
            "figures/annotation/05_umap_celltypes", 12, 8)
  save_plot(DimPlot(object, reduction = "umap_raw", group.by = "cell_type", raster = FALSE, pt.size = .1) +
              DimPlot(object, reduction = "umap", group.by = "cell_type", raster = FALSE, pt.size = .1),
            "figures/annotation/05_celltypes_before_after_harmony", 16, 6)
}
stage_end()
