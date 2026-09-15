# 输入：singlet counts；输出：标准化数据、PCA/Harmony/UMAP、聚类对象与诊断图。
# Harmony 仅校正低维表示；组间 DE 始终使用原始 counts。未把 side 作为待消除批次。
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source("scripts/00_common.R", encoding = "UTF-8")
suppressPackageStartupMessages({library(Seurat); library(harmony); library(ggplot2); library(patchwork)})
stage_start("04_clustering")
ensure_dirs("results/clustering", "figures/clustering", "results/markers")
object <- readRDS("data/processed/03_seurat_singlets.rds")
set.seed(20260915)
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
# 剔除线粒体基因作为高变特征，保留其原 counts，不回归掉可能具有生物意义的所有 QC 指标。
m <- gene_mapping(); mt_ids <- m$seurat_feature[grepl("^MT-", m$gene_symbol)]
VariableFeatures(object) <- setdiff(VariableFeatures(object), mt_ids)
object <- ScaleData(object, features = VariableFeatures(object), verbose = FALSE)
object <- RunPCA(object, features = VariableFeatures(object), npcs = 40, seed.use = 20260915, verbose = FALSE)
write.csv(data.frame(feature = VariableFeatures(object)), "results/clustering/variable_features.csv", row.names = FALSE)
save_plot(ElbowPlot(object, ndims = 40), "figures/clustering/04_pca_elbow")
log_msg("PCA 完成，生成未校正 UMAP")
object <- RunUMAP(object, reduction = "pca", dims = 1:30, reduction.name = "umap_raw",
                  seed.use = 20260915, verbose = FALSE)
object <- FindNeighbors(object, reduction = "pca", dims = 1:30, graph.name = c("raw_nn", "raw_snn"), verbose = FALSE)
object <- FindClusters(object, graph.name = "raw_snn", resolution = .5, random.seed = 20260915, verbose = FALSE)
object$cluster_raw <- as.character(Idents(object))
# 采用温和 Harmony（theta=1）；保留 raw 表示作对照，并在后续注释后检查细胞类型保持情况。
log_msg("Harmony 按 sample_id 校正低维表示")
emb <- harmony::RunHarmony(Embeddings(object, "pca")[, 1:30],
                          meta_data = object[[]], vars_use = "sample_id", theta = 1,
                          max_iter = 20, ncores = 2, verbose = TRUE)
object[["harmony"]] <- CreateDimReducObject(embeddings = emb, key = "harmony_", assay = "RNA")
object <- RunUMAP(object, reduction = "harmony", dims = 1:30, reduction.name = "umap",
                  seed.use = 20260915, verbose = FALSE)
object <- FindNeighbors(object, reduction = "harmony", dims = 1:30, verbose = FALSE)
object <- FindClusters(object, resolution = c(.3, .5, .8), random.seed = 20260915, verbose = FALSE)
object$cluster <- as.character(object$RNA_snn_res.0.5)
Idents(object) <- "cluster"
write.csv(object[[]], "results/clustering/cell_clusters.csv")
save_plot((DimPlot(object, reduction = "umap_raw", group.by = "sample_id", raster = FALSE, pt.size = .1) + ggtitle("Before Harmony")) +
            (DimPlot(object, reduction = "umap", group.by = "sample_id", raster = FALSE, pt.size = .1) + ggtitle("After Harmony")),
          "figures/clustering/04_batch_comparison", 13, 5)
save_plot(DimPlot(object, reduction = "umap", group.by = "cluster", label = TRUE, raster = FALSE, pt.size = .1),
          "figures/clustering/04_umap_clusters", 10, 7)
save_plot(DimPlot(object, reduction = "umap", group.by = "side", raster = FALSE, pt.size = .1, cols = side_colors),
          "figures/clustering/04_umap_side", 9, 6)
# 移除可重建的大型 scale.data 降低 RDS 体积；PCA embeddings/loadings 保留。
if (inherits(object[["RNA"]], "Assay5")) {
  object[["RNA"]]$scale.data <- NULL
} else {
  object[["RNA"]]@scale.data <- matrix(numeric(), nrow = 0, ncol = 0)
}
saveRDS(object, "data/processed/04_seurat_clustered.rds")
log_msg("寻找簇 marker（仅用于注释，不能作为左右侧组间检验）")
# 全部细胞参与聚类；每簇最多抽取 400 个细胞做描述性 marker 检验，固定种子。
markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = .20, logfc.threshold = .35,
                          max.cells.per.ident = 400, random.seed = 20260915, verbose = FALSE)
markers$gene_id <- markers$gene
markers$gene_symbol <- m$gene_symbol[match(markers$gene, m$seurat_feature)]
write.csv(markers, "results/markers/cluster_markers.csv", row.names = FALSE)
saveRDS(markers, "data/processed/04_cluster_markers.rds")
stage_end()
