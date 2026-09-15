# 输入：已保存聚类对象；输出：重新绘制降维诊断图，不重跑计算。
# 当对象同时有 umap_raw/umap 时，所有 DimPlot 都显式指定 reduction，避免默认选择 raw。
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source("scripts/00_common.R", encoding = "UTF-8")
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(patchwork)})
stage_start("04_plot_diagnostics")
object <- readRDS("data/processed/04_seurat_clustered.rds")
p1 <- DimPlot(object, reduction = "umap_raw", group.by = "sample_id", raster = FALSE, pt.size = .1) + ggtitle("Before Harmony")
p2 <- DimPlot(object, reduction = "umap", group.by = "sample_id", raster = FALSE, pt.size = .1) + ggtitle("After Harmony")
save_plot(p1 + p2, "figures/clustering/04_batch_comparison", 13, 5)
save_plot(DimPlot(object, reduction = "umap", group.by = "cluster", label = TRUE, raster = FALSE, pt.size = .1),
          "figures/clustering/04_umap_clusters", 10, 7)
save_plot(DimPlot(object, reduction = "umap", group.by = "side", raster = FALSE, pt.size = .1, cols = side_colors),
          "figures/clustering/04_umap_side", 9, 6)
stage_end()
