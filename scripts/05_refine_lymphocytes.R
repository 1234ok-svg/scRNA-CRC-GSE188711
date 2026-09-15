# 输入：主聚类对象；输出：细胞毒性簇 2 和增殖簇 10 的局部再聚类及 marker。
# 目的：区分 T/NK 身份、B/T 谱系与增殖状态；不依赖左右侧标签。
# 本步骤针对已审阅的固定主聚类结果；如上游聚类变化，应重新审阅 parent cluster。
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source("scripts/00_common.R", encoding = "UTF-8")
suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
stage_start("05_refine_lymphocytes")
ensure_dirs("results/annotation", "figures/annotation")
object <- readRDS("data/processed/04_seurat_clustered.rds"); m <- gene_mapping()
cc <- unique(c(Seurat::cc.genes.updated.2019$s.genes, Seurat::cc.genes.updated.2019$g2m.genes))
exclude <- m$seurat_feature[m$gene_symbol %in% cc | grepl("^MT-", m$gene_symbol)]
cells <- list(); marker_tables <- list(); diagnostics <- list()
panel <- c("CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "NKG7", "GNLY", "KLRD1", "KLRF1",
           "MS4A1", "CD79A", "CD79B", "CD37", "MZB1", "JCHAIN", "MKI67", "TOP2A", "FOXP3")
ids <- m$seurat_feature[match(panel, m$gene_symbol)]; names(ids) <- panel
ids <- ids[!is.na(ids)]
for (parent in c("2", "10")) {
  set.seed(20260915)
  sub <- subset(object, cells = colnames(object)[object$cluster == parent])
  sub <- FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)
  VariableFeatures(sub) <- setdiff(VariableFeatures(sub), exclude)
  sub <- ScaleData(sub, features = VariableFeatures(sub), verbose = FALSE)
  sub <- RunPCA(sub, npcs = 20, verbose = FALSE, seed.use = 20260915)
  sub <- FindNeighbors(sub, dims = 1:15, verbose = FALSE)
  sub <- FindClusters(sub, resolution = .3, random.seed = 20260915, verbose = FALSE)
  sub$sub_key <- paste(parent, as.character(Idents(sub)), sep = "_")
  Idents(sub) <- "sub_key"
  cells[[parent]] <- data.frame(cell_id = colnames(sub), parent_cluster = parent, sub_key = sub$sub_key)
  mk <- FindAllMarkers(sub, only.pos = TRUE, min.pct = .15, logfc.threshold = .25,
                       max.cells.per.ident = 400, random.seed = 20260915, verbose = FALSE)
  mk$gene_symbol <- m$gene_symbol[match(mk$gene, m$seurat_feature)]
  marker_tables[[parent]] <- mk
  p <- DotPlot(sub, features = unname(ids), group.by = "sub_key") +
    scale_x_discrete(labels = setNames(names(ids), ids)) + RotatedAxis() +
    labs(title = paste("Lineage evidence in parent cluster", parent), x = NULL, y = "Subcluster")
  save_plot(p, paste0("figures/annotation/05_refinement_", parent), 12, 4)
  avg <- AverageExpression(sub, features = unname(ids), group.by = "sub_key", verbose = FALSE)$RNA
  diagnostics[[parent]] <- data.frame(gene_symbol = names(ids)[match(rownames(avg), ids)], as.matrix(avg), check.names = FALSE)
  write.csv(diagnostics[[parent]], paste0("results/annotation/refinement_", parent, "_expression.csv"), row.names = FALSE)
}
write.csv(do.call(rbind, cells), "results/annotation/lymphocyte_refinement_cells.csv", row.names = FALSE)
write.csv(do.call(rbind, marker_tables), "results/annotation/lymphocyte_refinement_markers.csv", row.names = FALSE)
saveRDS(list(cells = cells, markers = marker_tables, expression = diagnostics), "data/processed/05_lymphocyte_refinement.rds")
stage_end()
