# 输入：注释对象与 pseudobulk 结果；输出：整合诊断、状态 marker、总体复现检查。
# marker 状态评分仅为预先定义基因的平均 log-normalized 表达，非临床分类或功能实验。
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source("scripts/00_common.R", encoding = "UTF-8")
suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
stage_start("09_validation")
ensure_dirs("results/validation", "figures/validation")
object <- readRDS("data/processed/05_seurat_annotated.rds")
md <- object[[]]; m <- gene_mapping(); expr <- data_of(object)
# 同细胞类型邻居比例评估是否保持类型结构；同样本邻居比例不是单独的整合质量标准。
diagnostics <- list()
set.seed(20260915)
selected <- sort(sample(seq_len(ncol(object)), min(5000, ncol(object))))
for (reduction in c("pca", "harmony")) {
  emb <- Embeddings(object, reduction)[, 1:30]
  nn <- RANN::nn2(data = emb, query = emb[selected, , drop = FALSE], k = 31)$nn.idx[, -1, drop = FALSE]
  same_type <- rowMeans(matrix(md$cell_type[nn], nrow = nrow(nn)) == md$cell_type[selected])
  same_sample <- rowMeans(matrix(md$sample_id[nn], nrow = nrow(nn)) == md$sample_id[selected])
  diagnostics[[reduction]] <- data.frame(cell_id = colnames(object)[selected], cell_type = md$cell_type[selected],
                                         reduction = reduction, same_type_fraction = same_type, same_sample_fraction = same_sample)
}
diag <- do.call(rbind, diagnostics)
write.csv(aggregate(cbind(same_type_fraction, same_sample_fraction) ~ reduction + cell_type, diag, mean),
          "results/validation/integration_neighbors.csv", row.names = FALSE)
saveRDS(diag, "data/processed/09_integration_diagnostics.rds")
signatures <- list(Cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH"),
                    Checkpoint_related = c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4", "TOX"),
                    Interferon_response = c("ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "STAT1"),
                    Cycling = c("MKI67", "TOP2A", "UBE2C", "CENPF"))
state <- list()
for (name in names(signatures)) {
  ids <- m$seurat_feature[m$gene_symbol %in% signatures[[name]]]
  score <- Matrix::colMeans(expr[intersect(ids, rownames(expr)), , drop = FALSE])
  state[[name]] <- aggregate(score, list(sample_id = md$sample_id, cell_type = md$cell_type), mean)
  names(state[[name]])[3] <- "mean_signature_expression"
  state[[name]]$signature <- name
}
state <- do.call(rbind, state); sm <- read.csv("data/metadata/sample_metadata.csv")
state$side <- sm$side[match(state$sample_id, sm$sample_id)]
write.csv(state, "results/validation/patient_marker_signatures.csv", row.names = FALSE)
write.csv(do.call(rbind, lapply(names(signatures), function(x) data.frame(signature = x, gene_symbol = signatures[[x]]))),
          "results/validation/signature_genes.csv", row.names = FALSE)
lymph <- state[state$cell_type %in% c("T_cells", "NK_cells", "T_NK"), ]
if (nrow(lymph)) {
  p <- ggplot(lymph, aes(side, mean_signature_expression, color = side, label = sample_id)) +
    geom_point(size = 2.5) + geom_text(nudge_x = .1, size = 2.5) +
    facet_grid(signature ~ cell_type, scales = "free_y") + scale_color_manual(values = side_colors) +
    theme_bw() + theme(legend.position = "none") + labs(title = "Descriptive marker programs | patient means", x = NULL, y = "Mean log-normalized expression")
  save_plot(p, "figures/validation/09_lymphocyte_programs", 10, 10)
}
# 最终对象的实际 counts 与原 QC 元数据相符，所有细胞都为预测 singlet。
x <- counts_of(object)
stopifnot(all(Matrix::colSums(x) == object$nCount_RNA),
            all(Matrix::colSums(x > 0) == object$nFeature_RNA),
            all(object$doublet_class == "singlet"), !anyNA(object$cell_type), !anyDuplicated(colnames(object)))
paths <- c("data/processed/02_seurat_basic_qc.rds", "data/processed/03_seurat_singlets.rds",
           "data/processed/04_seurat_clustered.rds", "data/processed/05_seurat_annotated.rds")
write.csv(data.frame(path = paths, md5 = unname(tools::md5sum(paths))), "results/validation/object_checksums.csv", row.names = FALSE)
writeLines(c("PASS: cell identities are unique", "PASS: all annotated cells are predicted singlets",
             "PASS: actual counts equal metadata", "PASS: annotation complete (including explicit unresolved labels if needed)"),
           "results/validation/final_checks.txt")
stage_end()
