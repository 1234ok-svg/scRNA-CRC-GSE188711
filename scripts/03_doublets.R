# 输入：基础 QC 对象；输出：按样本的双细胞分数/判定、统计图、singlet Seurat RDS。
# 标准 10x 数据按捕获样本分别建模；预测不是实验验证，尤其难识别同类型双细胞。
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source("scripts/00_common.R", encoding = "UTF-8")
suppressPackageStartupMessages({library(SeuratObject); library(scDblFinder); library(SingleCellExperiment); library(ggplot2)})
stage_start("03_doublets")
ensure_dirs("results/qc/doublets", "figures/qc", "data/processed")
object <- readRDS("data/processed/02_seurat_basic_qc.rds")
original <- read.csv("results/data_loading/sample_loading_summary.csv")
input_hash <- unname(tools::md5sum("data/processed/02_seurat_basic_qc.rds"))
scores <- list(); summary <- list()
for (i in seq_along(unique(object$sample_id))) {
  s <- unique(object$sample_id)[i]
  checkpoint <- paste0("results/qc/doublets/", s, ".rds")
  # 用原公开矩阵细胞数给出 10x 标准捕获的预期双细胞率，非已知真实发生率。
  dbr <- 0.008 * original$n_cells[match(s, original$sample_id)] / 1000
  cached <- if (file.exists(checkpoint)) readRDS(checkpoint) else NULL
  if (!is.null(cached) && identical(cached$input_hash, input_hash) && identical(cached$dbr, dbr)) {
    out <- cached$scores; log_msg("复用", s, "双细胞检查点")
  } else {
    set.seed(20260915 + i)
    mat <- counts_of(object)[, object$sample_id == s, drop = FALSE]
    log_msg("检测", s, ncol(mat), "细胞；预期率", round(dbr, 4))
    sce <- SingleCellExperiment(assays = list(counts = mat))
    sce <- scDblFinder(sce, clusters = FALSE, dbr = dbr, dbr.sd = .015,
                       nfeatures = 2000, dims = 20, BPPARAM = BiocParallel::SerialParam())
    out <- data.frame(cell_id = colnames(sce), sample_id = s,
                      doublet_score = sce$scDblFinder.score,
                      doublet_class = as.character(sce$scDblFinder.class))
    saveRDS(list(input_hash = input_hash, dbr = dbr, scores = out), checkpoint)
  }
  scores[[s]] <- out
  summary[[s]] <- data.frame(sample_id = s, n_before = nrow(out), expected_rate = dbr,
                              n_doublet = sum(out$doublet_class == "doublet"),
                              n_singlet = sum(out$doublet_class == "singlet"))
}
scores <- do.call(rbind, scores); summary <- do.call(rbind, summary)
stopifnot(!anyDuplicated(scores$cell_id), setequal(scores$cell_id, colnames(object)))
scores <- scores[match(colnames(object), scores$cell_id), ]
object$doublet_score <- scores$doublet_score; object$doublet_class <- scores$doublet_class
singlets <- subset(object, cells = scores$cell_id[scores$doublet_class == "singlet"])
stopifnot(all(singlets$doublet_class == "singlet"), ncol(singlets) == sum(summary$n_singlet))
saveRDS(singlets, "data/processed/03_seurat_singlets.rds")
saveRDS(scores, "data/processed/03_doublet_scores.rds")
write.csv(scores, "results/qc/doublet_scores.csv", row.names = FALSE)
write.csv(summary, "results/qc/doublet_summary.csv", row.names = FALSE)
p <- ggplot(scores, aes(doublet_score, fill = doublet_class)) + geom_histogram(bins = 40, position = "identity", alpha = .7) +
  facet_wrap(~sample_id, scales = "free_y") + theme_bw() +
  scale_fill_manual(values = c("singlet" = "#287D8E", "doublet" = "#CC7546")) +
  labs(title = "scDblFinder predictions by sample", x = "Doublet score", y = "Cells", fill = "Prediction")
save_plot(p, "figures/qc/03_doublet_scores", 10, 6)
log_msg("保留 singlets：", ncol(singlets)); print(summary); stage_end()
