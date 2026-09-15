# 通用函数；输入：项目根目录、独立 R 包库；输出：日志、稳定图形设备和计数访问器。
# 被后续脚本 source，不单独运行；不进行网络安装。
stopifnot(file.exists("scRNA_CRC_project.Rproj"))
.libPaths(c("renv/library-local", .libPaths()))
set.seed(20260915)
options(stringsAsFactors = FALSE, future.globals.maxSize = 4 * 1024^3)
Sys.setenv(OMP_NUM_THREADS = "2", OPENBLAS_NUM_THREADS = "2")
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")
dir.create("logs", showWarnings = FALSE)
stage_start <- function(stage) {
  assign(".stage", stage, envir = .GlobalEnv)
  cat("", file = paste0("logs/", stage, ".log"))
  log_msg("开始", stage)
}
log_msg <- function(...) {
  z <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " "))
  cat(z, "\n"); cat(z, "\n", file = paste0("logs/", .stage, ".log"), append = TRUE)
}
stage_end <- function() {
  writeLines(capture.output(sessionInfo()), paste0("logs/", .stage, "_sessionInfo.txt"))
  log_msg("成功完成")
}
ensure_dirs <- function(...) invisible(lapply(c(...), dir.create, recursive = TRUE, showWarnings = FALSE))
save_plot <- function(p, path, w = 9, h = 6) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(path, ".png"), p, device = grDevices::png, width = w, height = h, dpi = 170)
  ggplot2::ggsave(paste0(path, ".pdf"), p, device = grDevices::pdf, width = w, height = h)
}
counts_of <- function(x) {
  if (inherits(x[["RNA"]], "Assay5")) {
    if (length(SeuratObject::Layers(x[["RNA"]], search = "^counts")) > 1)
      x <- SeuratObject::JoinLayers(x, assay = "RNA")
    SeuratObject::LayerData(x, assay = "RNA", layer = "counts")
  } else SeuratObject::GetAssayData(x, assay = "RNA", slot = "counts")
}
data_of <- function(x) {
  if (inherits(x[["RNA"]], "Assay5")) SeuratObject::LayerData(x, assay = "RNA", layer = "data")
  else SeuratObject::GetAssayData(x, assay = "RNA", slot = "data")
}
gene_mapping <- function() {
  m <- readRDS("data/processed/01_gene_feature_mapping.rds")
  unique(m[, c("gene_id", "seurat_feature", "gene_symbol")])
}
side_colors <- c("left-sided" = "#287D8E", "right-sided" = "#CC7546")
