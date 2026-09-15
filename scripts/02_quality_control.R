# 阶段 02：基础细胞质量控制（QC）；不进行归一化、聚类或细胞注释。
# 输入：阶段 01 合并 Seurat RDS、基因 ID/symbol 映射；应用模式另读阈值 CSV。
# 输出：逐细胞 QC 指标、分位数、阈值敏感性表、QC 图及基础 QC 后对象。
# 项目根目录运行：Rscript --vanilla scripts/02_quality_control.R
# 首次默认 inspect 模式，只计算指标/绘图，不筛选；查看图表后填写 qc_thresholds.csv。
# 应用阈值：Rscript --vanilla scripts/02_quality_control.R --apply
# RStudio 中 source 默认 inspect；可设置环境变量 CRC_QC_MODE=apply 后 source。

main <- function() {
  if (!file.exists("scRNA_CRC_project.Rproj")) stop("请从含 .Rproj 的项目根目录运行。")
  set.seed(20260915)
  if (dir.exists("renv/library-local")) .libPaths(c("renv/library-local", .libPaths()))
  for (pkg in c("Matrix", "SeuratObject", "ggplot2", "patchwork")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("缺少依赖：", pkg)
  }
  suppressPackageStartupMessages(library(ggplot2))
  dirs <- c("results/qc", "figures/qc", "data/processed", "logs")
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  apply_mode <- "--apply" %in% commandArgs(trailingOnly = TRUE) ||
    identical(Sys.getenv("CRC_QC_MODE"), "apply")
  mode <- if (apply_mode) "apply" else "inspect"
  logfile <- paste0("logs/02_qc_", mode, ".log")
  cat("", file = logfile)
  log_msg <- function(...) {
    line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " "))
    cat(line, "\n"); cat(line, "\n", file = logfile, append = TRUE)
  }
  on.exit(writeLines(capture.output(sessionInfo()), paste0("logs/02_qc_", mode, "_sessionInfo.txt")))
  log_msg("运行模式：", mode)

  # 1. 输入：初始对象、基因映射。输出：从 counts 独立计算的每细胞 QC 指标。
  object <- readRDS("data/processed/01_seurat_merged.rds")
  map <- readRDS("data/processed/01_gene_feature_mapping.rds")
  # ID 为行名，不能直接对 rownames(object) 搜索 ^MT-；先用 symbol 识别再映射到 ID。
  gene_map <- unique(map[, c("gene_id", "seurat_feature", "gene_symbol")])
  if (anyDuplicated(gene_map$seurat_feature)) stop("同一 feature 存在冲突映射，需要人工检查。")
  mt <- gene_map[grepl("^MT-", gene_map$gene_symbol), ]
  if (!nrow(mt) || !all(mt$seurat_feature %in% rownames(object))) stop("线粒体基因映射失败。")
  write.csv(mt, "results/qc/mitochondrial_genes.csv", row.names = FALSE)
  # v4 从 counts slot 读取；v5 逐 counts layer 读取，避免不必要的稠密化。
  if (inherits(object[["RNA"]], "Assay5")) {
    layers <- SeuratObject::Layers(object[["RNA"]], search = "^counts")
    get_counts <- function(layer) SeuratObject::LayerData(object, assay = "RNA", layer = layer)
  } else {
    layers <- "counts"
    get_counts <- function(layer) SeuratObject::GetAssayData(object, assay = "RNA", slot = "counts")
  }
  metrics <- lapply(layers, function(layer) {
    x <- get_counts(layer)
    umi <- Matrix::colSums(x)
    detected <- Matrix::colSums(x > 0)
    mt_umi <- Matrix::colSums(x[intersect(mt$seurat_feature, rownames(x)), , drop = FALSE])
    data.frame(cell_id = colnames(x), counts_recomputed = umi,
               features_recomputed = detected, mt_counts = mt_umi,
               percent.mt = ifelse(umi > 0, 100 * mt_umi / umi, NA_real_))
  })
  metrics <- do.call(rbind, metrics)
  stopifnot(!anyDuplicated(metrics$cell_id), setequal(metrics$cell_id, colnames(object)))
  metrics <- metrics[match(colnames(object), metrics$cell_id), ]
  qc <- cbind(data.frame(cell_id = colnames(object)), object[[]],
              metrics[, c("mt_counts", "percent.mt")])
  stopifnot(all(metrics$counts_recomputed == qc$nCount_RNA),
            all(metrics$features_recomputed == qc$nFeature_RNA),
            all(is.finite(qc$percent.mt)), all(qc$percent.mt >= 0 & qc$percent.mt <= 100))
  sample_order <- unique(qc$sample_id)
  qc$sample_id <- factor(qc$sample_id, levels = sample_order)
  # 高 counts/高基因数仅作异常提示：肿瘤或大细胞也可能很高，不能直接视为 doublet。
  qc$high_counts_flag <- qc$high_features_flag <- FALSE
  diagnostics <- list()
  for (s in sample_order) {
    k <- qc$sample_id == s
    upper <- function(x) 10^(median(log10(x + 1)) + 3 * mad(log10(x + 1))) - 1
    count_high <- upper(qc$nCount_RNA[k]); feature_high <- upper(qc$nFeature_RNA[k])
    qc$high_counts_flag[k] <- qc$nCount_RNA[k] > count_high
    qc$high_features_flag[k] <- qc$nFeature_RNA[k] > feature_high
    diagnostics[[s]] <- data.frame(sample_id = s, counts_upper_3mad = count_high,
                                   features_upper_3mad = feature_high,
                                   mt_upper_3mad = median(qc$percent.mt[k]) + 3 * mad(qc$percent.mt[k]))
  }
  write.csv(do.call(rbind, diagnostics), "results/qc/mad_diagnostics.csv", row.names = FALSE)
  write.csv(qc, "results/qc/cell_qc_metrics.csv", row.names = FALSE)
  saveRDS(qc, "data/processed/02_qc_metrics.rds")

  # 2. 输入：QC 指标。输出：分样本分位数和不同候选阈值下的细胞保留率。
  probs <- c(0, .01, .05, .25, .5, .75, .95, .99, 1)
  stats <- list(); sensitivity <- list()
  for (s in sample_order) {
    d <- qc[qc$sample_id == s, ]
    for (metric in c("nCount_RNA", "nFeature_RNA", "percent.mt")) {
      stats[[length(stats) + 1L]] <- data.frame(sample_id = s, metric = metric,
                                               quantile = probs, value = as.numeric(quantile(d[[metric]], probs)))
    }
    grid <- expand.grid(min_features = c(200, 300, 500), max_mt = c(10, 15, 20, 25, 30))
    for (j in seq_len(nrow(grid))) {
      keep <- d$nFeature_RNA >= grid$min_features[j] & d$percent.mt <= grid$max_mt[j]
      sensitivity[[length(sensitivity) + 1L]] <- data.frame(
        sample_id = s, grid[j, ], n_before = nrow(d), n_after = sum(keep), retained_pct = 100 * mean(keep))
    }
  }
  write.csv(do.call(rbind, stats), "results/qc/qc_quantiles.csv", row.names = FALSE)
  write.csv(do.call(rbind, sensitivity), "results/qc/threshold_sensitivity.csv", row.names = FALSE)

  # 3. 输入：全部细胞 QC 指标。输出：PNG（查看）+ PDF（矢量导出）。
  # 图内使用英文以避免机器中文字体差异；注释、解读和参数说明使用中文。
  palette <- c("left-sided" = "#287D8E", "right-sided" = "#CC7546")
  theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                             plot.title = element_text(face = "bold")))
  save_plot <- function(plot, stem, width, height) {
    # 显式使用 R 自带设备，避免自动选中与本机 R 图形 API 不兼容的 ragg。
    ggsave(paste0("figures/qc/", stem, ".png"), plot, device = grDevices::png,
           width = width, height = height, dpi = 160)
    ggsave(paste0("figures/qc/", stem, ".pdf"), plot, device = grDevices::pdf,
           width = width, height = height)
  }
  violin <- function(data, metric, label) {
    p <- ggplot(data, aes(x = sample_id, y = .data[[metric]], fill = side)) +
      geom_violin(scale = "width", linewidth = .25) +
      geom_boxplot(width = .13, outlier.shape = NA, fill = "white", linewidth = .3) +
      scale_fill_manual(values = palette) + labs(x = NULL, y = label, fill = "Tumor side")
    if (metric != "percent.mt") p <- p + scale_y_log10()
    p
  }
  make_panel <- function(data, title) {
    patchwork::wrap_plots(violin(data, "nCount_RNA", "UMI counts (log10 axis)"),
                         violin(data, "nFeature_RNA", "Detected genes (log10 axis)"),
                         violin(data, "percent.mt", "Mitochondrial UMI (%)"), ncol = 3, guides = "collect") +
      patchwork::plot_annotation(title = title)
  }
  save_plot(make_panel(qc, paste0("GSE188711 | Before QC | ", nrow(qc), " cells")), "02_qc_before", 13, 4.6)
  # 所有细胞参与绘图，无随机抽样；透明度降低重叠点遮挡。
  scatter <- ggplot(qc, aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
    geom_point(size = .3, alpha = .35) + scale_x_log10() + scale_y_log10() +
    scale_color_viridis_c(limits = c(0, 100)) + facet_wrap(~sample_id, ncol = 3) +
    labs(x = "UMI counts (log10 axis)", y = "Detected genes (log10 axis)", color = "MT UMI (%)",
         title = "Counts, gene complexity and mitochondrial fraction | Before QC")
  save_plot(scatter, "02_qc_scatter", 11, 7)
  log_msg("指标检查通过：", nrow(qc), "细胞；线粒体基因数：", nrow(mt))
  if (!apply_mode) {
    log_msg("inspect 完成。先查看分布图与敏感性表，再填写阈值并使用 --apply。")
    return(invisible(qc))
  }

  # 4. 输入：明确记录的按样本阈值表。输出：逐细胞决定、过滤统计、基础 QC 后对象。
  # 参数 CSV 必须每样本一行，可统一或按样本设置；不得按左右侧直接设不同规则。
  thresholds <- read.csv("results/qc/qc_thresholds.csv", stringsAsFactors = FALSE)
  required_columns <- c("sample_id", "min_features", "min_counts", "max_mt")
  stopifnot(all(required_columns %in% names(thresholds)),
            !anyDuplicated(thresholds$sample_id), setequal(thresholds$sample_id, sample_order))
  stopifnot(all(is.finite(as.matrix(thresholds[, required_columns[-1]]))),
            all(thresholds$min_features >= 0), all(thresholds$min_counts >= 0),
            all(thresholds$max_mt >= 0 & thresholds$max_mt <= 100))
  t <- thresholds[match(qc$sample_id, thresholds$sample_id), ]
  qc$fail_low_features <- qc$nFeature_RNA < t$min_features
  qc$fail_low_counts <- qc$nCount_RNA < t$min_counts
  qc$fail_high_mt <- qc$percent.mt > t$max_mt
  qc$pass_basic_qc <- !(qc$fail_low_features | qc$fail_low_counts | qc$fail_high_mt)
  # 原因列允许重叠；n_removed 用联合条件计算，不可把各原因数量相加。
  decisions <- do.call(rbind, lapply(sample_order, function(s) {
    d <- qc[qc$sample_id == s, ]
    data.frame(sample_id = s, side = d$side[1], n_before = nrow(d), n_after = sum(d$pass_basic_qc),
               n_removed = sum(!d$pass_basic_qc), retained_pct = 100 * mean(d$pass_basic_qc),
               fail_low_features = sum(d$fail_low_features), fail_low_counts = sum(d$fail_low_counts),
               fail_high_mt = sum(d$fail_high_mt),
               high_counts_flag_retained = sum(d$high_counts_flag & d$pass_basic_qc),
               high_features_flag_retained = sum(d$high_features_flag & d$pass_basic_qc))
  }))
  if (any(decisions$n_after == 0)) stop("至少一个样本全部被删除，停止保存过滤对象。")
  rownames(qc) <- qc$cell_id
  object <- SeuratObject::AddMetaData(object, qc[, c("mt_counts", "percent.mt", "high_counts_flag",
                        "high_features_flag", "fail_low_features", "fail_low_counts", "fail_high_mt", "pass_basic_qc")])
  filtered <- subset(object, cells = qc$cell_id[qc$pass_basic_qc])
  stopifnot(ncol(filtered) == sum(qc$pass_basic_qc), all(filtered$pass_basic_qc),
            identical(colnames(filtered), qc$cell_id[qc$pass_basic_qc]))
  filtered@misc$stage02 <- list(seed = 20260915L, thresholds = thresholds,
    doublet_detection = "not yet performed; high-count flags are not doublet calls",
    gene_filtering = FALSE, normalized = FALSE, integrated = FALSE, annotated = FALSE,
    mitochondrial_features = mt$seurat_feature)
  # 保留阶段 01 原件；未过滤完整细胞指标/决策另存，避免复制一个大型未过滤对象。
  saveRDS(qc, "data/processed/02_qc_decisions.rds")
  saveRDS(thresholds, "data/processed/02_qc_thresholds.rds")
  write.csv(qc, "results/qc/cell_qc_decisions.csv", row.names = FALSE)
  write.csv(decisions, "results/qc/qc_filtering_summary.csv", row.names = FALSE)
  saveRDS(filtered, "data/processed/02_seurat_basic_qc.rds")
  restored <- readRDS("data/processed/02_seurat_basic_qc.rds")
  stopifnot(inherits(restored, "Seurat"), identical(dim(restored), dim(filtered)),
            identical(restored[[]], filtered[[]]), isTRUE(methods::validObject(restored)))
  save_plot(make_panel(qc[qc$pass_basic_qc, ], paste0("GSE188711 | After basic QC | ", ncol(filtered), " cells")),
            "02_qc_after", 13, 4.6)
  retention <- ggplot(decisions, aes(factor(sample_id, levels = sample_order), retained_pct, fill = side)) +
    geom_col(width = .65) + geom_text(aes(label = sprintf("%s / %s\n%.1f%%", n_after, n_before, retained_pct)),
                                    vjust = -.3, size = 3.3) +
    scale_fill_manual(values = palette) + scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 20)) +
    labs(x = NULL, y = "Cells retained (%)", fill = "Tumor side", title = "Basic QC retention by sample")
  save_plot(retention, "02_qc_retention", 8, 5)
  log_msg("基础 QC 完成：", ncol(object), "→", ncol(filtered), "细胞；RDS 回读通过。")
  print(decisions, row.names = FALSE)
}

main()
