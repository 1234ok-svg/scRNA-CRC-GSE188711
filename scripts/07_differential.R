# 输入：注释对象的原始 counts；输出：患者×细胞类型 pseudobulk、edgeR QL 结果、敏感性和图。
# 对比方向统一为 right-sided / left-sided；细胞不是独立生物学重复。
# 每类型必须六个患者都至少有 30 个细胞，否则明确跳过，不补造生物学重复。
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source("scripts/00_common.R", encoding = "UTF-8")
suppressPackageStartupMessages({library(SeuratObject); library(edgeR); library(ggplot2)})
stage_start("07_differential")
ensure_dirs("results/differential", "figures/differential")
object <- readRDS("data/processed/05_seurat_annotated.rds")
counts <- counts_of(object); meta <- object[[]]; sm <- read.csv("data/metadata/sample_metadata.csv")
m <- gene_mapping(); all_results <- list(); statuses <- list(); pseudobulks <- list(); mds <- list(); stability <- list()
types <- setdiff(sort(unique(meta$cell_type)), c("Unknown", "Mixed", "Unresolved"))
for (scenario in c("primary_mt25", "sensitivity_mt20")) {
  for (ct in types) {
    key <- paste(scenario, ct, sep = "__")
    sel <- meta$cell_type == ct & (scenario == "primary_mt25" | meta$percent.mt <= 20)
    n <- table(factor(meta$sample_id[sel], levels = sm$sample_id))
    status <- data.frame(scenario = scenario, cell_type = ct, min_cells_per_patient = min(n),
                          status = if (any(n < 30)) "SKIPPED_less_than_30_cells_in_at_least_one_patient" else "OK")
    if (any(n < 30)) { statuses[[key]] <- status; next }
    # 每个患者一列，直接相加 UMI，不能使用均值或 Harmony 数据。
    pb <- sapply(sm$sample_id, function(s) Matrix::rowSums(counts[, sel & meta$sample_id == s, drop = FALSE]))
    colnames(pb) <- sm$sample_id
    group <- factor(sm$side, levels = c("left-sided", "right-sided"))
    design <- model.matrix(~group)
    y <- DGEList(pb, group = group); keep <- filterByExpr(y, design = design)
    y <- y[keep, , keep.lib.sizes = FALSE]
    if (nrow(y) < 100) { status$status <- "SKIPPED_less_than_100_expressed_genes"; statuses[[key]] <- status; next }
    y <- calcNormFactors(y, method = "TMM"); y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design, robust = TRUE); test <- glmQLFTest(fit, coef = 2)
    result <- topTags(test, n = Inf, sort.by = "none")$table
    result$gene_id <- rownames(result); result$gene_symbol <- m$gene_symbol[match(result$gene_id, m$seurat_feature)]
    result$cell_type <- ct; result$scenario <- scenario
    result$significant <- result$FDR < .05 & abs(result$logFC) >= 1
    result$direction <- ifelse(result$logFC > 0, "higher_right", "higher_left")
    all_results[[key]] <- result; pseudobulks[[key]] <- list(counts = pb, cells_per_patient = n,
                                                           tested_genes = rownames(y), design = design,
                                                           norm_factors = y$samples, fit = fit)
    status$n_tested <- nrow(result); status$n_significant <- sum(result$significant)
    statuses[[key]] <- status
    log_msg(key, "tested", nrow(result), "DE genes", sum(result$significant))
    if (scenario == "primary_mt25") {
      lcpm <- cpm(y, log = TRUE, prior.count = 2)
      top <- order(apply(lcpm, 1, var), decreasing = TRUE)[seq_len(min(500, nrow(lcpm)))]
      pc <- prcomp(t(lcpm[top, ]), scale. = FALSE)$x[, 1:2]
      mds[[ct]] <- data.frame(sample_id = rownames(pc), PC1 = pc[, 1], PC2 = pc[, 2], cell_type = ct, side = sm$side)
      # 描述性 leave-one-patient-out 方向检查；不是重新拟合的显著性检验。
      eff <- sapply(seq_len(ncol(lcpm)), function(drop) {
        l <- which(group == "left-sided" & seq_along(group) != drop)
        r <- which(group == "right-sided" & seq_along(group) != drop)
        rowMeans(lcpm[, r, drop = FALSE]) - rowMeans(lcpm[, l, drop = FALSE])
      })
      ref <- result$logFC[match(rownames(lcpm), result$gene_id)]
      stability[[ct]] <- data.frame(gene_id = rownames(lcpm), cell_type = ct,
                                     loo_direction_consistency = rowMeans(sign(eff) == sign(ref)),
                                     min_loo_logcpm_difference = apply(eff, 1, min),
                                     max_loo_logcpm_difference = apply(eff, 1, max))
    }
  }
}
# 不同状态行补齐列后再拼接，保留被跳过的分析和具体原因。
cols <- unique(unlist(lapply(statuses, names)))
status_table <- do.call(rbind, lapply(statuses, function(x) {for (z in setdiff(cols, names(x))) x[[z]] <- NA; x[, cols]}))
write.csv(status_table, "results/differential/analysis_status.csv", row.names = FALSE)
if (!length(all_results)) stop("没有细胞类型满足 pseudobulk 样本门槛，见状态表。")
res <- do.call(rbind, all_results)
# 提供每类型 FDR 以及跨主分析全部 gene×type 假设的 FDR，防止选择性报告。
res$FDR_global <- NA_real_
for (sc in unique(res$scenario)) {
  idx <- res$scenario == sc; res$FDR_global[idx] <- p.adjust(res$PValue[idx], "BH")
}
write.csv(res, "results/differential/pseudobulk_all_results.csv", row.names = FALSE)
write.csv(res[res$significant, ], "results/differential/pseudobulk_significant.csv", row.names = FALSE)
saveRDS(pseudobulks, "data/processed/07_pseudobulk_models.rds")
saveRDS(res, "data/processed/07_differential_results.rds")
loo <- do.call(rbind, stability)
write.csv(loo, "results/differential/leave_one_patient_out_direction.csv", row.names = FALSE)
primary <- res[res$scenario == "primary_mt25", ]; sensitive <- res[res$scenario == "sensitivity_mt20", ]
comparison <- merge(primary[, c("gene_id", "cell_type", "logFC", "FDR", "significant")],
                    sensitive[, c("gene_id", "cell_type", "logFC", "FDR", "significant")],
                    by = c("gene_id", "cell_type"), suffixes = c("_mt25", "_mt20"))
write.csv(comparison, "results/differential/mt20_comparison.csv", row.names = FALSE)
p <- ggplot(primary, aes(logFC, -log10(pmax(FDR, 1e-300)), color = significant)) +
  geom_point(size = .5, alpha = .5) + facet_wrap(~cell_type, ncol = 3, scales = "free_y") +
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey60") +
  geom_hline(yintercept = -log10(.05), linetype = 2, color = "grey60") +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#B54D42")) + theme_bw() +
  labs(title = "Patient-level pseudobulk | right versus left", x = "log2 fold change (right / left)", y = "-log10 within-type FDR")
save_plot(p, "figures/differential/07_volcano", 12, max(5, ceiling(length(unique(primary$cell_type))/3)*3))
if (length(mds)) {
  p <- ggplot(do.call(rbind, mds), aes(PC1, PC2, color = side, label = sample_id)) + geom_point(size = 3) +
    geom_text(vjust = -1, size = 3) + facet_wrap(~cell_type, scales = "free", ncol = 3) +
    scale_color_manual(values = side_colors) + theme_bw() + labs(title = "Pseudobulk PCA: each point is one patient")
  save_plot(p, "figures/differential/07_pseudobulk_pca", 12, 8)
}
stage_end()
