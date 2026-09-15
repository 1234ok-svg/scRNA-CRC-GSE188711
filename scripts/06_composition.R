# 输入：注释后的 singlets；输出：样本×细胞类型组成、精确置换检验、阈值敏感性。
# 统计重复是 6 位患者。3 vs 3 两侧精确置换仅 20 种分组，P 值分辨率有限。
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source("scripts/00_common.R", encoding = "UTF-8")
suppressPackageStartupMessages({library(SeuratObject); library(ggplot2)})
stage_start("06_composition")
ensure_dirs("results/composition", "figures/composition")
object <- readRDS("data/processed/05_seurat_annotated.rds")
metadata <- object[[]]; sample_meta <- read.csv("data/metadata/sample_metadata.csv")
composition <- function(md, label) {
  tab <- table(factor(md$sample_id, levels = sample_meta$sample_id), factor(md$cell_type, levels = sort(unique(metadata$cell_type))))
  props <- prop.table(tab, 1)
  out <- as.data.frame(tab); names(out) <- c("sample_id", "cell_type", "n_cells")
  out$proportion <- as.vector(props); out$side <- sample_meta$side[match(out$sample_id, sample_meta$sample_id)]
  out$scenario <- label
  tests <- lapply(colnames(props), function(ct) {
    y <- props[, ct]; left <- which(sample_meta$side == "left-sided"); right <- setdiff(1:6, left)
    observed <- mean(y[right]) - mean(y[left])
    permutations <- combn(6, 3, FUN = function(idx) mean(y[idx]) - mean(y[-idx]))
    data.frame(cell_type = ct, mean_left = mean(y[left]), mean_right = mean(y[right]),
               difference_right_minus_left = observed,
               p_value = mean(abs(permutations) >= abs(observed) - 1e-12), scenario = label)
  })
  tests <- do.call(rbind, tests); tests$FDR <- p.adjust(tests$p_value, "BH")
  list(composition = out, tests = tests)
}
primary <- composition(metadata, "mt25_primary")
sensitive <- composition(metadata[metadata$percent.mt <= 20, ], "mt20_sensitivity")
write.csv(primary$composition, "results/composition/sample_composition.csv", row.names = FALSE)
write.csv(primary$tests, "results/composition/composition_tests.csv", row.names = FALSE)
write.csv(rbind(primary$tests, sensitive$tests), "results/composition/mt20_sensitivity.csv", row.names = FALSE)
saveRDS(list(primary = primary, sensitivity = sensitive), "data/processed/06_composition.rds")
p <- ggplot(primary$composition, aes(sample_id, proportion, fill = cell_type)) + geom_col() +
  scale_y_continuous(labels = scales::percent) + theme_bw() + labs(title = "Cell composition by patient", x = NULL, y = "Fraction of retained singlets", fill = "Cell type")
save_plot(p, "figures/composition/06_sample_composition", 11, 6)
p <- ggplot(primary$composition, aes(side, proportion, color = side)) +
  geom_point(position = position_jitter(width = .08, seed = 20260915), size = 2.5) +
  facet_wrap(~cell_type, scales = "free_y", ncol = 4) + scale_color_manual(values = side_colors) +
  scale_y_continuous(labels = scales::percent) + theme_bw() + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none") +
  labs(title = "Each point is one patient | n = 3 per side", x = NULL, y = "Cell fraction")
save_plot(p, "figures/composition/06_patient_comparison", 13, 8)
print(primary$tests); stage_end()
