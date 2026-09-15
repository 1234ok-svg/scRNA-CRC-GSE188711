# 输入：项目根目录、已配置的 R 环境和数据；输出：顺序重建各阶段产物。
# 不自动安装软件。完整运行前先执行 00_setup_environment.R。
# 默认从阶段 01 开始；--from=3 可复用已有基础 QC 对象，从双细胞检测开始。
# 如果重新聚类得到新簇，05_annotation 会拒绝不匹配的注释表，应先审阅更新证据。
stopifnot(file.exists("scRNA_CRC_project.Rproj"))
set.seed(20260915)
args <- commandArgs(trailingOnly = TRUE)
arg <- grep("^--from=", args, value = TRUE)
from <- if (length(arg)) as.integer(sub("--from=", "", arg[1])) else 1L
stopifnot(!is.na(from), from %in% 1:10)
steps <- c("01_download_and_load.R", "02_quality_control.R", "03_doublets.R",
           "04_clustering.R", "05_annotation.R", "06_composition.R", "07_differential.R",
           "08_enrichment.R", "09_validation.R", "10_results.R")
for (i in seq.int(from, length(steps))) {
  if (i == 5) {
    # 再聚类先生成细胞映射和独立 marker 证据，再应用版本化的人工审阅标签。
    for (helper in c("05_refine_lymphocytes.R", "05_refinement_audit.R")) {
      status <- system2(file.path(R.home("bin"), "Rscript"),
                       c("--vanilla", shQuote(file.path("scripts", helper))))
      if (status != 0) stop("注释细化失败：", helper)
    }
  }
  cat("\nRunning", steps[i], "\n")
  flags <- if (i %in% c(2, 5)) "--apply" else character()
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--vanilla", shQuote(file.path("scripts", steps[i])), flags))
  if (status != 0) stop("阶段失败：", steps[i], "，请查看 logs 后从该阶段恢复。")
  if (i == 8) {
    status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote("scripts/08b_ranked_enrichment.R")))
    if (status != 0) stop("排序富集失败，请查看 logs/08b_ranked_enrichment.log")
  }
}
cat("All analysis stages completed. See results/ and figures/.\n")
