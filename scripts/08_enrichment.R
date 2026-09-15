# 输入：pseudobulk DE 结果；输出：按类型/方向的 GO BP ORA、RDS 和富集图。
# 背景为该类型实际通过表达过滤并参与检验的基因，不使用全部人类基因作背景。
# org.Hs.eg.db 记录于 sessionInfo；转换 ID 后去重，保存转换覆盖率。
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source("scripts/00_common.R", encoding = "UTF-8")
suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)})
stage_start("08_enrichment")
ensure_dirs("results/enrichment", "figures/enrichment")
res <- readRDS("data/processed/07_differential_results.rds")
res <- res[res$scenario == "primary_mt25", ]
idmap <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(res$gene_id), keytype = "ENSEMBL", columns = "ENTREZID")
idmap <- idmap[!is.na(idmap$ENTREZID), ]
write.csv(idmap, "results/enrichment/ensembl_entrez_mapping.csv", row.names = FALSE)
entrez <- function(ids) unique(idmap$ENTREZID[idmap$ENSEMBL %in% ids])
all_enrich <- list(); objects <- list(); statuses <- list()
for (ct in unique(res$cell_type)) {
  d <- res[res$cell_type == ct, ]; universe <- entrez(d$gene_id)
  for (direction in c("higher_right", "higher_left")) {
    key <- paste(ct, direction, sep = "__")
    sig <- d$gene_id[d$significant & d$direction == direction]; selected <- entrez(sig)
    status <- data.frame(cell_type = ct, direction = direction, n_selected_ensembl = length(sig),
                          n_selected_entrez = length(selected), n_background = length(universe),
                          n_GO_FDR05 = 0, status = "OK")
    if (length(selected) < 5) {
      status$status <- "SKIPPED_fewer_than_5_mapped_DE_genes"; statuses[[key]] <- status; next
    }
    ego <- enrichGO(gene = selected, universe = universe, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                    minGSSize = 10, maxGSSize = 500, readable = TRUE)
    objects[[key]] <- ego
    tab <- as.data.frame(ego)
    if (nrow(tab)) {
      tab$cell_type <- ct; tab$direction <- direction; all_enrich[[key]] <- tab
      status$n_GO_FDR05 <- sum(tab$p.adjust < .05)
    } else status$status <- "NO_GO_TERMS"
    statuses[[key]] <- status
  }
}
write.csv(do.call(rbind, statuses), "results/enrichment/enrichment_status.csv", row.names = FALSE)
saveRDS(objects, "data/processed/08_go_enrichment.rds")
if (length(all_enrich)) {
  tab <- do.call(rbind, all_enrich); tab$FDR_global <- p.adjust(tab$pvalue, "BH")
  write.csv(tab, "results/enrichment/go_bp_all.csv", row.names = FALSE)
  write.csv(tab[tab$p.adjust < .05, ], "results/enrichment/go_bp_significant.csv", row.names = FALSE)
  # 每类型方向最多 5 项，保留原始完整表；不将相似 GO 项当作独立生物学证据。
  sig <- tab[tab$p.adjust < .05, ]
  if (nrow(sig)) {
    sig$key <- paste(sig$cell_type, sig$direction, sep = " | ")
    top <- do.call(rbind, lapply(split(sig, sig$key), function(x) head(x[order(x$p.adjust), ], 5)))
    top$label <- paste(top$Description, top$key, sep = "\n")
    p <- ggplot(top, aes(-log10(p.adjust), reorder(label, -log10(p.adjust)), color = direction, size = Count)) +
      geom_point() + theme_bw(base_size = 10) + labs(x = "-log10 within-analysis GO FDR", y = NULL,
         title = "GO Biological Process | patient-level DE gene sets") +
      scale_color_manual(values = c("higher_right" = "#CC7546", "higher_left" = "#287D8E"))
    save_plot(p, "figures/enrichment/08_go_bp", 13, max(5, min(20, nrow(top)*.32)))
  }
} else {
  empty <- data.frame(ID=character(), Description=character(), pvalue=numeric(), p.adjust=numeric(), cell_type=character(), direction=character(), FDR_global=numeric())
  write.csv(empty, "results/enrichment/go_bp_all.csv", row.names = FALSE)
  write.csv(empty, "results/enrichment/go_bp_significant.csv", row.names = FALSE)
  # 空结果也给出明确可展示的状态图，不伪造通路条目。
  p <- ggplot() + annotate("text", x=0, y=0,
    label="GO ORA: no eligible input gene sets\nEach cell type/direction has fewer than 5 mapped DE genes\nSee enrichment_status.csv", size=5) +
    xlim(-1,1) + ylim(-1,1) + theme_void()
  save_plot(p, "figures/enrichment/08_go_bp", 11, 4)
}
stage_end()
