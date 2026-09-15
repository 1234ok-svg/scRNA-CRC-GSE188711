# 输入：所有通过 edgeR 表达过滤的基因及 QL F 统计量；输出：探索性 GO BP GSEA。
# 此补充分析在看到稀少 DE 基因后增加，不是预注册验证，不降低 DE 阈值来凑通路。
# 排名 = sign(logFC) * sqrt(F)。采用基因排序置换，不能替代患者重复或独立队列验证。
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source("scripts/00_common.R", encoding="UTF-8")
suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)})
stage_start("08b_ranked_enrichment")
BiocParallel::register(BiocParallel::SerialParam())
ensure_dirs("results/enrichment", "figures/enrichment", "results/enrichment/gsea_checkpoints")
res <- readRDS("data/processed/07_differential_results.rds")
res <- res[res$scenario == "primary_mt25", ]
idmap <- read.csv("results/enrichment/ensembl_entrez_mapping.csv", colClasses="character")
# 一个 Ensembl 对多个 Entrez 的模糊映射不进入排名；多个 Ensembl 对一个 Entrez 时
# 保留平均表达量最高的基因，避免按最显著结果选择代表基因。
idmap <- unique(idmap)
ambiguous <- names(which(table(idmap$ENSEMBL)>1))
idmap <- idmap[!idmap$ENSEMBL %in% ambiguous, ]
results <- list(); objects <- list(); ranks <- list(); status <- list()
for (ct in sort(unique(res$cell_type))) {
 d <- res[res$cell_type == ct, ]
 d$ENTREZID <- idmap$ENTREZID[match(d$gene_id,idmap$ENSEMBL)]
 d <- d[!is.na(d$ENTREZID), ]; d <- d[order(-d$logCPM,d$gene_id), ]
 d <- d[!duplicated(d$ENTREZID), ]
 d$rank_stat <- sign(d$logFC)*sqrt(pmax(d$F,0))
 d <- d[is.finite(d$rank_stat), ]; d <- d[order(-d$rank_stat,d$ENTREZID), ]
 ranks[[ct]] <- d[,c("gene_id","gene_symbol","ENTREZID","cell_type","rank_stat")]
 gl <- setNames(d$rank_stat,d$ENTREZID)
 set.seed(20260915)
 g <- gseGO(geneList=gl, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="BP",
   minGSSize=15, maxGSSize=500, pvalueCutoff=1, pAdjustMethod="BH", eps=1e-10,
   seed=TRUE, by="fgsea", verbose=FALSE, nPermSimple=1000)
 objects[[ct]] <- g; saveRDS(g,paste0("results/enrichment/gsea_checkpoints/",ct,".rds"))
 tab <- as.data.frame(g)
 if(nrow(tab)) {tab$cell_type <- ct; results[[ct]] <- tab}
 status[[ct]] <- data.frame(cell_type=ct,n_ranked=nrow(d),n_terms=nrow(tab),
    n_FDR05=if(nrow(tab))sum(tab$p.adjust<.05,na.rm=TRUE) else 0)
 log_msg(ct,"ranked genes",nrow(d),"GO terms",nrow(tab))
}
write.csv(do.call(rbind,ranks),"results/enrichment/gsea_gene_ranks.csv",row.names=FALSE)
write.csv(do.call(rbind,status),"results/enrichment/gsea_status.csv",row.names=FALSE)
saveRDS(objects,"data/processed/08b_gsea_objects.rds")
tab <- do.call(rbind,results)
if(!is.null(tab) && nrow(tab)) {
 tab$FDR_global <- p.adjust(tab$pvalue,"BH")
 write.csv(tab,"results/enrichment/gsea_go_bp_all.csv",row.names=FALSE)
 write.csv(tab[!is.na(tab$FDR_global) & tab$FDR_global<.05, ],"results/enrichment/gsea_go_bp_global_significant.csv",row.names=FALSE)
 # 各类型最多两项正向和两项负向，概览最多 20 项；报告全局 FDR，不隐去完整结果。
 sig <- tab[is.finite(tab$FDR_global) & tab$FDR_global<.05, ]
 if(nrow(sig)) {
   top <- do.call(rbind,lapply(split(sig,paste(sig$cell_type,sign(sig$NES))),function(x) head(x[order(x$FDR_global), ],2)))
   top <- head(top[order(top$FDR_global), ],20)
   top$label <- paste(top$Description,top$cell_type,sep="\n")
   p <- ggplot(top,aes(NES,reorder(label,NES),size=setSize,color=-log10(pmax(FDR_global,1e-300)))) +
     geom_point() + geom_vline(xintercept=0,linetype=2,color="grey60") + theme_bw(base_size=10) +
     scale_color_viridis_c() + labs(title="Exploratory ranked GO GSEA | top terms",subtitle="NES > 0: right-sided; NES < 0: left-sided | gene-set test, not independent validation",y=NULL,color="-log10 global FDR")
   save_plot(p,"figures/enrichment/08b_ranked_go_bp",14,max(5,nrow(top)*.43))
 } else {
   p <- ggplot() + annotate("text",x=0,y=0,label="Ranked GO GSEA: no terms with global FDR < 0.05",size=5) + theme_void()
   save_plot(p,"figures/enrichment/08b_ranked_go_bp",11,4)
 }
}
stage_end()
