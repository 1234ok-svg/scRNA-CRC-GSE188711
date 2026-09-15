# 输入：各阶段实际结果；输出：结果图、统计摘要及 RDS。
# 所有摘要数字从结果文件计算；不把探索性候选写成已验证生物学机制。
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source("scripts/00_common.R", encoding="UTF-8")
suppressPackageStartupMessages({library(Seurat); library(edgeR); library(ggplot2)})
stage_start("10_results")
ensure_dirs("results/summary","figures/overview")
object <- readRDS("data/processed/05_seurat_annotated.rds")
sm <- read.csv("data/metadata/sample_metadata.csv")
loading <- read.csv("results/data_loading/sample_loading_summary.csv")
dbl <- read.csv("results/qc/doublet_summary.csv")
comp <- read.csv("results/composition/composition_tests.csv")
res <- readRDS("data/processed/07_differential_results.rds")
primary <- res[res$scenario=="primary_mt25", ]
de <- primary[primary$significant, ]
loo <- read.csv("results/differential/leave_one_patient_out_direction.csv")
mt <- res[res$scenario=="sensitivity_mt20", ]
key <- function(d) paste(d$cell_type,d$gene_id,sep="__")
de$loo_direction_consistency <- loo$loo_direction_consistency[match(key(de),key(loo))]
de$mt20_logFC <- mt$logFC[match(key(de),key(mt))]
de$mt20_FDR <- mt$FDR[match(key(de),key(mt))]
de$mt20_pass <- mt$significant[match(key(de),key(mt))]
de <- de[order(de$cell_type,de$FDR), ]
write.csv(de,"results/summary/candidate_evidence.csv",row.names=FALSE)
types <- sort(unique(object$cell_type))
type_counts <- as.data.frame(table(object$cell_type)); names(type_counts)<-c("cell_type","n_cells")
write.csv(type_counts,"results/summary/cell_type_counts.csv",row.names=FALSE)
sample_counts <- data.frame(sample_id=sm$sample_id,side=sm$side,
 n_input=loading$n_cells[match(sm$sample_id,loading$sample_id)],
 n_basic_qc=dbl$n_before[match(sm$sample_id,dbl$sample_id)],
 n_doublet=dbl$n_doublet[match(sm$sample_id,dbl$sample_id)],
 n_final=dbl$n_singlet[match(sm$sample_id,dbl$sample_id)])
write.csv(sample_counts,"results/summary/sample_cell_flow.csv",row.names=FALSE)
# 固定并区分 13 种颜色，用同一颜色对应 UMAP 和组成图。
pal <- c(B_cells="#E69F00",Dendritic="#B79F00",Endothelial="#009E73",Epithelial="#8C6D31",
 Fibroblasts="#66A61E",Macrophages="#0072B2",Mast_cells="#A6761D",Monocytes="#56B4E9",
 NK_cells="#D73027",Neutrophils="#7570B3",Pericytes="#E78AC3",Plasma_cells="#CC79A7",T_cells="#4055A8")
p <- DimPlot(object,reduction="umap",group.by="cell_type",cols=pal,label=TRUE,repel=TRUE,
 raster=FALSE,pt.size=.1,shuffle=TRUE,seed=20260915) +
 labs(title="GSE188711 | 25,765 retained cells",subtitle="Six patients; marker-supported major lineages | NK-like annotation remains tentative")
save_plot(p,"figures/overview/10_cell_atlas",12,8)
cs <- read.csv("results/composition/sample_composition.csv")
p <- ggplot(cs,aes(sample_id,proportion,fill=cell_type)) + geom_col(width=.7) +
 facet_grid(~side,scales="free_x",space="free_x") + scale_fill_manual(values=pal) +
 scale_y_continuous(labels=scales::percent) + theme_bw() + labs(title="Cell composition by patient",x=NULL,y="Fraction of retained cells",fill="Cell type")
save_plot(p,"figures/overview/10_patient_composition",11,6)
# 只展示真实候选基因；每点为一个患者，不使用细胞级小提琴制造样本量错觉。
pbs <- readRDS("data/processed/07_pseudobulk_models.rds"); ex <- list()
for(i in seq_len(nrow(de))) {
 ct <- de$cell_type[i]; pb <- pbs[[paste0("primary_mt25__",ct)]]
 y <- DGEList(pb$counts[pb$tested_genes, ],lib.size=pb$norm_factors$lib.size,norm.factors=pb$norm_factors$norm.factors)
 v <- cpm(y,log=TRUE,prior.count=2)[de$gene_id[i], ]
 ex[[i]] <- data.frame(sample_id=names(v),side=sm$side[match(names(v),sm$sample_id)],
   logCPM=as.numeric(v),gene_symbol=de$gene_symbol[i],cell_type=ct,panel=paste(ct,de$gene_symbol[i],sep=" | "))
}
if(length(ex)) {
 ex <- do.call(rbind,ex); write.csv(ex,"results/summary/candidate_patient_expression.csv",row.names=FALSE)
 p <- ggplot(ex,aes(side,logCPM,color=side,label=sample_id)) + geom_point(size=2.5) +
   ggrepel::geom_text_repel(nudge_x=.14,direction="y",seed=20260915,max.overlaps=Inf,
                           box.padding=.35,min.segment.length=0,size=3,show.legend=FALSE) +
   facet_wrap(~panel,scales="free_y",ncol=3) +
   scale_color_manual(values=side_colors) + theme_bw() + theme(legend.position="none",axis.text.x=element_text(angle=20,hjust=1)) +
   labs(title="Exploratory expression candidates | each point is a patient",subtitle="Within-type FDR < 0.05; none pass global FDR < 0.05",x=NULL,y="TMM log2 CPM (prior count = 2)")
 save_plot(p,"figures/overview/10_candidate_expression",12,7)
}
gsea <- read.csv("results/enrichment/gsea_go_bp_all.csv")
gsea_sig <- gsea[is.finite(gsea$FDR_global)&gsea$FDR_global<.05, ]
gsea_top <- head(gsea_sig[order(gsea_sig$FDR_global),c("cell_type","Description","NES","FDR_global")],12)
status <- read.csv("results/differential/analysis_status.csv")
summary <- data.frame(metric=c("patients","raw_cells","basic_qc_cells","predicted_doublets","final_singlets","major_cell_types","DE_types_eligible","candidate_gene_type_pairs","global_DE_FDR05","global_GSEA_FDR05"),
 value=c(nrow(sm),sum(loading$n_cells),sum(dbl$n_before),sum(dbl$n_doublet),ncol(object),length(types),sum(status$scenario=="primary_mt25"&status$status=="OK"),nrow(de),sum(primary$FDR_global<.05),nrow(gsea_sig)))
write.csv(summary,"results/summary/project_summary.csv",row.names=FALSE)
saveRDS(list(summary=summary,samples=sample_counts,types=type_counts,candidates=de,composition=comp,gsea_top=gsea_top),"data/processed/10_project_summary.rds")
stage_end()
