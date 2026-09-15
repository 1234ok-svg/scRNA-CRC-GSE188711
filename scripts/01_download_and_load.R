# 项目：基于单细胞转录组解析左右侧结直肠癌肿瘤微环境差异
# 阶段 01：公开 count matrix 下载、读取和 Seurat 对象构建。
# 输入：GSE188711 的 18 个 MTX/TSV gzip 文件及下面的显式样本映射。
# 输出：data/metadata 中的来源与校验记录；data/processed 中的 RDS；
#       results/data_loading 中的读取统计；logs 中的日志和运行环境。
# 运行：从含 .Rproj 的项目根目录执行 Rscript --vanilla scripts/01_download_and_load.R
# 本阶段不归一化、不筛选细胞、不聚类、不注释；合并不代表批次校正。

main <- function() {
  # 1. 输入：项目工作目录。输出：所需目录、固定随机数状态及运行日志。
  # 不依赖个人计算机的绝对路径；错误工作目录立即报错，避免输出到未知位置。
  if (!file.exists("scRNA_CRC_project.Rproj")) {
    stop("请先将工作目录设为含 scRNA_CRC_project.Rproj 的项目根目录。")
  }
  set.seed(20260915)
  options(timeout = max(600, getOption("timeout")), stringsAsFactors = FALSE)
  if (dir.exists("renv/library-local")) .libPaths(c("renv/library-local", .libPaths()))
  dirs <- c("data/raw/GSE188711", "data/metadata", "data/processed",
            paste0("results/", c("data_loading", "qc", "clustering", "markers",
                                  "annotation", "composition", "differential", "enrichment")),
            "figures", "docs", "logs")
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  log_path <- "logs/01_run.log"
  cat("", file = log_path)
  log_msg <- function(...) {
    line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " "))
    cat(line, "\n")
    cat(line, "\n", file = log_path, append = TRUE)
  }
  on.exit(writeLines(capture.output(sessionInfo()), "logs/01_sessionInfo.txt"), add = TRUE)
  # 本阶段只需要稀疏矩阵和 Seurat 核心对象；优先用 Seurat::ReadMtx。
  # 若完整 Seurat 因与本阶段无关的依赖缺失而不可加载，使用 Matrix::readMM
  # 读取同一 MTX，并通过 SeuratObject 构建真正的 Seurat 对象，记录实际路径。
  required <- c("Matrix", "SeuratObject")
  for (pkg in required) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("缺少必要 R 包：", pkg)
  }
  use_seurat <- requireNamespace("Seurat", quietly = TRUE)
  log_msg("读取方式：", if (use_seurat) "Seurat::ReadMtx" else "Matrix::readMM + SeuratObject")

  # 2. 输入：GEO 样本标题、组织部位、补充文件名。输出：显式样本信息表。
  # patient_id 是本项目匿名标识；一位患者一个肿瘤样本，不是配对设计。
  # R2 的文件标识为 R_CRC3，R3 为 R_CRC4，不能按文件数字推断分组。
  samples <- data.frame(
    gsm = paste0("GSM", 5688706:5688711),
    sample_id = c("L1", "L2", "L3", "R1", "R2", "R3"),
    patient_id = paste0("P", 1:6),
    geo_title = c(paste("Left-sided CRC", 1:3), paste("Right-sided CRC", 1:3)),
    side = rep(c("left-sided", "right-sided"), each = 3),
    tumor_location = rep(c("Sigmoid", "Ascending"), each = 3),
    age = c(57L, 70L, 65L, 69L, 80L, 59L),
    sex = c("Female", "Male", "Male", "Female", "Female", "Male"),
    file_tag = c("WGC", "JCA", "LS-CRC3", "RS-CRC1", "R_CRC3", "R_CRC4")
  )
  samples$source_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", samples$gsm)
  write.csv(samples, "data/metadata/sample_metadata.csv", row.names = FALSE)
  saveRDS(samples, "data/processed/01_sample_metadata.rds")

  # 3. 输入：HTTPS 文件链接。输出：原始 gzip 文件及 MD5 下载清单。
  # 先完整解压扫描以检测截断文件，再替换目标文件；失败自动重试三次。
  # 只有 gzip 检查通过且匹配上次本地 MD5 才跳过已有文件。
  # MD5 是本地复现记录，不宣称来自 GEO 官方，也不能证明来源真实性。
  check_gzip <- function(path) {
    if (!file.exists(path) || file.info(path)$size < 20) return(FALSE)
    tryCatch({
      signature <- readBin(path, what = "raw", n = 2L)
      if (!identical(signature, as.raw(c(31, 139)))) return(FALSE)
      con <- gzfile(path, "rb")
      on.exit(close(con), add = TRUE)
      withCallingHandlers({
        repeat { if (!length(readBin(con, what = "raw", n = 1048576L))) break }
      }, warning = function(w) stop(conditionMessage(w)))
      TRUE
    }, error = function(e) FALSE)
  }
  manifest_path <- "data/metadata/download_manifest.csv"
  previous <- if (file.exists(manifest_path)) read.csv(manifest_path) else data.frame()
  manifest <- list()
  download_one <- function(url, path) {
    prior <- if (nrow(previous)) match(path, previous$local_path) else NA_integer_
    unchanged <- !is.na(prior) && file.exists(path) &&
      identical(unname(tools::md5sum(path)), as.character(previous$md5[prior]))
    if (unchanged && check_gzip(path)) {
      log_msg("复用已校验文件：", basename(path))
      return("cached")
    }
    tmp <- paste0(path, ".part")
    for (attempt in 1:3) {
      log_msg("下载", basename(path), "尝试", attempt)
      ok <- tryCatch({
        status <- download.file(url, tmp, mode = "wb", method = "libcurl", quiet = TRUE)
        identical(status, 0L) && check_gzip(tmp)
      }, error = function(e) { log_msg("下载错误：", conditionMessage(e)); FALSE })
      if (ok) {
        if (!file.copy(tmp, path, overwrite = TRUE)) stop("无法保存文件：", path)
        unlink(tmp)
        return("downloaded")
      }
      if (file.exists(tmp)) unlink(tmp)
    }
    stop("下载失败，重新运行可复用已完成文件：", url)
  }
  types <- c("matrix", "features", "barcodes")
  extensions <- c("mtx.gz", "tsv.gz", "tsv.gz")
  paths_by_sample <- list()
  for (i in seq_len(nrow(samples))) {
    sample <- samples[i, ]
    sample_dir <- file.path("data/raw/GSE188711", sample$gsm)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    paths <- setNames(character(3), types)
    for (j in seq_along(types)) {
      filename <- paste0(sample$gsm, "_", types[j], "_", sample$file_tag, ".", extensions[j])
      path <- file.path(sample_dir, filename)
      url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5688nnn/",
                    sample$gsm, "/suppl/", filename)
      status <- download_one(url, path)
      manifest[[length(manifest) + 1L]] <- data.frame(
        gsm = sample$gsm, file_type = types[j], url = url, local_path = path,
        bytes = file.info(path)$size, md5 = unname(tools::md5sum(path)),
        validated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), status = status)
      # 每完成一个文件立即保存，网络中断后也可复用已完成部分。
      write.csv(do.call(rbind, manifest), manifest_path, row.names = FALSE)
      paths[j] <- path
    }
    paths_by_sample[[sample$sample_id]] <- paths
  }

  # 4. 输入：每样本 MTX + features + barcodes。输出：经过一致性检查的 Seurat 列表。
  # 以第一列稳定基因 ID 为行名，防止重复 gene symbol 造成跨样本错误合并。
  # 原始 symbol 与 feature type 另存映射；后续 marker/富集需通过该表转换。
  objects <- list()
  feature_maps <- list()
  summaries <- list()
  for (i in seq_len(nrow(samples))) {
    sample <- samples[i, ]
    paths <- paths_by_sample[[sample$sample_id]]
    log_msg("读取样本：", sample$sample_id, sample$gsm)
    features <- read.delim(paths[["features"]], header = FALSE, quote = "", comment.char = "",
                           colClasses = "character", check.names = FALSE)
    barcodes <- read.delim(paths[["barcodes"]], header = FALSE, quote = "", comment.char = "",
                           colClasses = "character")
    if (ncol(features) < 2 || ncol(barcodes) != 1) stop("未知 TSV 列结构：", sample$gsm)
    ids <- features[[1]]
    cells <- barcodes[[1]]
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("基因 ID 缺失或重复：", sample$gsm)
    if (anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) stop("barcode 缺失或重复：", sample$gsm)
    if (ncol(features) >= 3 && any(features[[3]] != "Gene Expression")) {
      stop("发现非 Gene Expression 特征，需要明确选择模态：", sample$gsm)
    }
    if (use_seurat) {
      counts <- Seurat::ReadMtx(mtx = paths[["matrix"]], cells = paths[["barcodes"]],
                                features = paths[["features"]], feature.column = 1,
                                unique.features = FALSE, strip.suffix = FALSE)
    } else {
      counts <- local({
        con <- gzfile(paths[["matrix"]], "rt")
        on.exit(close(con))
        Matrix::readMM(con)
      })
      if (!identical(dim(counts), c(length(ids), length(cells)))) stop("矩阵与 TSV 维度不符")
      dimnames(counts) <- list(ids, cells)
    }
    counts <- methods::as(counts, "CsparseMatrix")
    stopifnot(identical(dim(counts), c(length(ids), length(cells))),
              identical(rownames(counts), ids), identical(colnames(counts), cells))
    if (any(!is.finite(counts@x)) || any(counts@x < 0) || any(counts@x != floor(counts@x))) {
      stop("矩阵不是有限非负整数 counts：", sample$gsm)
    }
    # min.cells/min.features 显式置零：本阶段不添加新的质控阈值。
    metadata <- sample[rep(1L, length(cells)), c("gsm", "sample_id", "patient_id", "side",
                                               "tumor_location", "age", "sex")]
    metadata$barcode_original <- cells
    cell_names <- paste(sample$sample_id, cells, sep = "_")
    rownames(metadata) <- cell_names
    colnames(counts) <- cell_names
    object <- SeuratObject::CreateSeuratObject(counts = counts, project = sample$sample_id,
                                               min.cells = 0, min.features = 0,
                                               meta.data = metadata)
    stopifnot(ncol(object) == length(cells), nrow(object) == length(ids))
    feature_maps[[sample$sample_id]] <- data.frame(
      gsm = sample$gsm, gene_id = ids, seurat_feature = rownames(object),
      gene_symbol = features[[2]],
      feature_type = if (ncol(features) >= 3) features[[3]] else "Gene Expression")
    summaries[[sample$sample_id]] <- data.frame(
      gsm = sample$gsm, sample_id = sample$sample_id, side = sample$side,
      n_genes = nrow(object), n_cells = ncol(object), nonzero_entries = Matrix::nnzero(counts),
      total_counts = sum(counts), median_counts = median(object$nCount_RNA),
      median_features = median(object$nFeature_RNA),
      zero_count_cells = sum(object$nCount_RNA == 0),
      duplicate_symbols = sum(duplicated(features[[2]])))
    objects[[sample$sample_id]] <- object
    log_msg("已构建：", nrow(object), "基因 ×", ncol(object), "细胞")
  }
     feature_map <- do.call(rbind, feature_maps)
  write.csv(feature_map, "data/metadata/gene_feature_mapping.csv", row.names = FALSE)
  saveRDS(feature_map, "data/processed/01_gene_feature_mapping.rds")
  saveRDS(objects, "data/processed/01_seurat_list.rds")

  # 5. 输入：六个样本对象。输出：保留样本标签和 counts 的合并对象。
  # merge.data=FALSE 不合并任何标准化数据；没有运行整合或批次校正。
  merged <- merge(objects[[1]], y = objects[-1], merge.data = FALSE,
                  project = "GSE188711_CRC")
  stopifnot(ncol(merged) == sum(vapply(objects, ncol, integer(1))),
            !anyDuplicated(colnames(merged)),
            !anyNA(merged$side), length(unique(merged$gsm)) == 6L,
            nrow(merged) == length(unique(unlist(lapply(objects, rownames)))))
  summary_table <- do.call(rbind, summaries)
  for (i in seq_len(nrow(summary_table))) {
    selected <- merged$sample_id == summary_table$sample_id[i]
    stopifnot(sum(selected) == summary_table$n_cells[i],
              sum(merged$nCount_RNA[selected]) == summary_table$total_counts[i])
  }
  # 保存 provenance，明确对象起点和本次未做的分析，避免误用初始对象。
  merged@misc$stage01 <- list(dataset = "GSE188711", seed = 20260915L,
                             input = "GEO public gene-barcode counts",
                             additional_qc_filtering = FALSE,
                             normalized = FALSE, integrated = FALSE, annotated = FALSE,
                             feature_key = "gene_id", reader = if (use_seurat) "ReadMtx" else "readMM")
  saveRDS(merged, "data/processed/01_seurat_merged.rds")
  saveRDS(summary_table, "data/processed/01_loading_summary.rds")
  write.csv(summary_table, "results/data_loading/sample_loading_summary.csv", row.names = FALSE)
  # 回读磁盘 RDS，确认可读取且核心尺寸/元数据在序列化后保持一致。
  restored <- readRDS("data/processed/01_seurat_merged.rds")
  stopifnot(inherits(restored, "Seurat"), identical(dim(restored), dim(merged)),
            identical(restored[[]], merged[[]]))
  # 不只比较元数据：直接核对磁盘对象的每个细胞 counts 总和及检测基因数。
  # 兼容 SeuratObject v4 的单 counts slot 和 v5 的多 counts layers。
  if (inherits(restored[["RNA"]], "Assay5")) {
    count_layers <- SeuratObject::Layers(restored[["RNA"]], search = "^counts")
    actual_umi <- actual_features <- numeric()
    for (layer in count_layers) {
      mat <- SeuratObject::LayerData(restored, assay = "RNA", layer = layer)
      actual_umi <- c(actual_umi, Matrix::colSums(mat))
      actual_features <- c(actual_features, Matrix::colSums(mat > 0))
    }
  } else {
    mat <- SeuratObject::GetAssayData(restored, assay = "RNA", slot = "counts")
    actual_umi <- Matrix::colSums(mat)
    actual_features <- Matrix::colSums(mat > 0)
  }
  stopifnot(!anyDuplicated(names(actual_umi)),
            setequal(names(actual_umi), colnames(restored)),
            all(actual_umi[colnames(restored)] == restored$nCount_RNA),
            all(actual_features[colnames(restored)] == restored$nFeature_RNA))
  log_msg("成功完成；合并对象：", nrow(merged), "基因 ×", ncol(merged), "细胞")
  print(summary_table, row.names = FALSE)
  invisible(summary_table)
}

main()
