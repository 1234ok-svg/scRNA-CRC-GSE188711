# 输入：R 4.3.x、CRAN/Bioconductor 3.18 软件仓库；输出：项目独立包库和包版本表。
# 从项目根目录运行。本脚本显式联网安装；其他分析脚本不自动安装软件。
# 使用与 R 4.3 对应的历史二进制仓库，避免旧全局包库混用问题。
stopifnot(file.exists("scRNA_CRC_project.Rproj"))
if (.Platform$OS.type != "windows" || getRversion() < "4.3.0" || getRversion() >= "4.4.0")
  stop("此安装脚本针对 Windows R 4.3.x。其他系统请使用 renv.lock 恢复，并重新验证结果。")
Sys.setlocale("LC_CTYPE", "English_United States.utf8")
set.seed(20260915)
options(timeout = 180, download.file.method = "libcurl")
dir.create("renv/library-local", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)
dir.create("data/software", recursive = TRUE, showWarnings = FALSE)
.libPaths(c("renv/library-local", .libPaths()))
repos <- c(CRAN = "https://cloud.r-project.org",
           BioCsoft = "https://bioconductor.statistik.tu-dortmund.de/packages/3.18/bioc",
           BioCann = "https://bioconductor.statistik.tu-dortmund.de/packages/3.18/data/annotation",
           BioCexp = "https://bioconductor.statistik.tu-dortmund.de/packages/3.18/data/experiment")
targets <- c("Seurat", "harmony", "scDblFinder", "SingleCellExperiment", "edgeR",
             "clusterProfiler", "org.Hs.eg.db", "ggplot2", "patchwork", "dplyr", "renv", "BiocManager", "BiocVersion")
ap <- if (file.exists("data/software/available_binary_packages.rds")) readRDS("data/software/available_binary_packages.rds") else
  available.packages(contriburl = contrib.url(repos, type = "win.binary"))
ap[, "Repository"] <- sub("https://bioconductor.org", "https://bioconductor.statistik.tu-dortmund.de", ap[, "Repository"], fixed = TRUE)
saveRDS(ap, "data/software/available_binary_packages.rds")
# 干净 R 环境可能尚无 curl；先用 R 自带下载器引导安装，避免安装器自身缺依赖。
if (!requireNamespace("curl", quietly=TRUE)) {
  curl_archive <- file.path("data/software",paste0("curl_",ap["curl","Version"],".zip"))
  download.file(paste0(ap["curl","Repository"],"/",basename(curl_archive)),curl_archive,mode="wb")
  install.packages(curl_archive,repos=NULL,type="win.binary",lib="renv/library-local")
}
# 注释数据库通常仅提供 source，后面单独安装；纯 R 数据包不需要编译器。
binary_targets <- intersect(targets, rownames(ap))
deps <- unique(c(binary_targets, unlist(tools::package_dependencies(binary_targets, db = ap,
                              which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE))))
deps <- intersect(deps, rownames(ap))
base_pkgs <- rownames(installed.packages(priority = "base"))
deps <- setdiff(deps, base_pkgs)
local <- installed.packages(lib.loc = "renv/library-local")
todo <- setdiff(deps, rownames(local))
cat("Installing", length(todo), "binary packages into project library\n")
if (length(todo)) {
  filenames <- paste0(todo, "_", ap[todo, "Version"], ".zip")
  urls <- paste0(ap[todo, "Repository"], "/", filenames)
  paths <- file.path("data/software", filenames)
  valid_zip <- function(p) file.exists(p) && isTRUE(tryCatch(nrow(unzip(p, list = TRUE)) > 0,
                                                            error = function(e) FALSE, warning = function(w) FALSE))
  # 小批次并发下载；已有完整 zip 复用，网络中断后可重跑。
  for (attempt in 1:3) {
    pending <- which(!vapply(paths, valid_zip, logical(1)))
    if (!length(pending)) break
    for (idx in split(pending, ceiling(seq_along(pending)/12))) {
      ans <- curl::multi_download(urls[idx], destfiles = paths[idx], resume = FALSE,
                                   timeout = 120, progress = FALSE)
      cat("Download batch:", sum(ans$success), "/", length(idx), "\n")
    }
  }
  good <- vapply(paths, valid_zip, logical(1))
  if (!all(good)) warning("Unavailable packages: ", paste(todo[!good], collapse = ", "))
  install.packages(paths[good], repos = NULL, type = "win.binary", lib = "renv/library-local")
}
source_versions <- c("GO.db" = "3.18.0", "org.Hs.eg.db" = "3.18.0", "GenomeInfoDbData" = "1.2.11", "HDO.db" = "0.99.1")
for (p in names(source_versions)) {
  if (!p %in% rownames(installed.packages(lib.loc = "renv/library-local"))) {
    archive <- paste0(p, "_", source_versions[[p]], ".tar.gz")
    path <- file.path("data/software", archive)
    # 可复用已完整下载的注释数据库，不在每次安装时重新传输大型文件。
    archive_ok <- file.exists(path) && isTRUE(tryCatch(length(untar(path, list = TRUE)) > 0,
                                             error = function(e) FALSE, warning = function(w) FALSE))
    if (!archive_ok) download.file(paste0(repos[["BioCann"]], "/src/contrib/", archive), path, mode = "wb")
    install.packages(path, repos = NULL, type = "source", lib = "renv/library-local")
  }
}
status <- lapply(targets, function(p) {
  err <- tryCatch({loadNamespace(p); "OK"}, error = function(e) conditionMessage(e))
  data.frame(package = p, status = err)
})
write.csv(do.call(rbind, status), "logs/environment_status.csv", row.names = FALSE)
write.csv(installed.packages(lib.loc = "renv/library-local")[, c("Package", "Version", "Built")],
          "logs/package_versions.csv", row.names = FALSE)
print(do.call(rbind, status))
if (any(vapply(status, function(x) x$status != "OK", logical(1)))) stop("Some dependencies require repair; see environment_status.csv")
