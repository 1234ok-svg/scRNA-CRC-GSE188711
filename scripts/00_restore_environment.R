# 输入：项目 renv.lock；输出：按照锁定版本恢复的项目独立包库。
# 需要网络；部分历史包可能需要 Rtools/系统编译依赖。未在第二台机器完成恢复验证。
stopifnot(file.exists("renv.lock"), file.exists("scRNA_CRC_project.Rproj"))
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
set.seed(20260915)
dir.create("renv/library-local",recursive=TRUE,showWarnings=FALSE)
.libPaths(c("renv/library-local",.libPaths()))
if(!requireNamespace("renv",quietly=TRUE))
  install.packages("renv",repos="https://cloud.r-project.org",lib="renv/library-local")
renv::restore(lockfile="renv.lock",library="renv/library-local",prompt=FALSE)
dir.create("logs",showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),"logs/restored_sessionInfo.txt")
