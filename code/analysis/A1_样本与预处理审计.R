# A1 样本与预处理审计（§3.1）
# 目的：对全部 bulk 数据集产出逐样本清单与逐数据集预处理清单，核实
#       分组规则、样本重叠、平台、原始尺度、探针映射、重复对数转换。
# 只读：不改动任何原始数据；输出写入 08_修订_2026-09-17/01_审计/
suppressPackageStartupMessages({library(Biobase)})

PROJ <- "/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/01.材料/01.博士材料/10.毕业/01.小论文/04.立项/滑膜基质刚度生成与刚度感知 0912"
RAW  <- file.path(PROJ, "04_分析数据/raw")
AUD  <- file.path(PROJ, "08_修订_2026-09-17/01_审计")
dir.create(AUD, showWarnings = FALSE, recursive = TRUE)

# ---- 与 03_preprocess.R 完全一致的分组规则（照抄，不改） ----
grp_by_title <- function(ttl, oa_pat, ra_pat, n_pat) {
  g <- rep("Other", length(ttl))
  g[grepl(oa_pat, ttl, ignore.case = TRUE)] <- "OA"
  g[grepl(ra_pat, ttl, ignore.case = TRUE)] <- "RA"
  g[grepl(n_pat,  ttl, ignore.case = TRUE)] <- "Normal"
  g
}
RULES <- list(
  GSE206848 = c("^OAS", "^RAS", "^NS"),
  GSE1919   = c("^OA_", "^RA_", "^ND_"),
  GSE55235  = c("osteoarthr", "rheumatoid", "healthy"),
  GSE55457  = c("osteoarthr", "rheumatoid", "Normal,"),
  GSE82107  = c("^OA-", "rheumatoid", "^HC-")
)
ROLE <- c(GSE89408 = "发现", GSE206848 = "发现", GSE1919 = "发现",
          GSE55235 = "复制", GSE55457 = "复制", GSE82107 = "复制")

# ================= 1. 逐样本清单 =================
samp_rows <- list()
for (gse in names(RULES)) {
  e <- readRDS(file.path(RAW, paste0(gse, "_1.rds")))
  pd <- pData(e)
  ttl <- pd$title
  g <- grp_by_title(ttl, RULES[[gse]][1], RULES[[gse]][2], RULES[[gse]][3])
  getf <- function(x) if (x %in% colnames(pd)) as.character(pd[[x]]) else rep(NA_character_, nrow(pd))
  samp_rows[[gse]] <- data.frame(
    dataset = gse, role = ROLE[[gse]], platform = annotation(e),
    gsm = rownames(pd), title = ttl,
    source = getf("source_name_ch1"),
    characteristics = getf("characteristics_ch1"),
    characteristics2 = getf("characteristics_ch1.1"),
    characteristics3 = getf("characteristics_ch1.2"),
    assigned_group = g, included = g %in% c("OA", "Normal"),
    stringsAsFactors = FALSE)
}

# GSE89408 是 RNA-seq，走 count 矩阵而非 series matrix
cnt <- read.delim(gzfile(file.path(RAW, "GSE89408/GSE89408_GEO_count_matrix_rename.txt.gz")),
                  check.names = FALSE, stringsAsFactors = FALSE)
cn <- sub("^X", "", colnames(cnt)[-1])
g89408 <- ifelse(grepl("^OA_tissue_", cn), "OA", ifelse(grepl("^normal_tissue_", cn), "Normal", "Other"))
samp_rows[["GSE89408"]] <- data.frame(
  dataset = "GSE89408", role = "发现", platform = "RNA-seq (count matrix)",
  gsm = cn, title = cn, source = NA, characteristics = NA, characteristics2 = NA,
  characteristics3 = NA, assigned_group = g89408,
  included = g89408 %in% c("OA", "Normal"), stringsAsFactors = FALSE)

samples <- do.call(rbind, samp_rows)
write.csv(samples, file.path(AUD, "样本清单_逐样本.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ================= 2. 逐数据集预处理清单 =================
ds_rows <- list()
for (gse in names(RULES)) {
  e <- readRDS(file.path(RAW, paste0(gse, "_1.rds")))
  x <- exprs(e)
  s <- samples[samples$dataset == gse, ]
  ds_rows[[gse]] <- data.frame(
    dataset = gse, role = ROLE[[gse]], platform = annotation(e),
    n_probes_raw = nrow(x), n_samples_raw = ncol(x),
    n_OA = sum(s$assigned_group == "OA"), n_Normal = sum(s$assigned_group == "Normal"),
    n_RA_excluded = sum(s$assigned_group == "RA"), n_other_excluded = sum(s$assigned_group == "Other"),
    scale_min = round(min(x), 3), scale_median = round(median(x), 3),
    scale_max = round(max(x), 3), frac_negative = round(100 * mean(x < 0), 2),
    already_logged = (max(x) < 25 && median(x) > 3),
    transform_applied = "log2(x+1)", stringsAsFactors = FALSE)
}
ds_rows[["GSE89408"]] <- data.frame(
  dataset = "GSE89408", role = "发现", platform = "RNA-seq count matrix",
  n_probes_raw = nrow(cnt) - 0, n_samples_raw = ncol(cnt) - 1,
  n_OA = sum(g89408 == "OA"), n_Normal = sum(g89408 == "Normal"),
  n_RA_excluded = 0, n_other_excluded = sum(g89408 == "Other"),
  scale_min = NA, scale_median = NA, scale_max = NA, frac_negative = NA,
  already_logged = FALSE,
  transform_applied = "CPM→log2(x+1)（03_preprocess）；DESeq2 另有 raw count 通路（05_de_meta）",
  stringsAsFactors = FALSE)

ds <- do.call(rbind, ds_rows)
write.csv(ds, file.path(AUD, "预处理清单_逐数据集.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ================= 3. 控制台审计报告 =================
cat("\n========== 逐数据集 ==========\n")
print(ds, row.names = FALSE)

cat("\n========== 分组规则命中检查（列出每队列全部 title 及其分组）==========\n")
for (gse in names(RULES)) {
  s <- samples[samples$dataset == gse, ]
  cat(sprintf("\n--- %s (%s, %s) : %d 样本, 纳入 %d (OA %d / Normal %d), 排除 %d\n",
      gse, ROLE[[gse]], annotation(readRDS(file.path(RAW, paste0(gse, "_1.rds")))),
      nrow(s), sum(s$included), sum(s$assigned_group == "OA"), sum(s$assigned_group == "Normal"),
      sum(!s$included)))
  tb <- table(s$title, s$assigned_group)
  print(tb)
}

cat("\n========== 样本重叠检查（跨数据集 GSM 是否有交集）==========\n")
gsm <- lapply(split(samples$gsm, samples$dataset), function(z) z)
for (i in seq_along(gsm)) for (j in seq_along(gsm)) {
  if (i < j) {
    ov <- intersect(gsm[[i]], gsm[[j]])
    if (length(ov)) cat(sprintf("  ⚠️ %s ∩ %s = %d 个: %s\n", names(gsm)[i], names(gsm)[j], length(ov), paste(head(ov,5), collapse=",")))
  }
}
cat("  （未列出即无 GSM 级重叠）\n")

cat("\n========== 输出 ==========\n")
cat("  01_审计/样本清单_逐样本.csv    ", nrow(samples), "行\n")
cat("  01_审计/预处理清单_逐数据集.csv", nrow(ds), "行\n")
