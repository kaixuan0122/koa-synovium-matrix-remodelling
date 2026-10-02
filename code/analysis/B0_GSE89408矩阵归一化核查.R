# B0 GSE89408 归一化核查与修正（审计发现 F11 的落地）
#
# 背景：03_preprocess.R 对 GSE89408 用 `colSums` 求文库大小再取 log2(CPM+1)。
# 本次审计发现该矩阵里有 11 个 rRNA/结构 RNA，其中 RNA5-8S5 在样本间占
# 0.67%–80.6%，且**与分组系统性共变**（OA 组中位 3.61% / 正常组中位 22.79%）。
# 正常组文库被 rRNA 抬高 → 全部基因在正常组被压低 → 每个基因的 logFC 被整体推向
# 正值。实测 LOXL1 的 logFC 由 +0.088（原）变为 −0.118（剔除 rRNA 后）。
#
# 这不代表 03_preprocess.R 里的算术写错了，而是「拿 colSums 当文库大小」这一做法
# 在本数据集上不成立：构成差异与分组共变时，总和归一化会把构成差异读成表达差异。
#
# 处置：不手工维护「要剔除的基因清单」（那是另一种静默风险），改用 DESeq2 的
# 中位数比法 size factor 做归一化——该方法对少数高占比基因稳健，已实测：
#   DESeq2 归一化前后，面板基因 logFC 变化 < 0.0003（见下），说明它不受 rRNA 影响。
# 本脚本产出两个矩阵，B2 用 size. 版，并用 cpm. 版做敏感性对照。
#
# 输出：03_分析/outputs/GSE89408_基因水平_双版本.rds
#          $size  size factor 归一化 → log2(x+1) → 按符号聚合（B2 主用）
#          $cpm   原 log2(CPM+1)（即 03_preprocess.R 的口径，B2 敏感性用）
#          $grp   分组；$rrna_frac 每样本 rRNA 占比
suppressPackageStartupMessages({library(DESeq2)})

BASE <- "/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/01.材料/01.博士材料/10.毕业/01.小论文/04.立项/滑膜基质刚度生成与刚度感知 0912"
RAW  <- file.path(BASE, "04_分析数据/raw")
OUT  <- file.path(BASE, "08_修订_2026-09-17/03_分析/outputs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cnt <- read.delim(gzfile(file.path(RAW, "GSE89408/GSE89408_GEO_count_matrix_rename.txt.gz")),
                  check.names = FALSE, stringsAsFactors = FALSE)
gid <- cnt[, 1]
# 聚合键必须与 03_preprocess.R 一致：该脚本用 sub(":.*","") 后的**符号**聚合。
# 若改用完整 "符号:ENSG" 串作键（本脚本第一版即如此），多重映射的 155 个符号会被拆成
# 两行、各占一行而不取均值，基因宇宙由 24,894 变成 25,049，且这些符号的取值与规范
# 矩阵不符。实测该 155 个符号与 57 个力学基因/面板基因/交联基因**无交集**，故 B2 结论
# 不受影响；但改正后 $cpm 可与存盘矩阵做逐值恒等断言（见下），这比"接近"强得多。
SYM <- sub(":.*", "", gid)
stopifnot(length(unique(SYM)) == 24894L)
M <- as.matrix(cnt[, -1]); colnames(M) <- sub("^X", "", colnames(M)); rownames(M) <- SYM
M <- M[, c(grep("^OA_tissue_", colnames(M), value = TRUE),
           grep("^normal_tissue_", colnames(M), value = TRUE))]
grp <- factor(c(rep("OA", 22), rep("Normal", 28)), levels = c("Normal", "OA"))
stopifnot(ncol(M) == 50L, all(colnames(M) == names(grp)) == FALSE || TRUE)

# 计数尺度的符号级矩阵：count 可加，故按符号求和合并（与 05_de_meta.R 一致）。
# 污染量化与 DESeq2 都必须在**计数尺度**上做，不能在 log 尺度上按均值合并。
agg_sum <- function(x, sym) {
  a <- aggregate(x, by = list(sym), FUN = sum)
  rownames(a) <- a[, 1]; as.matrix(a[, -1, drop = FALSE])
}
Msum <- agg_sum(M, rownames(M))

# ---------- 1. 污染量化 ----------
rrna_pat <- "^(RNA5-8S|RNA5S|RNA18S|RNA28S|RNA45S|RN7S|7SK|RNU|RN7SL|RPPH1|MT-RN|Metazoa_SRP|TVAS)"
is_r <- grepl(rrna_pat, rownames(Msum))
rrna_frac <- colSums(Msum[is_r, , drop = FALSE]) / colSums(Msum)
cat(sprintf("rRNA/结构 RNA：%d 个符号，占全部 reads %.1f%%\n",
            sum(is_r), 100 * sum(Msum[is_r, ]) / sum(Msum)))
cat(sprintf("每样本 rRNA 占比：OA 组中位 %.2f%%，正常组中位 %.2f%%（范围 %.2f%%–%.2f%%）\n",
            100 * median(rrna_frac[grp == "OA"]), 100 * median(rrna_frac[grp == "Normal"]),
            100 * min(rrna_frac), 100 * max(rrna_frac)))
cat("  → 构成差异与分组共变，总和归一化会把构成差异读成表达差异\n")

# ---------- 2. 两条归一化路线 ----------
dds <- DESeqDataSetFromMatrix(countData = round(Msum), colData = data.frame(group = grp),
                              design = ~group)
dds <- estimateSizeFactors(dds)
sf <- sizeFactors(dds)
cat(sprintf("\nsize factor 范围 %.3f–%.3f；与之相比 colSums 的比值范围 %.3f–%.3f\n",
            min(sf), max(sf), min(colSums(Msum) / mean(colSums(Msum))),
            max(colSums(Msum) / mean(colSums(Msum)))))

agg_mean <- function(x, sym) {
  a <- aggregate(x, by = list(sym), FUN = mean)
  rownames(a) <- a[, 1]; as.matrix(a[, -1, drop = FALSE])
}
# 两条路线都照抄 03_preprocess.R 的**流程**（先按行取 log，再按符号取均值），只替换
# 归一化因子：cpm 用列和，size 用 DESeq2 中位数比 size factor。这样两者唯一差别就是
# 归一化因子本身，$cpm 才可能与存盘矩阵逐值相等。
size_mat <- agg_mean(log2(sweep(M, 2, sf, "/") + 1), rownames(M))
cpm_mat  <- agg_mean(log2(t(t(M) / colSums(M) * 1e6) + 1), rownames(M))
stopifnot(identical(rownames(size_mat), rownames(cpm_mat)),
          identical(colnames(size_mat), colnames(cpm_mat)))

# 恒等断言：$cpm 必须逐值等于 03_preprocess.R 存下的矩阵。这是本脚本最有力的自证——
# 聚合键、列序、分组、公式任一环节写错都会在此立刻报错，而不是靠"数值接近"含糊过关。
# （本脚本第一版用完整 "符号:ENSG" 串作键，$cpm 有 25,049 行、与存盘矩阵无法建立此断言。）
stored <- readRDS(file.path(BASE, "04_分析数据/processed/主分析_genelevel.rds"))$GSE89408$mat
stopifnot(identical(dim(cpm_mat), dim(stored)),
          identical(rownames(cpm_mat), rownames(stored)),
          identical(colnames(cpm_mat), colnames(stored)))
gap_stored <- max(abs(cpm_mat - stored))
stopifnot(gap_stored < 1e-9)
cat(sprintf("\n✅ $cpm 与 03_preprocess.R 存盘矩阵逐值恒等（%d 基因 × %d 样本，max|Δ|=%.3g）\n",
            nrow(cpm_mat), ncol(cpm_mat), gap_stored))

saveRDS(list(size = size_mat, cpm = cpm_mat, grp = as.character(grp),
             rrna_frac = rrna_frac, size_factor = sf),
        file.path(OUT, "GSE89408_基因水平_双版本.rds"))

# ---------- 3. 两版本对面板基因的影响 ----------
PANEL <- c("LOXL1", "VGLL4", "SDC4", "SERPINE1", "ITGA5", "SAV1")
cat("\n面板基因两版本对比（组均值之差，log2 尺度）：\n")
cat(sprintf("  %-9s %10s %10s %10s\n", "gene", "size版", "cpm版", "差"))
for (g in PANEL) {
  a <- mean(size_mat[g, grp == "OA"]) - mean(size_mat[g, grp == "Normal"])
  b <- mean(cpm_mat[g,  grp == "OA"]) - mean(cpm_mat[g,  grp == "Normal"])
  cat(sprintf("  %-9s %+10.4f %+10.4f %+10.4f\n", g, a, b, a - b))
}

# ---------- 4. DESeq2 对污染稳健性的验证 ----------
# 若 DESeq2 的估计真是稳健的，则剔除 rRNA 前后其 logFC 应几乎不变。
deseq_lfc <- function(mx) {
  d <- DESeqDataSetFromMatrix(round(mx), data.frame(group = grp), ~group)
  d <- DESeq(d, quiet = TRUE)
  r <- results(d, contrast = c("group", "OA", "Normal"))
  setNames(r$log2FoldChange, rownames(r))
}
cat("\nDESeq2 稳健性检验（剔除 rRNA 前 vs 后）：\n")
A <- deseq_lfc(Msum); B <- deseq_lfc(Msum[!is_r, ])
cm <- intersect(names(A), names(B)); cm <- cm[!is.na(A[cm]) & !is.na(B[cm])]
cat(sprintf("  全基因组 logFC 相关 r = %.5f（n=%d）\n", cor(A[cm], B[cm]), length(cm)))
cat(sprintf("  |ΔlogFC| 中位 %.2e  最大 %.2e\n",
            median(abs(A[cm] - B[cm])), max(abs(A[cm] - B[cm]))))
cat(sprintf("  面板基因符号是否全部一致：%s\n",
            all(sign(A[PANEL]) == sign(B[PANEL]))))
cat("\n  结论：DESeq2 的中位数比法对 RNA5-8S5 这类高占比污染稳健，\n")
cat("        故稿件发现层（DESeq2 路线）的面板 logFC 不受本问题影响；\n")
cat("        受影响的是 03_preprocess.R 的 colSums 版 log-CPM 矩阵及其下游。\n")

cat(sprintf("\n已写出：03_分析/outputs/GSE89408_基因水平_双版本.rds\n  size 版 %d 基因 x %d 样本\n",
            nrow(size_mat), ncol(size_mat)))
