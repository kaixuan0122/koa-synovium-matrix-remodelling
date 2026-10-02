# B1 逐队列效应量（§4.1）：57 个力学来源基因 × 6 个队列，给出 logFC / SE / 95%CI / P / adjP
#
# 为什么补这一步：现有产物只对 6 个面板基因存过 logFC，且发现层
# （主分析_6刚度基因_meta.csv）根本没有 SE/CI。§4.1 要求给全部力学基因一个概览，
# §4.2 的森林图与留一法又都需要逐队列 SE 作为输入。本脚本一次性补齐。
#
# 设定与既有脚本严格一致，以便逐值回归：
#   * 微阵列：limma，与 13_replication_per_cohort.R 同设定。
#     **limma 必须跑在完整矩阵上**——adj.P.Val 是跨该队列全部基因做 BH 得到的，
#     若先子集到 65 个基因再跑，adjP 会整体变小，与表 4 对不上。
#   * GSE89408：DESeq2，与 05_de_meta.R 同设定（count 按基因求和 → round → ~group）。
#
# 符号别名：基因集沿用旧符号（CTGF / CYR61 / MKL1），而 GSE206848 / GSE55235 /
#   GSE55457 的注释表里只有新符号（CCN2 / CCN1 / MRTFA）。不显式解析别名，
#   %in% 会把它们静默判成「该平台未检出」，写出事实错误的结论。
#   本脚本对每个基因逐个队列尝试 [原名 → 别名]，并把实际取用的符号落盘备查。
#
# 输入：04_分析数据/processed/{主分析,验证层}_genelevel.rds
#       04_分析数据/raw/GSE89408/GSE89408_GEO_count_matrix_rename.txt.gz
#       04_分析数据/力学基因集/力学来源基因集.csv
# 输出：02_结果表/逐队列效应_力学基因.csv      逐队列效应量（含 SE）
#       02_结果表/基因符号映射核对.csv         每个基因在每个队列实际用到的符号
suppressPackageStartupMessages({
  library(Biobase); library(limma); library(DESeq2)
})

BASE <- "/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/01.材料/01.博士材料/10.毕业/01.小论文/04.立项/滑膜基质刚度生成与刚度感知 0912"
P    <- file.path(BASE, "04_分析数据/processed")
RAW  <- file.path(BASE, "04_分析数据/raw")
OUT  <- file.path(BASE, "08_修订_2026-09-17/02_结果表")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ROLE <- c(GSE89408 = "发现", GSE206848 = "发现", GSE1919 = "发现",
          GSE55235 = "复制", GSE55457 = "复制", GSE82107 = "复制")

# 别名表：基因集里的旧符号 → 平台上可能出现的新符号
ALIAS <- list(CTGF = "CCN2", CYR61 = "CCN1", MKL1 = "MRTFA", MRTFB = "MKL2")

## ---------- 1. 基因清单 ----------
gs <- read.csv(file.path(BASE, "04_分析数据/力学基因集/力学来源基因集.csv"),
               fileEncoding = "UTF-8", stringsAsFactors = FALSE)
GENES <- gs$Symbol
CROSSLINK <- c("LOX","LOXL1","LOXL2","LOXL3","LOXL4","PLOD1","PLOD2","TGM2")
cat(sprintf("基因集：%d 个基因；其中刚度生成-交联 %d 个\n", length(GENES), length(CROSSLINK)))
cat(sprintf("交联基因是否都在基因集内：%s\n", all(CROSSLINK %in% GENES)))

## ---------- 2. 微阵列：limma 跑完整矩阵 ----------
main <- readRDS(file.path(P, "主分析_genelevel.rds"))
valid <- readRDS(file.path(P, "验证层_genelevel.rds"))
# 注意：主分析_genelevel.rds 里**同时含 GSE89408**（03_preprocess.R 第 83 行三个元素一起存）。
# 若直接 c(main, valid)，GSE89408 会既被当微阵列跑一遍 limma、又被下面的 DESeq2 跑一遍，
# 同一基因产出两行 logFC（两者相差可达 0.6），回归测试的 max|Δ| 取到 limma 那一行即报不一致。
MATS <- c(main[setdiff(names(main), "GSE89408")], valid)
stopifnot(length(MATS) == 5L, !"GSE89408" %in% names(MATS))
cat("微阵列队列（limma）:", paste(names(MATS), collapse = ", "), "\n")

resolve_sym <- function(gene, avail) {
  # 返回实际可在 avail 中取到的符号；原符号优先，其次别名
  if (gene %in% avail) return(gene)
  alt <- ALIAS[[gene]]
  if (!is.null(alt) && alt %in% avail) return(alt)
  NA_character_
}

rows <- list(); map_rows <- list()
for (nm in names(MATS)) {
  m <- MATS[[nm]]$mat; g <- MATS[[nm]]$grp
  gf <- factor(ifelse(g == "OA", "OA", "Normal"), levels = c("Normal", "OA"))
  design <- model.matrix(~0 + gf); colnames(design) <- c("Normal", "OA")
  fit <- eBayes(contrasts.fit(lmFit(m, design), makeContrasts(OA - Normal, levels = design)))
  # confint=TRUE 才给 CI.L/CI.R；默认 topTable 无 SE 列，取 r$SE 会得 NULL 并静默退化
  tt <- topTable(fit, number = Inf, sort.by = "none", confint = TRUE)
  # df.residual 是长度为 nrow(m) 的向量。若直接参与算术，SE 会被广播成上万行，
  # data.frame() 随之把单基因的一行撑成 nrow(m) 行——本脚本第一版就栽在这里，
  # 产出 563 万行且不报任何错。故先断言其恒定，再取标量。
  dfres <- fit$df.residual[1] + fit$df.prior
  stopifnot(length(unique(fit$df.residual)) == 1L)
  q975 <- qt(0.975, df = dfres)
  cat(sprintf("\n=== %s (%s, %s): %d 基因 x %d 样本 (OA %d / Normal %d), 残差 df=%.1f ===\n",
              nm, ROLE[[nm]], "microarray", nrow(m), ncol(m),
              sum(g == "OA"), sum(g == "Normal"), dfres))
  n_miss <- 0
  for (gene in GENES) {
    s <- resolve_sym(gene, rownames(m))
    map_rows[[length(map_rows) + 1]] <- data.frame(
      dataset = nm, gene_set_symbol = gene,
      symbol_used = if (is.na(s)) NA_character_ else s,
      resolved_by = if (is.na(s)) "未检出" else if (s == gene) "原名" else "别名",
      stringsAsFactors = FALSE)
    if (is.na(s)) { n_miss <- n_miss + 1; next }
    r <- tt[s, ]
    se <- (r$CI.R - r$CI.L) / (2 * q975)
    rows[[length(rows) + 1]] <- data.frame(
      dataset = nm, role = ROLE[[nm]], gene = gene, symbol_used = s,
      n_OA = sum(g == "OA"), n_Normal = sum(g == "Normal"),
      logFC = r$logFC, SE = se, CI_low = r$CI.L, CI_high = r$CI.R,
      P = r$P.Value, adjP = r$adj.P.Val,
      test = "limma", stringsAsFactors = FALSE)
  }
  cat(sprintf("  基因集内未检出 %d / %d\n", n_miss, length(GENES)))
}

## ---------- 3. GSE89408：DESeq2 跑完整 count 矩阵 ----------
cnt <- read.delim(gzfile(file.path(RAW, "GSE89408/GSE89408_GEO_count_matrix_rename.txt.gz")),
                  check.names = FALSE, stringsAsFactors = FALSE)
genes_id <- sub(":.*", "", cnt[, 1])
M <- as.matrix(cnt[, -1])
oa_cols <- grep("^OA_tissue_",  colnames(M), value = TRUE)
n_cols  <- grep("^normal_tissue_", colnames(M), value = TRUE)
M <- M[, c(oa_cols, n_cols), drop = FALSE]
Magg <- aggregate(M, by = list(genes_id), FUN = sum)   # count 可加，按基因求和
rownames(Magg) <- Magg[, 1]; Magg <- as.matrix(Magg[, -1, drop = FALSE])
Mround <- round(Magg)
grp <- factor(c(rep("OA", length(oa_cols)), rep("Normal", length(n_cols))),
              levels = c("Normal", "OA"))
dds <- DESeqDataSetFromMatrix(countData = Mround, colData = data.frame(group = grp),
                              design = ~group)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("group", "OA", "Normal"))
stopifnot("lfcSE" %in% colnames(res))   # 不给 SE 就说明取错了对象，宁可报错不静默
res <- as.data.frame(res); res$symbol <- rownames(res)
# 自洽性：Wald 统计量须等于 logFC/lfcSE。若不等，说明 logFC 与 lfcSE 不同源
# （例如一个来自收缩、一个来自原始拟合），此时二者的比值与 P 都对不上。
stat_gap <- max(abs(res$stat - res$log2FoldChange / res$lfcSE), na.rm = TRUE)
cat(sprintf("  DESeq2 Wald 自洽 max|stat − logFC/SE| = %.3g\n", stat_gap))
stopifnot(stat_gap < 1e-8)
cat(sprintf("\n=== GSE89408 (发现, RNA-seq DESeq2): %d 基因 x %d 样本 (OA %d / Normal %d) ===\n",
            nrow(res), length(grp), sum(grp == "OA"), sum(grp == "Normal")))
n_miss <- 0
for (gene in GENES) {
  s <- resolve_sym(gene, res$symbol)
  map_rows[[length(map_rows) + 1]] <- data.frame(
    dataset = "GSE89408", gene_set_symbol = gene,
    symbol_used = if (is.na(s)) NA_character_ else s,
    resolved_by = if (is.na(s)) "未检出" else if (s == gene) "原名" else "别名",
    stringsAsFactors = FALSE)
  if (is.na(s)) { n_miss <- n_miss + 1; next }
  r <- res[match(s, res$symbol), ]
  rows[[length(rows) + 1]] <- data.frame(
    dataset = "GSE89408", role = "发现", gene = gene, symbol_used = s,
    n_OA = sum(grp == "OA"), n_Normal = sum(grp == "Normal"),
    logFC = r$log2FoldChange, SE = r$lfcSE,
    CI_low = r$log2FoldChange - qnorm(0.975) * r$lfcSE,
    CI_high = r$log2FoldChange + qnorm(0.975) * r$lfcSE,
    P = r$pvalue, adjP = r$padj,
    test = "DESeq2", stringsAsFactors = FALSE)
}
cat(sprintf("  基因集内未检出 %d / %d\n", n_miss, length(GENES)))

## ---------- 4. 落盘 ----------
eff <- do.call(rbind, rows)
# 结构断言：每个 (队列, 基因) 至多一行。少了说明有基因被静默丢弃，多了说明又发生了广播。
stopifnot(nrow(eff) <= length(GENES) * 6L)
stopifnot(!any(duplicated(paste(eff$dataset, eff$gene))))
stopifnot(all(is.finite(eff$SE)), all(eff$SE > 0))
cat(sprintf("\n表中 %d 行（上限 %d），无重复键，SE 全部有限且为正\n",
            nrow(eff), length(GENES) * 6L))
eff <- eff[order(eff$gene, factor(eff$dataset, levels = c(names(MATS), "GSE89408"))), ]
write.csv(eff, file.path(OUT, "逐队列效应_力学基因.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

symmap <- do.call(rbind, map_rows)
write.csv(symmap, file.path(OUT, "基因符号映射核对.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## ---------- 5. 回归测试：与既有产物逐值比对 ----------
cat("\n========== 回归测试：与既有产物比对 ==========\n")
ok <- TRUE
chk <- function(lbl, a, b, tol = 1e-8) {
  # 必须都是标量：若任一侧是向量，max() 只取最差的那个，会掩盖其余的不一致
  stopifnot(length(a) == 1L, length(b) == 1L)
  d <- abs(a - b)
  flag <- if (is.finite(d) && d < tol) "✅" else "❌"
  if (flag == "❌") ok <<- FALSE
  cat(sprintf("  %s %-46s |Δ| = %.3g\n", flag, lbl, d))
}

# SE 的对账基准不是朴素 OLS：topTable(confint=TRUE) 的 CI 用的是 eBayes 收缩后的方差，
# 故 SE = stdev.unscaled × sqrt(s2.post)，与 OLS 的 sigma 本来就不同（GSE55235/LOXL1 相差
# 0.04，方向合理）。对账要验的是「CI 折回 SE」这一反演是否成立，而不是它与 OLS 是否相等——
# 后者会被误读成数值错误，实际只是选错了基准。
nm_chk <- "GSE55235"; g_chk <- "LOXL1"
m  <- MATS[[nm_chk]]$mat
gf <- factor(ifelse(MATS[[nm_chk]]$grp == "OA", "OA", "Normal"), levels = c("Normal", "OA"))
des2 <- model.matrix(~0 + gf); colnames(des2) <- c("Normal", "OA")
fit2 <- eBayes(contrasts.fit(lmFit(m, des2), makeContrasts(OA - Normal, levels = des2)))
se_mod <- unname(fit2$stdev.unscaled[g_chk, 1] * sqrt(fit2$s2.post[g_chk]))
se_lim <- eff$SE[eff$dataset == nm_chk & eff$gene == g_chk]
chk(sprintf("SE 反演对账 %s %s（moderated）", nm_chk, g_chk), se_lim, se_mod)
cat(sprintf("     （参考：朴素 OLS 的 SE = %.4f，与 moderated SE 相差 %.4f，属正常）\n",
            summary(lm(m[g_chk, ] ~ gf))$coefficients["gfOA", "Std. Error"],
            abs(se_mod - summary(lm(m[g_chk, ] ~ gf))$coefficients["gfOA", "Std. Error"])))
# CI 与 P 必须来自同一次检验：CI 不含 0 ⟺ P < 0.05（同一 df 的 t 检验）
k <- eff$dataset == nm_chk & eff$gene == g_chk
stopifnot((eff$CI_low[k] > 0 | eff$CI_high[k] < 0) == (eff$P[k] < 0.05))
cat("     ✅ GSE55235/LOXL1：CI 是否含 0 与 P<0.05 判定一致\n")

# 5a. 复制层面板基因 ↔ 复制层_逐队列_面板基因.csv（limma 同设定）
old <- read.csv(file.path(P, "复制层_逐队列_面板基因.csv"), stringsAsFactors = FALSE)
for (i in seq_len(nrow(old))) {
  k <- eff$dataset == old$cohort[i] & eff$gene == old$gene[i]
  if (!any(k)) { cat(sprintf("  ❌ 缺 %s / %s\n", old$cohort[i], old$gene[i])); ok <- FALSE; next }
  chk(sprintf("%s %s logFC", old$cohort[i], old$gene[i]),
      eff$logFC[k], old$logFC[i])
  chk(sprintf("%s %s P", old$cohort[i], old$gene[i]), eff$P[k], old$P[i])
}

# 5b. 发现层面板基因 ↔ 主分析_6刚度基因_meta.csv（89408 走 DESeq2，其余走 limma）
old2 <- read.csv(file.path(P, "主分析_6刚度基因_meta.csv"), stringsAsFactors = FALSE)
for (i in seq_len(nrow(old2))) {
  for (ds in c("GSE89408", "GSE206848", "GSE1919")) {
    col <- paste0("logFC_", sub("GSE", "", ds))
    v <- old2[[col]][i]
    if (is.na(v)) next
    k <- eff$dataset == ds & eff$gene == old2$gene[i]
    if (!any(k)) { cat(sprintf("  ❌ 缺 %s / %s\n", ds, old2$gene[i])); ok <- FALSE; next }
    chk(sprintf("%s %s logFC", ds, old2$gene[i]), eff$logFC[k], v)
  }
}
cat(if (ok) "\n✅ 全部回归测试通过——本次逐队列数与既有产物一致\n" else "\n❌ 存在不一致，不得继续\n")

## ---------- 6. 「未检出」原因核查 ----------
# 「未检出」有两种性质完全不同的含义，不得混写：
#   (a) 平台根本没有该基因的探针 —— 平台限制，可写进局限；
#   (b) 平台有探针但注释表未给符号 —— 注释缺口，可换注释来源补救。
# 判据：直接查平台注释表里是否存在该符号（原名或别名）。
suppressPackageStartupMessages({
  library(AnnotationDbi); library(hgu133a.db); library(hgu133plus2.db); library(GEOquery)
})
miss <- symmap[symmap$resolved_by == "未检出", ]
if (nrow(miss)) {
  # 平台 → 注释来源
  tab91 <- Table(getGEO(filename = file.path(RAW, "GPL91.soft.gz")))
  sym91 <- unique(sub("///.*", "", tab91$`Gene Symbol`))
  sym91 <- sym91[!is.na(sym91) & sym91 != "" & sym91 != "---"]
  DB <- list(GSE55235 = hgu133a.db, GSE55457 = hgu133a.db,
             GSE206848 = hgu133plus2.db, GSE82107 = hgu133plus2.db, GSE1919 = NULL)
  has_sym <- function(ds, s) {
    cand <- unique(c(s, unlist(ALIAS[[s]])))
    if (ds == "GSE1919") return(any(cand %in% sym91))
    db <- DB[[ds]]
    any(!is.na(suppressWarnings(
      AnnotationDbi::select(db, keys = cand, keytype = "SYMBOL", columns = "PROBEID")$PROBEID)))
  }
  miss$平台有该符号的探针 <- mapply(has_sym, miss$dataset, miss$gene_set_symbol)
  miss$结论 <- ifelse(miss$`平台有该符号的探针`, "注释缺口（可换注释来源补救）", "平台无探针（平台限制）")
  write.csv(miss, file.path(OUT, "未检出原因核查.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("\n========== 「未检出」原因核查 ==========\n")
  print(miss[, c("dataset", "gene_set_symbol", "结论")], row.names = FALSE)
  cat(sprintf("  平台无探针 %d 条；注释缺口 %d 条\n",
              sum(!miss$`平台有该符号的探针`), sum(miss$`平台有该符号的探针`)))
}

cat(sprintf("\n已写出：\n  02_结果表/逐队列效应_力学基因.csv  %d 行\n  02_结果表/基因符号映射核对.csv     %d 行\n",
            nrow(eff), nrow(symmap)))
