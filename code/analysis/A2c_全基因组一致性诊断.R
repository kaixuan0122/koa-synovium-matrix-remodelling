# A-2.3 复制层全基因组 logFC 一致性诊断
#
# 目的：把「GSE82107 低功效」与「GSE82107 真异质性」分开——**这是个待检验的命题，
# 不是既定解释**。方案 §八 明确要求现有解释不得被自动套用。
#
# 三个判据，都不依赖任何先验结论：
#   (1) 全基因组一致性：三个复制队列两两的 logFC 相关与同号率（取基因交集）。
#   (2) 置换零假设：打乱 GSE82107 的组别标签（保持其 n 与平台不变），重做全基因组检验，
#       得到「纯噪声下」的一致性分布。观测值若落在零分布内 → 该队列对这批基因没提供信息。
#   (3) 功效下界（MDE）：逐基因最小可检效应 = (z_{1-α/2}+z_{power}) × SE。
#       方向相左的基因分三类，不能一律写成「不提供信息」：
#         a. |logFC| > MDE            → 有功效仍反向，是异质性证据；
#         b. |logFC| ≤ MDE 但 P<0.05  → 标称显著却反向，须单独点出（P<0.05 不等于 80% 功效，
#                                        两者可同时成立：|t|≈2.1 落在 2.1–2.8 之间）；
#         c. |logFC| ≤ MDE 且 P≥0.05  → 测不出，不提供信息。
#       跨队列比 MDE 前必须换算到公共尺度（d = logFC / 后验 SD），
#       否则两队列残差 SD 不同（本例 GPL570 0.690 vs GPL96 0.512）会反转结论。
#
# **符号解析**：表格用旧名 CYR61/CTGF/MKL1，而 GPL96 平台上叫 CCN1/CCN2/MRTFA。
# 不做别名映射会把它们当成「该平台没测」而静默丢掉（本项目踩过），故下面显式解析并记录所用符号。
#
# 检验设定逐字照抄 B1_逐队列效应_力学基因.R，否则诊断与本表不可比：
#   limma 跑**完整矩阵**、design = ~0+group、contrast OA-Normal、eBayes、
#   SE 由 confint 反演、dfres = df.residual[1] + df.prior。
#
# 输入：04_分析数据/processed/验证层_genelevel.rds（工程根目录下）
# 输出：02_结果表/A2c_队列级诊断.csv、A2c_全基因组两两一致性.csv、
#       A2c_57基因_逐队列.csv、A2c_置换零分布.csv
#       03_分析/rds/A2c_全基因组拟合.rds
#       01_审计/A2c_全基因组一致性诊断_结果.md
suppressPackageStartupMessages({ library(limma) })

ROOT <- paste0("/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/",
               "Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/",
               "01.材料/01.博士材料/10.毕业/01.小论文/04.立项/",
               "滑膜基质刚度生成与刚度感知 0912")
REV  <- file.path(ROOT, "08_修订_2026-09-17")
RES  <- file.path(REV, "02_结果表"); AUD <- file.path(REV, "01_审计")
RDSD <- file.path(REV, "03_分析/rds")
for (d in c(RES, AUD, RDSD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

ALPHA <- 0.05; POWER <- 0.80; N_PERM <- 200L; SEED <- 20260923L
Z <- qnorm(1 - ALPHA / 2) + qnorm(POWER)          # 1.959964 + 0.8416212 = 2.801585

ln <- character(0)
rec <- function(s = "") { cat(s, "\n"); ln <<- c(ln, s) }

# 旧名 ↔ 新名互为别名。只做单向会漏一半，故两个方向都写。
ALIAS <- list(CTGF = "CCN2", CYR61 = "CCN1", MKL1 = "MRTFA", MRTFB = "MKL2",
              CCN2 = "CTGF", CCN1 = "CYR61", MRTFA = "MKL1", MKL2 = "MRTFB")

ov  <- read.csv(file.path(RES, "57基因总览.csv"), fileEncoding = "UTF-8", check.names = FALSE)
G57 <- ov$gene

VAL <- readRDS(file.path(ROOT, "04_分析数据/processed/验证层_genelevel.rds"))
stopifnot(setequal(names(VAL), c("GSE55235", "GSE55457", "GSE82107")))

## ---------- 0. 工具：检验函数 + 符号解析 ----------
fit_cohort <- function(mat, grp) {
  stopifnot(identical(sort(unique(grp)), c("Normal", "OA")))
  gf <- factor(ifelse(grp == "OA", "OA", "Normal"), levels = c("Normal", "OA"))
  design <- model.matrix(~0 + gf); colnames(design) <- c("Normal", "OA")
  fit <- eBayes(contrasts.fit(lmFit(mat, design), makeContrasts(OA - Normal, levels = design)))
  tt  <- topTable(fit, number = Inf, sort.by = "none", confint = TRUE)
  stopifnot(length(unique(fit$df.residual)) == 1L)
  dfres <- fit$df.residual[1] + fit$df.prior
  se <- (tt$CI.R - tt$CI.L) / (2 * qt(0.975, df = dfres))
  data.frame(gene = rownames(mat), logFC = tt$logFC, SE = se, P = tt$P.Value,
             s_post = sqrt(fit$s2.post), dfres = dfres, stringsAsFactors = FALSE)
}

# 把显示名解析到某队列实际使用的 rowname；解析不到返回 NA（**不静默丢弃，由调用方报出**）
resolve <- function(genes, rn) {
  sym <- rep(NA_character_, length(genes)); ali <- rep(FALSE, length(genes))
  for (i in seq_along(genes)) {
    if (genes[i] %in% rn) { sym[i] <- genes[i]; next }
    a <- ALIAS[[genes[i]]]
    if (!is.null(a) && a %in% rn) { sym[i] <- a; ali[i] <- TRUE }
  }
  list(sym = sym, aliased = ali)
}

## ---------- 1. 逐队列全基因组拟合 ----------
fits <- lapply(names(VAL), function(n) {
  o <- VAL[[n]]; f <- fit_cohort(o$mat, o$grp)
  f$cohort <- n
  f$n_OA <- sum(o$grp == "OA"); f$n_Normal <- sum(o$grp == "Normal")
  f$n_gene_platform <- nrow(o$mat)
  f
})
names(fits) <- names(VAL)
saveRDS(fits, file.path(RDSD, "A2c_全基因组拟合.rds"))

cohort_tab <- do.call(rbind, lapply(fits, function(f) data.frame(
  cohort = f$cohort[1], n_OA = f$n_OA[1], n_Normal = f$n_Normal[1],
  n_gene_platform = f$n_gene_platform[1],
  med_s_post = median(f$s_post, na.rm = TRUE), med_SE = median(f$SE, na.rm = TRUE),
  dfres = f$dfres[1],
  MDE_logFC_med = Z * median(f$SE, na.rm = TRUE),
  MDE_d_med = Z * median(f$SE / f$s_post, na.rm = TRUE),
  stringsAsFactors = FALSE)))
cohort_tab$n <- cohort_tab$n_OA + cohort_tab$n_Normal
cohort_tab <- cohort_tab[match(c("GSE55235", "GSE55457", "GSE82107"), cohort_tab$cohort), ]
rownames(cohort_tab) <- NULL

rec("# A-2.3 复制层全基因组 logFC 一致性诊断（实测）")
rec()
rec("> 生成自 `03_分析/scripts/A2c_全基因组一致性诊断.R`。检验设定与 `B1_逐队列效应_力学基因.R` 逐字相同")
rec("> （limma 完整矩阵 / `~0+group` / `contrast OA-Normal` / `eBayes` / SE 由 confint 反演）。")
rec()
rec("## 一、三个复制队列的基本量与功效下界")
rec()
rec(sprintf("逐基因最小可检效应 MDE = (z_{1-α/2} + z_power) × SE = %.3f × SE，α=%.2f、power=%.2f。",
            Z, ALPHA, POWER))
rec("MDE_d 为换算到公共尺度 d = logFC / 后验 SD 后的下界；跨队列比功效必须用 d 尺度。")
rec()
rec("| 队列 | n (OA/Normal) | 平台基因数 | 后验 SD 中位 | SE 中位 | 残差自由度 | MDE(logFC) 中位 | MDE(d) 中位 |")
rec("|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(cohort_tab))) { z <- cohort_tab[i, ]
  rec(sprintf("| %s | %d (%d/%d) | %d | %.3f | %.3f | %.1f | %.3f | %.2f |",
              z$cohort, z$n, z$n_OA, z$n_Normal, z$n_gene_platform,
              z$med_s_post, z$med_SE, z$dfres, z$MDE_logFC_med, z$MDE_d_med)) }
rec()

## ---------- 2. 全基因组两两一致性（取交集） ----------
as_map <- function(f) { v <- f$logFC; names(v) <- f$gene; v }
M <- lapply(fits, as_map)
pair_tab <- do.call(rbind, lapply(combn(names(M), 2, simplify = FALSE), function(p) {
  g <- intersect(names(M[[p[1]]]), names(M[[p[2]]]))
  x <- M[[p[1]]][g]; y <- M[[p[2]]][g]
  ok <- is.finite(x) & is.finite(y)
  data.frame(pair = paste(p[1], "vs", p[2]), n_common = length(g), n_used = sum(ok),
             pearson = cor(x[ok], y[ok]), spearman = cor(x[ok], y[ok], method = "spearman"),
             sign_agree_pct = 100 * mean(sign(x[ok]) == sign(y[ok])), stringsAsFactors = FALSE)
}))

rec("## 二、全基因组两两一致性（仅用两队列共有的基因）")
rec()
rec("| 队列对 | 交集基因 | 参与计算 | Pearson r | Spearman ρ | 同号率 |")
rec("|---|---|---|---|---|---|")
for (i in seq_len(nrow(pair_tab))) { z <- pair_tab[i, ]
  rec(sprintf("| %s | %d | %d | %.3f | %.3f | %.1f%% |",
              z$pair, z$n_common, z$n_used, z$pearson, z$spearman, z$sign_agree_pct)) }
rec()

## ---------- 3. 57 基因的符号解析（先解析，再谈一致性） ----------
REF_COHORTS <- c("GSE55235", "GSE55457")
res_ref <- lapply(fits[REF_COHORTS], function(f) resolve(G57, f$gene))
res_82  <- resolve(G57, fits[["GSE82107"]]$gene)
in_ref_both <- !is.na(res_ref[[1]]$sym) & !is.na(res_ref[[2]]$sym)
g57_use <- G57[in_ref_both]

symtab <- data.frame(gene = G57,
  GSE55235 = res_ref[[1]]$sym, GSE55457 = res_ref[[2]]$sym, GSE82107 = res_82$sym,
  stringsAsFactors = FALSE)
symtab$aliased <- symtab$gene != symtab$GSE82107
rec("## 三、57 基因在三个复制队列上的符号解析")
rec()
rec(sprintf("- 三个队列都能解析到：%d 个基因；两个 GPL96 队列都能解析到：%d 个。",
            sum(!is.na(res_82$sym) & in_ref_both), length(g57_use)))
miss96 <- G57[!in_ref_both]
rec(sprintf("- 在 GPL96 平台上**解析不到**（不是「无差异」，是该平台没测）：%d 个——%s。",
            length(miss96), paste(miss96, collapse = "、")))
al <- symtab[symtab$aliased & !is.na(symtab$GSE55235), ]
if (nrow(al)) {
  rec(sprintf("- 需别名才能匹配的：%s。",
              paste(sprintf("%s→%s", al$gene, al$GSE55235), collapse = "、")))
  rec("  （不做这一步会把这几个基因误判成「该平台没测」而静默丢掉。）")
}
# 别名匹配的唯一性核对：若某队列**同时**有旧名与新名两行，别名匹配就会在两者中二选一，
# 而任何结果都不会报错——必须显式查一遍（本项目踩过同类静默丢基因/取错行的坑）。
AMB <- do.call(rbind, lapply(names(fits), function(n) {
  rn <- rownames(VAL[[n]]$mat)
  do.call(rbind, lapply(names(ALIAS), function(a) {
    b <- ALIAS[[a]]
    if (a %in% rn && b %in% rn)
      data.frame(cohort = n, old = a, new = b, stringsAsFactors = FALSE) else NULL
  }))
}))
if (is.null(AMB)) {
  rec("- 唯一性核对：三个队列里**没有任何**旧名与新名并存的行，故别名匹配无二选一歧义。")
} else {
  rec("- ⚠️ 以下队列同时存在旧名与新名两行，别名匹配可能取错：")
  for (i in seq_len(nrow(AMB))) rec(sprintf("    %s：%s 与 %s 并存", AMB$cohort[i], AMB$old[i], AMB$new[i]))
  stop("别名匹配存在歧义，须先人工定向")
}
rec()

## ---------- 4. 置换零假设：打乱 GSE82107 标签 ----------
refs <- lapply(fits[REF_COHORTS], as_map)
g_ref <- intersect(names(refs[[1]]), names(refs[[2]]))
ref_lfc <- rowMeans(cbind(refs[[1]][g_ref], refs[[2]][g_ref]))     # 两队列均值作参照

# 57 基因子集上的参照（按各队列实际符号取）
r1 <- res_ref[[1]]$sym[in_ref_both]; r2 <- res_ref[[2]]$sym[in_ref_both]
ref57 <- (refs[[1]][r1] + refs[[2]][r2]) / 2
names(ref57) <- g57_use
sym82 <- res_82$sym[in_ref_both]; names(sym82) <- g57_use
stopifnot(!any(is.na(sym82)))            # GSE82107 应能解析全部 57 个；解析不到就不是「不可评」而是丢基因

g_eval <- intersect(g_ref, names(M[["GSE82107"]]))
obs82 <- M[["GSE82107"]]

set.seed(SEED)
o <- VAL[["GSE82107"]]
perm <- t(vapply(seq_len(N_PERM), function(i) {
  m <- as_map(fit_cohort(o$mat, sample(o$grp)))
  c(all = sum(sign(m[g_eval]) == sign(ref_lfc[g_eval]), na.rm = TRUE),
    g57 = sum(sign(m[sym82]) == sign(ref57), na.rm = TRUE),
    r_all = cor(m[g_eval], ref_lfc[g_eval]))
}, numeric(3)))
perm <- as.data.frame(perm)
n_all <- length(g_eval); n_g57 <- length(g57_use)
obs_all <- sum(sign(obs82[g_eval]) == sign(ref_lfc[g_eval]), na.rm = TRUE)
obs_57  <- sum(sign(obs82[sym82]) == sign(ref57), na.rm = TRUE)
p_emp <- function(obs, null) (1 + sum(null >= obs)) / (1 + length(null))

rec("## 四、置换零假设：GSE82107 的标签被打乱")
rec()
rec(sprintf("打乱 GSE82107 组别标签 %d 次（种子 %d），每次重做全基因组 limma，统计其 logFC 与", N_PERM, SEED))
rec(sprintf("「%s 两队列均值」同号的基因数。零分布即「82107 无真实组间差异」下的取值。",
            paste(REF_COHORTS, collapse = " + ")))
rec()
rec("| 基因集 | 可评基因 | 实测同号数 | 实测同号率 | 零分布中位 (范围) | 经验 P(≥实测) |")
rec("|---|---|---|---|---|---|")
for (spec in list(list(nm = "全基因组交集", n = n_all, obs = obs_all, null = perm$all),
                  list(nm = "57 基因子集", n = n_g57, obs = obs_57, null = perm$g57))) {
  rec(sprintf("| %s | %d | %d | %.1f%% | %.1f (%.0f–%.0f) | %.3f |", spec$nm, spec$n, spec$obs,
              100 * spec$obs / spec$n, median(spec$null), min(spec$null), max(spec$null),
              p_emp(spec$obs, spec$null))) }
rec()
rec(sprintf("同为参照的全基因组 Pearson r：实测 %.3f；零分布中位 %.3f（范围 %.3f–%.3f）。",
            cor(obs82[g_eval], ref_lfc[g_eval]), median(perm$r_all), min(perm$r_all), max(perm$r_all)))
rec()

## ---------- 5. 57 基因逐队列 + 方向相左的分类 ----------
per <- do.call(rbind, lapply(names(fits), function(n) {
  f <- fits[[n]]; rr <- resolve(G57, f$gene); k <- match(rr$sym, f$gene)
  data.frame(gene = G57, cohort = n, symbol_used = rr$sym, aliased = rr$aliased,
             logFC = f$logFC[k], SE = f$SE[k], P = f$P[k], s_post = f$s_post[k],
             stringsAsFactors = FALSE)
}))
per$measured <- !is.na(per$logFC)
per$MDE <- Z * per$SE
per$powered <- per$measured & abs(per$logFC) > per$MDE
per$nominal_sig <- per$measured & per$P < ALPHA
per$sign <- ifelse(per$measured, sign(round(per$logFC, 10)), NA_integer_)

dis <- per[per$cohort == "GSE82107" & per$gene %in% g57_use, ]
dis$ref_logFC <- ref57[dis$gene]
dis$ref_sign  <- sign(round(dis$ref_logFC, 10))
dis$discordant <- dis$sign != dis$ref_sign & dis$sign != 0 & dis$ref_sign != 0
# 三类，互斥且穷尽
dis$class <- ifelse(!dis$discordant, "方向一致",
             ifelse(dis$powered, "a_有功效仍反向",
             ifelse(dis$nominal_sig, "b_标称显著却反向", "c_测不出_不提供信息")))
disc <- dis[dis$discordant, ]

rec("## 五、GSE82107 与另两个复制队列方向相左的基因")
rec()
rec(sprintf("参照 = %s 两队列 logFC 的均值；可评子集 = 两个 GPL96 队列都能解析的 %d 个基因。",
            paste(REF_COHORTS, collapse = " + "), length(g57_use)))
rec()
tb <- table(factor(dis$class, levels = c("a_有功效仍反向", "b_标称显著却反向", "c_测不出_不提供信息", "方向一致")))
rec(sprintf("- 方向相左共 **%d 个**（占可评 %d 的 %.0f%%）：有功效仍反向 %d 个；标称显著却反向 %d 个；测不出 %d 个。",
            nrow(disc), nrow(dis), 100 * nrow(disc) / nrow(dis), tb[1], tb[2], tb[3]))
rec()
rec("| 基因 | 参照 logFC | 82107 logFC | 82107 SE | 82107 MDE | 82107 P | 82107 标称显著? | 分类 |")
rec("|---|---|---|---|---|---|---|---|")
for (i in order(dis$P)[seq_len(nrow(dis))]) { z <- dis[i, ]
  if (!z$discordant) next
  rec(sprintf("| %s | %+.3f | %+.3f | %.3f | %.3f | %s | %s | %s |", z$gene,
              z$ref_logFC, z$logFC, z$SE, z$MDE,
              ifelse(is.na(z$P), "—", formatC(z$P, format = "g", digits = 3)),
              ifelse(z$nominal_sig, "是", "否"),
              sub("^[abc]_", "", z$class))) }
rec()

# 全基因组背景：同样的分类占比，用来判断 57 基因是否特殊
bg <- data.frame(g = g_ref, ref = ref_lfc, q = obs82[g_ref],
                 se = fits[["GSE82107"]]$SE[match(g_ref, names(obs82))],
                 p = fits[["GSE82107"]]$P[match(g_ref, names(obs82))])
bg <- bg[is.finite(bg$ref) & is.finite(bg$q), ]
bg$disc <- sign(round(bg$q, 10)) != sign(round(bg$ref, 10)) & sign(round(bg$ref, 10)) != 0
bg$powered <- abs(bg$q) > Z * bg$se
bg$nominal_sig <- bg$p < ALPHA
rec(sprintf("全基因组背景（可评 %d 基因）：方向相左 %.0f%%；其中「有功效仍反向」占全部 %.0f%%、「标称显著却反向」%.0f%%。",
            nrow(bg), 100 * mean(bg$disc), 100 * mean(bg$disc & bg$powered),
            100 * mean(bg$disc & bg$nominal_sig & !bg$powered)))
rec(sprintf("57 基因子集（可评 %d 基因）：方向相左 %.0f%%；与零分布（%.1f 同号）相比经验 P = %.3f。",
            nrow(dis), 100 * nrow(disc) / nrow(dis), median(perm$g57), p_emp(obs_57, perm$g57)))
rec()

## ---------- 6. 落盘 ----------
write.csv(cohort_tab, file.path(RES, "A2c_队列级诊断.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(pair_tab,   file.path(RES, "A2c_全基因组两两一致性.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(per,        file.path(RES, "A2c_57基因_逐队列.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(symtab,     file.path(RES, "A2c_57基因_符号解析.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(dis[, c("gene", "ref_logFC", "logFC", "SE", "MDE", "P", "powered",
                  "nominal_sig", "discordant", "class")],
          file.path(RES, "A2c_82107_方向相左分类.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(perm[, c("all", "g57", "r_all")], file.path(RES, "A2c_置换零分布.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
writeLines(c(ln, "", "---",
  sprintf("置换 %d 次，种子 %d。检验设定与 B1 逐字相同。", N_PERM, SEED),
  sprintf("符号解析：GPL96 平台上解析不到 %d 个（%s），非「无差异」。",
          length(miss96), paste(miss96, collapse = "、")),
  "**结论只写进 §五 的实测句，未预置任何解释。**"),
  file.path(AUD, "A2c_全基因组一致性诊断_结果.md"))
cat("\n已写出 02_结果表/A2c_*.csv（6 份）、03_分析/rds/A2c_全基因组拟合.rds、01_审计/A2c_全基因组一致性诊断_结果.md\n")
