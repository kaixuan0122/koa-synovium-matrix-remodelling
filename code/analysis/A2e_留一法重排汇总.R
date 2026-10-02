# A2e 留一法重排汇总（供正文 §3.2「队列依赖性」引用）
#
# 目的：把留一法转成**可直接写进正文的重排名次**。本脚本只做两件事：
#   (1) 发现层：用**该层主估计**（IVW-FE Meta）重做留一法并排名次；
#   (2) 复制层：直接用既有 `留一法_合并效应.csv`（B2 产出）排名次——因为该层主估计
#       就是 pooled(~0+group+cohort)，与留一法用的是同一个估计量，无需重算。
#
# 为什么必须分成两条路径（这是本脚本存在的理由）：
#   B2 的 `留一法_合并效应.csv` 全部由 `fit_pooled` 产出（`lm(y ~ 0+group+cohort)` + HC3）。
#   该估计量在**复制层**是主估计，但在**发现层**不是——发现层主估计是 IVW-FE，
#   而 pooled 在发现层要跨 RNA-seq 与两套微阵列平台池化**原始表达值**，
#   仅靠 cohort 固定效应只能吸收加性平台位移、不能吸收尺度差异。
#   实测后果（本脚本的断言之一）：发现层 LOXL1 的 pooled=+0.0101（P=0.976），
#   而其逐队列 DESeq2/limma 估计的 IVW-FE=+0.6382（P=3.9e-4）；
#   逐队列均值差为 1919=+0.585、206848=+0.638、89408=**−0.285**，
#   即 pooled 与 B1 的 DESeq2 估计在 GSE89408 上**符号相反**（B1 为 +0.742）。
#   若拿这份 pooled 留一法去谈「发现层的队列依赖性」，报出的将是一个
#   不属于该层主估计的性质——正是记忆 sensitivity-script-must-copy-primary-settings-exactly
#   所禁止的（敏感性分析必须逐字照抄主分析的设定）。
#
# 判据：与 A2_交联家族_判据敏感性.R 一致，名次两套——按 |logFC| 降序、按 P 升序，
#   「1 = 最靠前」。名次一律在**该层全体可评基因**内计算，不在家族子集内计算
#   （记忆 rank-within-filtered-subset-always-yields-1-to-n：先过滤再排名会造出构造性假象）。
#
# 输入：02_结果表/逐队列效应_力学基因.csv（B1）、02_结果表/合并效应_{发现,复制}层.csv（B2）、
#       02_结果表/留一法_合并效应.csv（B2）
# 输出：02_结果表/A2e_留一法重排汇总.csv、02_结果表/A2e_留一法_逐基因名次变动.csv

suppressPackageStartupMessages(library(metafor))

ROOT <- "/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/01.材料/01.博士材料/10.毕业/01.小论文/04.立项/滑膜基质刚度生成与刚度感知 0912/08_修订_2026-09-17"
TAB <- file.path(ROOT, "02_结果表")

PANEL <- c("LOXL1", "VGLL4", "SDC4", "SERPINE1", "ITGA5", "SAV1")
CROSSLINK <- c("LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4", "PLOD1", "PLOD2", "TGM2")
DISC <- c("GSE89408", "GSE206848", "GSE1919")
REP  <- c("GSE55235", "GSE55457", "GSE82107")

eff <- read.csv(file.path(TAB, "逐队列效应_力学基因.csv"), fileEncoding = "UTF-8",
                check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("dataset", "role", "gene", "logFC", "SE", "P") %in% names(eff)))
eff <- eff[is.finite(eff$logFC) & is.finite(eff$SE) & eff$SE > 0, ]

ivw_fe <- function(gene, cohorts) {          # 与 B2 的 meta_ivw 同一设定：metafor、method="FE"
  s <- eff[eff$gene == gene & eff$dataset %in% cohorts, ]
  if (nrow(s) < 2) return(NULL)
  r <- try(metafor::rma(yi = s$logFC, sei = s$SE, method = "FE"), silent = TRUE)
  if (inherits(r, "try-error")) return(NULL)
  data.frame(logFC = as.numeric(r$beta), P = as.numeric(r$pval), k = nrow(s))
}

## ---------- 1. 发现层：IVW-FE 留一法（自洽）；复制层：读既有 pooled 留一法 ----------
loo_new <- list()
for (L in c("发现", "复制")) {
  cohorts <- if (L == "发现") DISC else REP
  genes <- unique(eff$gene[eff$role == L])
  for (g in genes) {
    full <- ivw_fe(g, cohorts)
    if (is.null(full) || full$k < 3) next
    for (drop in cohorts) {
      f2 <- ivw_fe(g, setdiff(cohorts, drop))
      if (is.null(f2)) next
      loo_new[[length(loo_new) + 1]] <- data.frame(
        layer = L, gene = g, dropped = drop, k = f2$k,
        logFC = f2$logFC, P = f2$P, estimator = "IVW-FE", stringsAsFactors = FALSE)
    }
  }
}
LOO_IVW <- do.call(rbind, loo_new)

full_ivw <- do.call(rbind, lapply(c("发现", "复制"), function(L) {
  cohorts <- if (L == "发现") DISC else REP
  gs <- unique(eff$gene[eff$role == L])
  d <- do.call(rbind, lapply(gs, function(g) {
    f <- ivw_fe(g, cohorts); if (is.null(f)) return(NULL)
    data.frame(layer = L, gene = g, dropped = "（全量）", k = f$k, logFC = f$logFC, P = f$P,
               estimator = "IVW-FE", stringsAsFactors = FALSE)
  }))
  d
}))

## 复制层主估计是 pooled，故其留一法与全量参照都用既有 B2 产物（同一估计量）
LOO_POOL <- read.csv(file.path(TAB, "留一法_合并效应.csv"), fileEncoding = "UTF-8",
                     check.names = FALSE, stringsAsFactors = FALSE)
LOO_POOL <- LOO_POOL[LOO_POOL$layer == "复制", ]
LOO_POOL <- LOO_POOL[is.finite(LOO_POOL$logFC) & is.finite(LOO_POOL$P) & LOO_POOL$P > 0,
                     c("layer", "gene", "dropped", "k", "logFC", "P")]
LOO_POOL$estimator <- "pooled(~0+group+cohort)"

full_pool <- read.csv(file.path(TAB, "合并效应_复制层.csv"), fileEncoding = "UTF-8",
                      check.names = FALSE, stringsAsFactors = FALSE)
full_pool <- full_pool[is.finite(full_pool$pooled_logFC) & is.finite(full_pool$pooled_P), ]
full_pool <- data.frame(layer = "复制", gene = full_pool$gene, dropped = "（全量）", k = full_pool$k,
                        logFC = full_pool$pooled_logFC, P = full_pool$pooled_P,
                        estimator = "pooled(~0+group+cohort)", stringsAsFactors = FALSE)

## ---------- 2. 对账：本脚本重算的 IVW-FE 全量必须复现 B2 的 ivw_* 列 ----------
chk <- read.csv(file.path(TAB, "合并效应_发现层.csv"), fileEncoding = "UTF-8",
                check.names = FALSE, stringsAsFactors = FALSE)
chk <- chk[is.finite(chk$ivw_logFC), ]
m <- merge(full_ivw[full_ivw$layer == "发现", c("gene", "logFC", "P")],
           chk[, c("gene", "ivw_logFC", "ivw_P")], by = "gene")
d_b <- max(abs(m$logFC - m$ivw_logFC)); d_p <- max(abs(m$P - m$ivw_P))
cat(sprintf("\n对账：本脚本重算 IVW-FE vs B2 的 ivw 列 —— n=%d，logFC 最大绝对差 %.3g，P 最大绝对差 %.3g\n",
            nrow(m), d_b, d_p))
stopifnot(nrow(m) == nrow(chk), d_b < 1e-9, d_p < 1e-12)
cat("对账通过：发现层留一法所用估计量与 B2 主估计逐位一致（非近似）。\n")

## ---------- 3. 暴露既有 pooled 留一法在发现层的不可用（记入日志，供审计引用）----------
pf <- read.csv(file.path(TAB, "合并效应_发现层.csv"), fileEncoding = "UTF-8",
               check.names = FALSE, stringsAsFactors = FALSE)
diag <- merge(pf[is.finite(pf$pooled_logFC), c("gene", "pooled_logFC", "pooled_P", "ivw_logFC", "ivw_P", "cohort_logmeandiff")],
              full_ivw[full_ivw$layer == "发现", c("gene", "logFC")], by = "gene")
diag$diff <- diag$pooled_logFC - diag$logFC
diag <- diag[order(-abs(diag$diff)), ]
cat("\n发现层 pooled 与其主估计 IVW-FE 分歧最大的 5 个基因（示 pooled 不能当发现层主估计用）：\n")
print(head(diag[, c("gene", "pooled_logFC", "pooled_P", "ivw_logFC", "ivw_P", "diff")], 5), row.names = FALSE)
cat(sprintf("发现层 sign(pooled) != sign(IVW-FE) 的基因数：%d / %d\n",
            sum(sign(diag$pooled_logFC) != sign(diag$logFC)), nrow(diag)))

## ---------- 4. 名次 ----------
# 名次一律在 **(层, 估计量)** 内计算。复制层有两个估计量（IVW-FE 与 pooled），
# 二者必须各自成集比较——混在一个集里排名会把「换估计量」的效果当成「换队列」的效果。
KEY <- c("layer", "estimator")
rank_set <- function(df) transform(df, rank_abs = rank(-abs(logFC), ties.method = "min"),
                                       rank_P   = rank(P,         ties.method = "min"))
LOO  <- rbind(LOO_IVW, LOO_POOL)
FULL <- rbind(full_ivw, full_pool)
FULL$grp <- do.call(paste, c(FULL[KEY], sep = " | "))
LOO$grp  <- do.call(paste, c(LOO[KEY],  sep = " | "))
# 排名单位是 **(层, 估计量, 剔除队列)**：一次留一只在一个剔除下排名；
# 若不带 dropped，三个剔除的 156 行会被排成 1..156，名次全部失真。
agg <- function(df) do.call(rbind, lapply(split(df, paste(df$grp, df$dropped, sep = " || ")), function(z) {
  # 不变量：同一个 (基因, 剔除队列) 在每个排名集里只能出现一次——若两个估计量又被并到一起，
  # 这里立刻失败（这正是上一版 n_set=104 那个 bug 的形态）。
  stopifnot(length(unique(z$gene)) == nrow(z))
  z <- rank_set(z); z$n_set <- nrow(z); z
}))
LOO <- agg(LOO); FULL <- agg(FULL)
sz <- vapply(split(FULL$n_set, FULL$grp), function(z) z[1], numeric(1))
print(sz)
stopifnot(sz[["发现 | IVW-FE"]] == 57L,
          sz[["复制 | pooled(~0+group+cohort)"]] == 52L)  # 复制层主估计的排名集必须是 52，不是 104

cmp <- merge(LOO[, c(KEY, "gene", "dropped", "k", "logFC", "P", "rank_abs", "rank_P")],
             FULL[, c(KEY, "gene", "rank_abs", "rank_P")],
             by = c(KEY, "gene"), suffixes = c("_loo", "_full"))
cmp$grp <- do.call(paste, c(cmp[KEY], sep = " | "))
cmp$d_abs <- cmp$rank_abs_loo - cmp$rank_abs_full
cmp$d_P   <- cmp$rank_P_loo   - cmp$rank_P_full
cmp$in_panel <- cmp$gene %in% PANEL
cmp$is_crosslink <- cmp$gene %in% CROSSLINK
cmp <- cmp[order(cmp$layer, cmp$gene, cmp$dropped), ]
write.csv(cmp, file.path(TAB, "A2e_留一法_逐基因名次变动.csv"), row.names = FALSE, fileEncoding = "UTF-8")

summ <- list()
for (G in unique(FULL$grp)) {
  f <- FULL[FULL$grp == G, ]
  z_all <- LOO[LOO$grp == G, ]
  base_abs <- f$gene[f$rank_abs <= 10]; base_P <- f$gene[f$rank_P <= 10]
  cf <- cmp[cmp$grp == G, ]
  for (dd in unique(z_all$dropped)) {
    z <- z_all[z_all$dropped == dd, ]
    a <- z$gene[z$rank_abs <= 10]; p <- z$gene[z$rank_P <= 10]
    fl <- sign(z$logFC) != sign(f$logFC[match(z$gene, f$gene)])
    summ[[length(summ) + 1]] <- data.frame(
      layer = f$layer[1], estimator = f$estimator[1], dropped = dd, n_set = nrow(z),
      top10_abs_overlap = length(intersect(a, base_abs)),
      top10_P_overlap   = length(intersect(p, base_P)),
      max_abs_shift_crosslink = max(abs(cf$d_abs[cf$dropped == dd & cf$is_crosslink]), na.rm = TRUE),
      max_P_shift_crosslink   = max(abs(cf$d_P  [cf$dropped == dd & cf$is_crosslink]), na.rm = TRUE),
      max_abs_shift_panel     = max(abs(cf$d_abs[cf$dropped == dd & cf$in_panel]), na.rm = TRUE),
      n_sign_flip = sum(fl),
      flips = paste(z$gene[fl], collapse = "、"),
      stringsAsFactors = FALSE)
  }
}
SUMM <- do.call(rbind, summ)
write.csv(SUMM, file.path(TAB, "A2e_留一法重排汇总.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("\n========== A2e 留一法重排汇总（每层用该层自己的主估计量）==========\n")
cat(sprintf("各排名集规模：%s\n",
            paste(sprintf("%s=%d", names(sz), sz), collapse = "、")))
cat("（每层只把该层主估计量放进名次集；复制层的 IVW-FE 行仅供对照，不参与该层主估计的名次。）\n")
print(SUMM, row.names = FALSE)

cat("\n---------- 面板与交联家族的名次变动（|变动| ≥ 3 者）----------\n")
big <- cmp[(cmp$in_panel | cmp$is_crosslink) & (abs(cmp$d_abs) >= 3 | abs(cmp$d_P) >= 3),
           c("layer", "estimator", "gene", "dropped", "rank_abs_full", "rank_abs_loo", "d_abs",
             "rank_P_full", "rank_P_loo", "d_P")]
if (nrow(big) == 0) cat("（无：受关注基因在任一剔除下的名次变动均 < 3）\n") else print(big, row.names = FALSE)

cat(sprintf("\n已写出：\n  %s  %d 行\n  %s  %d 行\n",
            file.path(TAB, "A2e_留一法重排汇总.csv"), nrow(SUMM),
            file.path(TAB, "A2e_留一法_逐基因名次变动.csv"), nrow(cmp)))
