# A-2b. 主估计（REML+HK）下两层各自的名义显著基因清单
#
# 为什么单独立一个脚本：主估计换成 REML+HK 后，正文的中心断言从
#   「PLOD1 是两层唯一显著者」变成「**没有任何候选在两层同时显著**」。
# 这是一个**全基因集**（57 / 52）上的计数断言，而 A2 只在交联家族 8 基因内做跨线计数
# （其 n_eval ≤ 8，是过滤后的明细表，不能当普查分母用）。故此处独立成表，
# 分母就是各层的可汇总基因数本身，并把「面板 6 基因」与「交联家族 8 基因」两个
# 受关注的子集单独标出，避免读者从全表里手挑。
#
# 纪律：本脚本**不做任何新估计**，只读 B2 写出的冻结合并表并计数。
#
# 输入：$A2B_RES/合并效应_发现层.csv、合并效应_复制层.csv（B2 产出）
# 输出：$A2B_RES/A2b_主估计显著基因.csv        逐基因一行（两层合并视图）
#       $A2B_AUD/A2b_主估计显著_结果.md        计数与点名
suppressPackageStartupMessages({ library(stats) })

BASE <- paste0("/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/",
               "Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/",
               "01.材料/01.博士材料/10.毕业/01.小论文/04.立项/",
               "滑膜基质刚度生成与刚度感知 0912/08_修订_2026-09-17")
RES <- Sys.getenv("A2B_RES", unset = file.path(BASE, "02_结果表"))
AUD <- Sys.getenv("A2B_AUD", unset = file.path(BASE, "01_审计"))
dir.create(AUD, showWarnings = FALSE, recursive = TRUE)

ALPHA     <- 0.05
PANEL     <- c("LOXL1", "VGLL4", "SDC4", "SERPINE1", "ITGA5", "SAV1")
CROSSLINK <- c("LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4", "PLOD1", "PLOD2", "TGM2")
PRIMARY   <- "REML+HK"

dis <- read.csv(file.path(RES, "合并效应_发现层.csv"), stringsAsFactors = FALSE)
rep <- read.csv(file.path(RES, "合并效应_复制层.csv"), stringsAsFactors = FALSE)

# 口径断言：主估计必须是 REML+HK，且 primary_* 确实取自 hk_*（不许挂了名字却取旧列）
stopifnot(all(dis$primary_method == PRIMARY), all(rep$primary_method == PRIMARY))
stopifnot(all(dis$primary_logFC == dis$hk_logFC), all(rep$primary_logFC == rep$hk_logFC))
stopifnot(all(dis$primary_P    == dis$hk_P),    all(rep$primary_P    == rep$hk_P))
# 四个估计量必须都在表里：位移归因要用到 rml_*（REML 但不含 HK）
stopifnot(all(c("ivw_P", "re_P", "rml_P", "hk_P") %in% names(dis)),
          all(c("ivw_P", "re_P", "rml_P", "hk_P") %in% names(rep)))

say <- function(...) cat(sprintf(...), "\n")
ln  <- character(0)
rec <- function(s = "") { cat(s, "\n"); ln <<- c(ln, s) }

# 逐层：可汇总 = primary_P 非 NA 且有限（k≥2 才合并得出来）
eval_of <- function(tab) tab[is.finite(tab$primary_P), ]

sig_names <- function(tab) {
  x <- eval_of(tab)
  x$gene[x$primary_P < ALPHA]
}
sig_named <- function(tab) {
  x <- eval_of(tab)
  x <- x[order(x$primary_P), ]
  x$gene[x$primary_P < ALPHA]
}

d_sig <- sig_named(dis); r_sig <- sig_named(rep)
d_eval <- eval_of(dis);   r_eval <- eval_of(rep)

rec("# A-2b 主估计（REML+HK）下的名义显著清单")
rec()
rec("> 由 `03_分析/scripts/A2b_主估计显著基因清单.R` 生成。**本脚本不做新估计**，")
rec("> 分母是各层**可汇总**的基因数（`primary_P` 有限者），不是 57 / 52 的名义候选数。")
rec()

## ---------- 一、四估计量下的显著数（位移归因） ----------
n_of <- function(tab, col) sum(eval_of(tab)[[col]] < ALPHA)
tab4 <- data.frame(
  layer = c("发现", "复制"),
  n_evaluable = c(nrow(d_eval), nrow(r_eval)),
  FE   = c(n_of(dis, "ivw_P"), n_of(rep, "ivw_P")),
  DL   = c(n_of(dis, "re_P"),  n_of(rep, "re_P")),
  REML = c(n_of(dis, "rml_P"), n_of(rep, "rml_P")),
  REML_HK = c(length(d_sig),   length(r_sig)),
  stringsAsFactors = FALSE)
rec("## 一、四个估计量各自的名义显著数（P<0.05）")
rec()
rec("| 层 | 可汇总 n | FE | DL | REML | REML+HK（主） |")
rec("|---|---|---|---|---|---|")
for (i in seq_len(nrow(tab4))) {
  rec(sprintf("| %s | %d | %d | %d | %d | **%d** |", tab4$layer[i], tab4$n_evaluable[i],
              tab4$FE[i], tab4$DL[i], tab4$REML[i], tab4$REML_HK[i]))
}
rec()
# REML 与 DL 逐格同数 ⇒ 位移全部来自 HK。这条不是修辞，是可直接核的等式。
same_rd <- all(tab4$REML == tab4$DL)
rec(sprintf("REML（不含 HK）与 DL 的显著数逐层%s；故位移%s来自 HK 那一层校正。",
            ifelse(same_rd, "**完全相同**", "**不同**"),
            ifelse(same_rd, "全部", "并非全部")))
rec()

## ---------- 二、两层显著清单与交集 ----------
rec("## 二、两层各自的名义显著基因，及其交集")
rec()
rec(sprintf("发现层（可汇总 %d 个）显著 %d 个：%s", nrow(d_eval), length(d_sig),
            ifelse(length(d_sig) == 0, "（无）", paste(d_sig, collapse = "、"))))
rec(sprintf("复制层（可汇总 %d 个）显著 %d 个：%s", nrow(r_eval), length(r_sig),
            ifelse(length(r_sig) == 0, "（无）", paste(r_sig, collapse = "、"))))
inter <- intersect(d_sig, r_sig)
rec()
rec(sprintf("**两层交集：%s**。", ifelse(length(inter) == 0, "空集 —— 无任何候选在两层同时名义显著",
                                          paste(inter, collapse = "、"))))
rec()
rec("两个受关注子集在这两层清单里的位置：")
for (nm in c("面板 6 基因", "交联家族 8 基因")) {
  S <- if (nm == "面板 6 基因") PANEL else CROSSLINK
  rec(sprintf("  · %s：发现层显著 %s；复制层显著 %s。", nm,
              ifelse(length(intersect(S, d_sig)) == 0, "（无）", paste(intersect(S, d_sig), collapse = "、")),
              ifelse(length(intersect(S, r_sig)) == 0, "（无）", paste(intersect(S, r_sig), collapse = "、"))))
}
rec()
# 「只有一层的」——这是「无跨层复制」这句话的实际内容，必须点名，不能只报交集为空
rec("**只在单层通过的候选**（交集为空时，这批就是全部的信号所在）：")
rec(sprintf("  · 仅发现层：%s", ifelse(length(setdiff(d_sig, r_sig)) == 0, "（无）",
                                        paste(setdiff(d_sig, r_sig), collapse = "、"))))
rec(sprintf("  · 仅复制层：%s", ifelse(length(setdiff(r_sig, d_sig)) == 0, "（无）",
                                        paste(setdiff(r_sig, d_sig), collapse = "、"))))
rec()

## ---------- 三、逐基因表（两层合并视图） ----------
join_by <- function(col_d, col_r) {
  v <- rep[[col_r]][match(dis$gene, rep$gene)]
  list(d = dis[[col_d]], r = v)
}
out <- data.frame(
  gene = dis$gene,
  in_panel     = dis$gene %in% PANEL,
  is_crosslink = dis$gene %in% CROSSLINK,
  disc_k = dis$k,
  disc_primary_logFC = dis$primary_logFC,
  disc_primary_CI_low = dis$primary_CI_low,
  disc_primary_CI_high = dis$primary_CI_high,
  disc_primary_P = dis$primary_P,
  disc_FE_P = dis$ivw_P, disc_DL_P = dis$re_P, disc_REML_P = dis$rml_P,
  disc_sig_primary = is.finite(dis$primary_P) & dis$primary_P < ALPHA,
  rep_k = rep$k[match(dis$gene, rep$gene)],
  rep_primary_logFC = rep$primary_logFC[match(dis$gene, rep$gene)],
  rep_primary_CI_low = rep$primary_CI_low[match(dis$gene, rep$gene)],
  rep_primary_CI_high = rep$primary_CI_high[match(dis$gene, rep$gene)],
  rep_primary_P = rep$primary_P[match(dis$gene, rep$gene)],
  rep_FE_P = rep$ivw_P[match(dis$gene, rep$gene)],
  rep_DL_P = rep$re_P[match(dis$gene, rep$gene)],
  rep_REML_P = rep$rml_P[match(dis$gene, rep$gene)],
  rep_sig_primary = is.finite(rep$primary_P[match(dis$gene, rep$gene)]) &
                    rep$primary_P[match(dis$gene, rep$gene)] < ALPHA,
  # 两层皆可汇总才谈得上「跨层」；只在一层可汇总的基因该列取 NA，不取 FALSE。
  # 取 FALSE 会把「不可评估」读成「评估了且未通过」（记忆 agreement-count-column-carries-no-direction）。
  sig_both_layers = ifelse(is.finite(dis$primary_P) & is.finite(rep$primary_P[match(dis$gene, rep$gene)]),
                           is.finite(dis$primary_P) & dis$primary_P < ALPHA &
                           rep$primary_P[match(dis$gene, rep$gene)] < ALPHA, NA),
  stringsAsFactors = FALSE)
stopifnot(nrow(out) == nrow(dis))
# 落盘列与上面点名的清单必须出自同一判定（不许各算一遍）
stopifnot(setequal(out$gene[!is.na(out$disc_sig_primary) & out$disc_sig_primary], d_sig),
          setequal(out$gene[!is.na(out$rep_sig_primary)  & out$rep_sig_primary],  r_sig))
write.csv(out, file.path(RES, "A2b_主估计显著基因.csv"), row.names = FALSE, fileEncoding = "UTF-8")

rec("## 三、逐基因表（两层合并视图）")
rec()
rec("| 基因 | 面板 | 家族 | 发现 k | 发现 logFC | 发现 95%CI | 发现 主 P | 发现 FE P | 发现 DL P | 复制 k | 复制 logFC | 复制 95%CI | 复制 主 P | 复制 FE P | 复制 DL P | 两层皆显著 |")
rec("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
fm  <- function(x) ifelse(is.na(x), "—", formatC(x, format = "f", digits = 3, flag = "+"))
fP  <- function(x) ifelse(is.na(x), "—", formatC(x, format = "g", digits = 3))
fCI <- function(lo, hi) ifelse(is.na(lo) | is.na(hi), "—", sprintf("(%s, %s)", fm(lo), fm(hi)))
fB  <- function(x) ifelse(is.na(x), "不可评估", ifelse(x, "**是**", "否"))
for (i in seq_len(nrow(out))) {
  z <- out[i, ]
  rec(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
              z$gene, ifelse(z$in_panel, "是", ""), ifelse(z$is_crosslink, "是", ""),
              ifelse(is.na(z$disc_k), "—", z$disc_k),
              fm(z$disc_primary_logFC), fCI(z$disc_primary_CI_low, z$disc_primary_CI_high),
              fP(z$disc_primary_P), fP(z$disc_FE_P), fP(z$disc_DL_P),
              ifelse(is.na(z$rep_k), "—", z$rep_k),
              fm(z$rep_primary_logFC), fCI(z$rep_primary_CI_low, z$rep_primary_CI_high),
              fP(z$rep_primary_P), fP(z$rep_FE_P), fP(z$rep_DL_P),
              fB(z$sig_both_layers)))
}
rec()
n_both <- sum(out$sig_both_layers, na.rm = TRUE)
rec(sprintf("两层皆显著（且两层皆可汇总）的基因数：**%d**。", n_both))
rec()

writeLines(c(ln, "", "---",
             sprintf("发现层 %d/%d 显著，复制层 %d/%d 显著，交集 %d 个。",
                     length(d_sig), nrow(d_eval), length(r_sig), nrow(r_eval), length(inter)),
             "**本脚本不做新估计；全部数值来自 B2 的冻结合并表，可逐格回查。**"),
           file.path(AUD, "A2b_主估计显著_结果.md"))

cat("\n已写出：\n")
cat("  ", file.path(RES, "A2b_主估计显著基因.csv"), "\n")
cat("  ", file.path(AUD, "A2b_主估计显著_结果.md"), "\n")
