# B8 候选范围与选择依据：57 基因总览（§4.1）
#
# 目的：把 6 基因面板放回它被挑出来的那 57 个基因里，让读者能看见面板的位置，
#   而不是只看见被挑中的那 6 个。这是「候选选择透明」那一块证据的实体。
#
# 纪律：
#  ① 只读 B1/B2 的落盘产物，不重算任何效应。
#  ② 一致性计数必须带方向。「3/3 一致」这个说法本身不含方向，据此标上下会把
#     三个下调基因标成「3/3↑」（见笔记 agreement-count-column-carries-no-direction）。
#     故本表同时给 n_agree 与 direction 两列，且方向取自逐队列值的符号而不另行判断。
#  ③ 一致性分子分母都要报：分母是该基因**实际有测量值**的队列数，不是 6。
#     分母为 0 或 1 时不报比例——一个队列谈不上「一致」。
#
# 输入：02_结果表/逐队列效应_力学基因.csv（B1）、合并效应_{发现,复制}层.csv（B2）
# 输出：02_结果表/57基因总览.csv、02_结果表/B8_57基因总览_结果.md
# 本脚本只用 base R，不加载任何包（原先那句空的 suppressPackageStartupMessages({})
# 会把 NULL 打进日志头一行）。
# 定义必须在最前：后面若干处 `paste(...) %||% "（无）"` 都是空串兜底，
# 放到文件末尾会在用到时还没定义。
`%||%` <- function(a, b) if (length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

BASE <- "/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/01.材料/01.博士材料/10.毕业/01.小论文/04.立项/滑膜基质刚度生成与刚度感知 0912"
# 2026-10-01：落点改由环境变量控制。默认仍是旧树（历史可复现），但本轮修订必须显式传
# B8_RES 指向修订工作树，否则会读旧口径的合并表、并把新口径的 57基因总览写回旧树。
RES  <- Sys.getenv("B8_RES", unset = file.path(BASE, "08_修订_2026-09-17/02_结果表"))

DISCOVERY   <- c("GSE89408", "GSE206848", "GSE1919")
REPLICATION <- c("GSE55235", "GSE55457", "GSE82107")
ALL_COH     <- c(DISCOVERY, REPLICATION)
PANEL       <- c("LOXL1", "VGLL4", "SDC4", "SERPINE1", "ITGA5", "SAV1")
CROSSLINK   <- c("LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4", "PLOD1", "PLOD2", "TGM2")

per <- read.csv(file.path(RES, "逐队列效应_力学基因.csv"), stringsAsFactors = FALSE)
dis <- read.csv(file.path(RES, "合并效应_发现层.csv"), stringsAsFactors = FALSE)
rep <- read.csv(file.path(RES, "合并效应_复制层.csv"), stringsAsFactors = FALSE)
# 2026-10-01：两层主估计统一为 REML+HK。断言字形须与 B2 写出的 primary_method 逐字相同 ——
# 这里对不上就不许出表（否则本表所有 *_logFC/*_P 列会挂着旧估计量的值继续流向下游）。
stopifnot(all(dis$primary_method == "REML+HK"),
          all(rep$primary_method == "REML+HK"))
stopifnot(all(dis$primary_logFC == dis$hk_logFC), all(rep$primary_logFC == rep$hk_logFC))

GENES <- sort(unique(per$gene))
cat(sprintf("基因集：%d 个基因；逐队列记录 %d 行\n", length(GENES), nrow(per)))
stopifnot(length(GENES) == 57L, all(CROSSLINK %in% GENES))

# 逐队列方向计数：分母是**实际有测量值**的队列，分子是符号与「该层主估计」相同的队列数。
# 与主估计比符号、而不与队列间多数比，是为了让 n_agree 与表里的 direction 出自同一处。
dir_count <- function(gene, cohorts, ref) {
  v <- per$logFC[per$gene == gene & per$dataset %in% cohorts]
  if (is.na(ref) || !length(v)) return(c(n = length(v), a = NA_integer_))
  c(n = length(v), a = sum(sign(v) == sign(ref)))
}

rows <- list()
for (g in GENES) {
  i <- match(g, dis$gene); j <- match(g, rep$gene)
  d <- if (is.na(i)) NULL else dis[i, ]
  r <- if (is.na(j)) NULL else rep[j, ]
  dc <- dir_count(g, DISCOVERY, if (is.null(d)) NA_real_ else d$primary_logFC)
  rc <- dir_count(g, REPLICATION, if (is.null(r)) NA_real_ else r$primary_logFC)
  rows[[g]] <- data.frame(
    gene = g,
    in_panel = g %in% PANEL,
    is_crosslink = g %in% CROSSLINK,
    n_coh_total = sum(per$gene == g),
    disc_k = if (is.null(d)) NA_integer_ else d$k,
    disc_logFC = if (is.null(d)) NA_real_ else d$primary_logFC,
    disc_CI_low = if (is.null(d)) NA_real_ else d$primary_CI_low,
    disc_CI_high = if (is.null(d)) NA_real_ else d$primary_CI_high,
    disc_P = if (is.null(d)) NA_real_ else d$primary_P,
    disc_I2 = if (is.null(d)) NA_real_ else d$I2,
    # 被降级的敏感性列（DL 随机效应 Meta）。自 2026-10-01 起，两层的 `primary_*` 主估计
    # 同为 REML+HK（`hk_*`），`re_*` 在两层都是 **IVW-FE 列**的 DL 随机效应版本
    # （见 B2 的 re_layer_note 列）——两层含义相同，不再是旧版的「一层是主估计的 RE 版本」。
    disc_re_logFC = if (is.null(d)) NA_real_ else d$re_logFC,
    disc_re_CI_low = if (is.null(d)) NA_real_ else d$re_CI_low,
    disc_re_CI_high = if (is.null(d)) NA_real_ else d$re_CI_high,
    disc_re_P = if (is.null(d)) NA_real_ else d$re_P,
    disc_agree_n = dc[["n"]], disc_agree_k = dc[["a"]],
    rep_k = if (is.null(r)) NA_integer_ else r$k,
    rep_logFC = if (is.null(r)) NA_real_ else r$primary_logFC,
    rep_CI_low = if (is.null(r)) NA_real_ else r$primary_CI_low,
    rep_CI_high = if (is.null(r)) NA_real_ else r$primary_CI_high,
    rep_P = if (is.null(r)) NA_real_ else r$primary_P,
    rep_I2 = if (is.null(r)) NA_real_ else r$I2,
    rep_re_logFC = if (is.null(r)) NA_real_ else r$re_logFC,
    rep_re_CI_low = if (is.null(r)) NA_real_ else r$re_CI_low,
    rep_re_CI_high = if (is.null(r)) NA_real_ else r$re_CI_high,
    rep_re_P = if (is.null(r)) NA_real_ else r$re_P,
    rep_agree_n = rc[["n"]], rep_agree_k = rc[["a"]],
    stringsAsFactors = FALSE)
}
ov <- do.call(rbind, rows)

# 方向用词直接由估计值的符号生成，不手写 —— 手写的方向词与取到的系数是两个字符串，
# 错配后数值看着仍合理（见笔记 summary-table-label-must-bind-to-taken-coefficient）。
word <- function(x) ifelse(is.na(x), NA_character_, ifelse(x > 0, "上调", "下调"))
ov$disc_dir <- word(ov$disc_logFC)
ov$rep_dir  <- word(ov$rep_logFC)
# 三列都必须区分「相反」与「不可评估」。写成 !is.na(a) & !is.na(b) & a==b 的话，
# 不可评估会落成 FALSE —— 与「方向相反」共用同一个取值，读 CSV 的人会把 5 个
# 复制层只有一个队列、根本算不出合并估计的基因读成「两层方向相反」。
ov$same_dir <- ifelse(is.na(ov$disc_logFC) | is.na(ov$rep_logFC), NA,
                      sign(ov$disc_logFC) == sign(ov$rep_logFC))
ov$both_sig <- ifelse(is.na(ov$disc_P) | is.na(ov$rep_P), NA,
                      ov$disc_P < 0.05 & ov$rep_P < 0.05)
# 判据 B：两层统一改用 DL 随机效应 Meta —— 本轮起这是**被降级的敏感性对照**。
# 2026-10-01 之前，A 的两层是两个不同的模型（发现 IVW-FE ／ 复制 ~0+group+cohort），
# 故 A 与 B 是两条不可混读的独立判据。本轮两层主估计统一为 REML+HK，
# **A 与 B 现在都是「两层同一个估计量」**，二者的差别只剩估计量本身（HK 校正的有无），
# 不再是「口径不一致 vs 一致」。这个区别必须在正文里写准。
ov$both_sig_re <- ifelse(is.na(ov$disc_re_P) | is.na(ov$rep_re_P), NA,
                         ov$disc_re_P < 0.05 & ov$rep_re_P < 0.05)
ov$disc_sig_re <- ifelse(is.na(ov$disc_re_P), NA, ov$disc_re_P < 0.05)
ov$rep_sig_re <- ifelse(is.na(ov$rep_re_P), NA, ov$rep_re_P < 0.05)
# 单层内方向是否全同：只有 ≥2 个队列时才有定义
ov$disc_unanimous <- ifelse(ov$disc_agree_n < 2, NA, ov$disc_agree_k == ov$disc_agree_n)
ov$rep_unanimous  <- ifelse(ov$rep_agree_n  < 2, NA, ov$rep_agree_k  == ov$rep_agree_n)

## 排名：把面板放回 57 个基因里看位置。按 |logFC| 与 P 各排一次，
## 两张名次都落盘，避免只挑好看的那张报。
ov$rank_abs_disc <- rank(-abs(ov$disc_logFC), na.last = "keep", ties.method = "min")
ov$rank_P_disc   <- rank(ov$disc_P,            na.last = "keep", ties.method = "min")
ov$rank_abs_rep  <- rank(-abs(ov$rep_logFC),   na.last = "keep", ties.method = "min")
ov$rank_P_rep    <- rank(ov$rep_P,             na.last = "keep", ties.method = "min")

ov <- ov[order(ov$disc_P, -abs(ov$disc_logFC)), ]
write.csv(ov, file.path(RES, "57基因总览.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- 汇总 ----------
say <- function(...) cat(sprintf(...), "\n")
say("\n================ 57 基因总览（§4.1）================")
say("基因数 %d；其中刚度生成—交联一臂 %d 个（%s）", nrow(ov), sum(ov$is_crosslink),
    paste(ov$gene[ov$is_crosslink], collapse = ", "))
say("\n两层都能合并（k≥2）的基因：%d / %d", sum(!is.na(ov$disc_logFC) & !is.na(ov$rep_logFC)), nrow(ov))
say("  仅发现层有：%s", paste(ov$gene[!is.na(ov$disc_logFC) & is.na(ov$rep_logFC)], collapse = ", "))
say("  仅复制层有：%s", paste(ov$gene[is.na(ov$disc_logFC) & !is.na(ov$rep_logFC)], collapse = ", ") %||% "（无）")

both <- ov[!is.na(ov$disc_logFC) & !is.na(ov$rep_logFC), ]
say("\n可评估两层方向一致性的：%d / %d", nrow(both), nrow(ov))
say("  不可评估的 %d 个：%s（复制层仅 1 个队列，合并估计不存在）",
    nrow(ov) - nrow(both), paste(ov$gene[is.na(ov$rep_logFC)], collapse = ", "))
say("  两层方向一致的：%d / %d", sum(both$same_dir), nrow(both))
say("    其中两层均 P<0.05：%d 个（%s）", sum(both$both_sig),
    paste(both$gene[both$both_sig], collapse = ", "))
say("  两层方向不一致的：%d 个", sum(!both$same_dir))
# 落盘列与打印计数来自同一次判定，必须对得上：两处各自算一遍是这类表最容易
# 「表里一个数、文字里另一个数」的地方。
stopifnot(sum(ov$same_dir, na.rm = TRUE) == sum(both$same_dir),
          sum(!ov$same_dir, na.rm = TRUE) == sum(!both$same_dir),
          sum(ov$same_dir, na.rm = TRUE) + sum(!ov$same_dir, na.rm = TRUE) == nrow(both))
say("    %s",
    paste(sprintf("%s(发现%s %+.2f / 复制%s %+.2f)", both$gene[!both$same_dir],
                  both$disc_dir[!both$same_dir], both$disc_logFC[!both$same_dir],
                  both$rep_dir[!both$same_dir], both$rep_logFC[!both$same_dir]),
          collapse = "; ") %||% "（无）")

say("\n层内逐队列方向全同（≥2 队列）的基因数：发现层 %d、复制层 %d",
    sum(ov$disc_unanimous, na.rm = TRUE), sum(ov$rep_unanimous, na.rm = TRUE))

## ---------- 判据 A 与判据 B：主估计（REML+HK） vs 被降级的 DL ----------
# 这一段是给审稿人看的：发现层 I² 普遍偏高，若换成随机效应，合格基因集会缩到多小。
# 本轮起两层的 A 与 B **都用同一个估计量跨层**，差别只在 HK 校正的有无 ——
# 这正是审稿意见要的那句「主估计是 REML+HK，DL 只作描述性敏感性」。
say("\n---------- 判据 A（主估计：两层统一 REML+HK）vs 判据 B（被降级：两层统一 DL）----------")
say("判据 A：两层统一 REML+HK（主估计，`method=\"REML\", test=\"knha\"`）")
say("  方向一致且两层均 P<0.05：%d 个（%s）", sum(both$both_sig),
    paste(both$gene[both$both_sig], collapse = ", "))
say("判据 B：两层统一 DL 随机效应 Meta（该层主估计的敏感性对照）")
say("  方向一致且两层均 P<0.05：%d 个（%s）", sum(both$both_sig_re),
    paste(both$gene[both$both_sig_re], collapse = ", "))
say("  仅 A 合格而 B 不合格：%s",
    paste(both$gene[both$both_sig & !both$both_sig_re], collapse = ", ") %||% "（无）")
say("  仅 B 合格而 A 不合格：%s",
    paste(both$gene[!both$both_sig & both$both_sig_re], collapse = ", ") %||% "（无）")
say("  方向在 A 与 B 之间改变的基因：%s",
    paste(both$gene[sign(both$disc_logFC) != sign(both$disc_re_logFC) |
                    sign(both$rep_logFC) != sign(both$rep_re_logFC)], collapse = ", ") %||% "（无）")
say("\n面板 6 基因：主估计（REML+HK）P 与 DL 对照 P（1 = 该层 P<0.05）：")
for (g in PANEL) {
  r <- ov[ov$gene == g, ]
  say("  %-9s 发现 %s（P=%.3g→%.3g）  复制 %s（P=%.3g→%.3g）",
      g, ifelse(isTRUE(r$disc_sig_re), "合格", "不合格"), r$disc_P, r$disc_re_P,
      ifelse(isTRUE(r$rep_sig_re), "合格", "不合格"), r$rep_P, r$rep_re_P)
}

say("\n---------- 面板 6 基因在 57 个基因里的位置 ----------")
say("%-9s %9s %9s %9s %9s  %s", "基因", "|logFC|名次", "P名次", "复制|logFC|名次", "复制P名次", "两层")
for (g in PANEL) {
  r <- ov[ov$gene == g, ]
  say("%-9s %6d/%-2d %6d/%-2d %6d/%-2d %6d/%-2d  %s",
      g, r$rank_abs_disc, nrow(ov), r$rank_P_disc, nrow(ov),
      r$rank_abs_rep, nrow(ov), r$rank_P_rep, nrow(ov),
      ifelse(r$same_dir, "方向一致", "方向相反"))
}

say("\n---------- 全 57 基因里排在最前的 10 个（按发现层 P）----------")
say("%-10s %8s %20s %10s %8s  %s", "基因", "log2FC", "95%CI(发现)", "P", "I2%", "面板")
for (i in 1:10) {
  r <- ov[i, ]
  say("%-10s %+8.3f [%+.3f,%+.3f] %10.3g %7.1f  %s", r$gene, r$disc_logFC,
      r$disc_CI_low, r$disc_CI_high, r$disc_P, r$disc_I2, ifelse(r$in_panel, "★", ""))
}

## 面板 Bonferroni：分母按**实际检验过的基因数**（7，含后来撤出报告的 CCN1/CYR61），
## 不按最终报告的 6 —— 撤出基因不减少已经做过的检验（见笔记
## panel-bonferroni-denominator-is-tests-actually-run）。
ALPHA_PANEL <- 0.05 / 7
say("\n面板 Bonferroni 阈值 α=0.05/7=%.5f（分母 7 = 实际检验过的基因数）", ALPHA_PANEL)
# 用 isTRUE() 兜底：该基因对该层没有合并行时 ov[[col]] 取到 NA，裸比较会得到 NA，
# 而 NA 当下标会往结果里插进一个 NA 基因名。
hit <- function(g, col) isTRUE(ov[[col]][ov$gene == g] < ALPHA_PANEL)
say("  发现层达阈值的面板基因：%s",
    paste(PANEL[sapply(PANEL, hit, col = "disc_P")], collapse = ", ") %||% "（无）")
say("  复制层达阈值的面板基因：%s",
    paste(PANEL[sapply(PANEL, hit, col = "rep_P")], collapse = ", ") %||% "（无）")

say("\n已写出 02_结果表/57基因总览.csv（%d 行 × %d 列）", nrow(ov), ncol(ov))
