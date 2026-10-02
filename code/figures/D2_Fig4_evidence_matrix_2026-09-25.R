# D2 — 图 4：八交联候选的证据矩阵
#
# 稿件 Figure 4 图注规定的列序（逐列照做，不自作增删）：**17 列** ——
#   discovery层主估计 REML+HK 的 log2FC / P / I²（3 列）、该层敏感性 IVW-FE 的 log2FC 与 P（2 列）、
#   replication层主估计 REML+HK 的 log2FC / P / I²（3 列）、该层敏感性 pooled 的 log2FC 与 P（2 列）、
#   两层各自 |log2FC| 与 P 判据下的名次（4 列）、三次单队列删除的结果（3 列）。
#
# 三处口径必须先定准，否则图会与正文打架（均已回源核实）：
#
#  ① **名次列取「主估计量下的名次」**。图注写的是 "each layer's |log2FC| and P criteria"。
#     ★ 2026-10-01：§2.3 已改声明**两层主估计 = REML+HK**（IVW-FE 与 pooled 降为敏感性）。
#     B8 的 rank_*_disc / rank_*_rep 取自 candidate_overview_57.csv 的 disc_logFC/disc_P 与
#     rep_logFC/rep_P，而 B8 已改为从 primary_*（即 hk_*）取列，故这两列现在**就是**
#     主估计下的名次；下面再加一条断言把它与 S4a 的 RE_rank 逐格钉死（S4a 的
#     re_estimator 列已标注 REML+HK）。列组名自陈估计量：`Rank (REML+HK, both layers)`。
#     PLOD1 在replication层是 22nd of 52 与 2nd of 52（图注与正文 §3.3 所引的主估计名次），
#     在 pooled 下 23rd/8th、在 IVW-FE 下 21st/15th —— 三者都对，只有主估计那一对进本图。
#
#  ② **删除列取discovery层、且用discovery层 IVW-FE**。★ 2026-10-01：IVW-FE 已不是该层主估计，
#     但它仍是本图这一列所用的估计量，列组名照旧自陈为 `Deletion (discovery, IVW-FE)`。
#     旧注释援引的 §3.2 句（「the layer's own primary estimator, IVW-FE…」）已被 F21 改写，
#     此处改为按估计量陈述，不再援引那句。图注说
#     「the three single-cohort deletions」只有三列，而discovery层恰有三个队列；replication层的三次
#     删除（§3.2 已述：无gene变号）留在 Table S4。本图组标签写明层与估计量，避免歧义。
#
#  ③ **LOXL4 / LOXL3 的删除列为空**：二者discovery层只可评估 2 个队列，删一剩 k=1，无法合并，
#     S4c 中本就没有它们的删除行 —— 记作 n/a，而不是 0 或空白（空白会被读成"未做"）。
#
# 着色：逐格按**该列自己的轴**连续着色（图注原文 "shaded continuously within each column"）。
#   log2FC 四列共用发散色标（同一量，共用才可跨层跨估计量读）；P 列（含删除列）按 −log10 P；
#   I² 列按 0–100%；名次列按列内百分位（1 名最深）。深浅只表示"在该列轴上的高低"，
#   不构成任何合成分 —— 与图注末句一致。每格印出原值，色仅作辅助。
#
# 输入：02_结果表/crosslink_family_criterion_sensitivity.csv、candidate_overview_57.csv、leave_one_out_rank_shift.csv
# 输出：04_图表/Fig_evidence_matrix.pdf / .png；02_结果表/Fig4_values.csv（审计用）
suppressPackageStartupMessages({library(ggplot2); library(ggnewscale)})

# ★ 2026-09-29 路径修正（与 B7/D3 同款）。原 BASE 写死上一棵树
#   （…/滑膜基质刚度生成与刚度感知 0912/08_修订_2026-09-17）：本脚本读的是旧树的输入、
#   图也落进旧树，本树 04_图表/Fig_evidence_matrix.* 与交付目录都不会被更新，而退出码是 0。
#   实测形态：本树的 Fig_evidence_matrix.pdf 与交付的 Figure4 都还是 09-25 的字节，
#   且 §Table 3 对账读的是旧树的 09-23 手写底稿 —— 对账对的是一个早已被 F 链改过的稿子。
#   现改为从本脚本位置推树根（scripts/ → 01_冻结结果与代码），并同时写进交付件位。
# ── 交付布局下的路径（导出时改写；原脚本见分析脚本目录）──────────────────────
# 交付包把绘图代码与源数据并列：本脚本在 绘图代码_EN/，源数据在 源数据_EN/，
# 输出写到同级的 重跑产出/ —— 不覆盖交付的源数据，读者可据此逐件对账。
.here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
.pick <- function(...) { for (p in c(...)) if (length(p) == 1L && !is.na(p) && file.exists(p)) return(normalizePath(p)); NA_character_ }
.out  <- file.path(.here, "重跑产出"); dir.create(.out, showWarnings = FALSE)
# 交付的 源数据_EN/ 是待核对的成品。本脚本**既从它读，也往它写**自己的审计表（原脚本里
# 读入与审计表落盘共用同一个 RES 根）。故先把整目录replication进 重跑产出/，之后一律只碰副本：
# 跑完交付的源数据一个字节都不动，把 重跑产出/源数据_EN/ 与它逐件对账即可。
# （旧写法把 RES 直接指向交付的源数据，读者一跑就把交付件就地覆盖 —— 若他的 R 版本
#   让 write.csv 多写一位小数，交付的源数据就被静默改掉了，而脚本照样退出码 0。）
.src0 <- .pick(file.path(.here, "源数据_EN"), file.path(.here, "..", "源数据_EN"))
if (is.na(.src0)) stop("找不到源数据目录 源数据_EN/（应与本脚本并列或位于其上一级）")
.src  <- file.path(.out, "源数据_EN")
if (!dir.exists(.src)) file.copy(.src0, .out, recursive = TRUE)
if (!dir.exists(.src)) stop("replication源数据到 重跑产出/ 失败：", .src0)
# 图注/表注的来源：**优先现役主稿**（D2 要取 Table 3 的表体逐格对账，图注文件里没有表体；
# E2 要取 Figure 1 图注，两者都有）。主稿名带日期标签，词典序即时间序，故取末一个。
.ms   <- .pick(tail(sort(list.files(file.path(.here, "..", "..", "02_共用主稿"),
                                    "^English_manuscript_.*[.]md$", full.names = TRUE)), 1L),
               file.path(.here, "Figure_Legends.md"), file.path(.here, "..", "Figure_Legends.md"))
# D2 的 Table 3 表体只在主稿里，图注文件没有它。只在 03_图表与源数据/ 一级被单独拷走时
# 才会走到这里 —— 要么把主稿一起带上，要么明确报错：静默跑出来的图 4 会缺掉对账那一半。
if (is.na(.ms)) stop("取不到主稿或图注文件：D2 需从主稿取 Table 3 表体、E2 需取 Figure 1 图注。",
                     "\n  主稿应在 ../../02_共用主稿/English_manuscript_*.md；本脚本在 ", .here)
RES <- .src
OUT <- .out
DELIV <- .out
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

GENES <- c("LOXL1", "LOX", "PLOD1", "PLOD2", "LOXL2", "TGM2", "LOXL4", "LOXL3")  # 同 Table 3 序
PANEL <- "LOXL1"                      # 八者中唯一进入报告面板者
N_DISC <- 57L; N_REP <- 52L           # 各层名次的分母（S4a 的 n_in_set，已核实）

fam  <- read.csv(file.path(RES, "crosslink_family_criterion_sensitivity.csv"), stringsAsFactors = FALSE)
ovw  <- read.csv(file.path(RES, "candidate_overview_57.csv"),            stringsAsFactors = FALSE)
loo  <- read.csv(file.path(RES, "leave_one_out_rank_shift.csv"), stringsAsFactors = FALSE)
fam$gene <- as.character(fam$gene); ovw$gene <- as.character(ovw$gene); loo$gene <- as.character(loo$gene)
fam <- fam[match(GENES, fam$gene), ]; ovw <- ovw[match(GENES, ovw$gene), ]
stopifnot(identical(fam$gene, GENES), identical(ovw$gene, GENES))

# ── 口径断言：名次列的来源必须真是各层主估计，且与图注点名的两处一致
n2 <- function(x) suppressWarnings(as.numeric(x))
p1 <- ovw[ovw$gene == "PLOD1", ]
stopifnot(n2(p1$rank_abs_rep) == 22, n2(p1$rank_P_rep) == 2)   # 图注 "22nd … 2nd"，主估计下
# 两层都取主估计 REML+HK，故 57gene总览 的名次应与 S4a 的 RE_rank **逐格相同**。
# S4a 的 re_estimator 列就是这件事的断言（不是推断）：它若退回 DL/FE，这里必须失败。
s4a <- read.csv(file.path(RES, "rank_shift_by_model.csv"), stringsAsFactors = FALSE)
stopifnot(all(s4a$re_estimator == "REML+HK"))
for (LG in c("discovery", "replication")) {
  s <- s4a[s4a$layer == LG & s4a$set == "all genes", ]
  stopifnot(all(n2(s$n_in_set) == if (LG == "discovery") N_DISC else N_REP))
  for (g in GENES) for (bs in c(abs_logFC = "abs_logFC", P = "P")) {
    r <- s[as.character(s$gene) == g & s$basis == bs, ]
    if (nrow(r) != 1) next
    want <- if (LG == "discovery") {
      if (bs == "abs_logFC") ovw$rank_abs_disc[ovw$gene == g] else ovw$rank_P_disc[ovw$gene == g]
    } else {
      if (bs == "abs_logFC") ovw$rank_abs_rep[ovw$gene == g] else ovw$rank_P_rep[ovw$gene == g]
    }
    stopifnot(n2(r$RE_rank) == n2(want))
  }
}

# ── 删除列：discovery层 × IVW-FE（该层的敏感性估计），三次删除
loo_p <- loo[loo$layer == "discovery" & loo$estimator == "IVW-FE" & loo$is_crosslink == "TRUE", ]
DROPS <- c("GSE1919", "GSE206848", "GSE89408")
stopifnot(setequal(unique(loo_p$dropped), DROPS))

# ── 列定义（顺序即图注顺序）。列组名一律**自陈估计量**：主估计（两层同名 REML+HK）与
#    敏感性（discovery层 IVW-FE、replication层 pooled）各成一组，读者不必回正文才知道名次按谁算。
#    组名只允许 ASCII（见下方非 ASCII 断言：基础 pdf 设备会静默丢字），故只换行、不用破折号。
#    每个组名内嵌 "\n"：两列宽的组（"sensitivity" 那两个）在 12.4 in 画布上只有约 105 pt，
#    单行会顶出组框；换行后最长一行约 86 pt。
G_PD <- "Discovery\nprimary (REML+HK)";        G_DS <- "Discovery\nsensitivity (IVW-FE)"
G_RP <- "Replication\nprimary (REML+HK)";      G_RS <- "Replication\nsensitivity (pooled)"
G_RK <- "Rank\n(REML+HK, both layers)";        G_DL <- "Deletion\n(discovery, IVW-FE)"
COLS <- list(
  list(grp = G_PD, lab = "log2FC", src = "disc_primary_logFC", kind = "lfc"),
  list(grp = G_PD, lab = "P",      src = "disc_primary_P",     kind = "p"),
  list(grp = G_PD, lab = "I2 (%)", src = "disc_I2",            kind = "i2"),
  list(grp = G_DS, lab = "log2FC", src = "disc_ivw_logFC",     kind = "lfc"),
  list(grp = G_DS, lab = "P",      src = "disc_ivw_P",         kind = "p"),
  list(grp = G_RP, lab = "log2FC", src = "rep_primary_logFC",  kind = "lfc"),
  list(grp = G_RP, lab = "P",      src = "rep_primary_P",      kind = "p"),
  list(grp = G_RP, lab = "I2 (%)", src = "rep_I2",             kind = "i2"),
  list(grp = G_RS, lab = "log2FC", src = "rep_pooled_logFC",   kind = "lfc"),
  list(grp = G_RS, lab = "P",      src = "rep_pooled_P",       kind = "p"),
  list(grp = G_RK, lab = "|log2FC| disc", src = "rank_abs_disc", kind = "rank", n = N_DISC, tab = "ovw"),
  list(grp = G_RK, lab = "P disc",        src = "rank_P_disc",   kind = "rank", n = N_DISC, tab = "ovw"),
  list(grp = G_RK, lab = "|log2FC| rep",  src = "rank_abs_rep",  kind = "rank", n = N_REP,  tab = "ovw"),
  list(grp = G_RK, lab = "P rep",         src = "rank_P_rep",    kind = "rank", n = N_REP,  tab = "ovw"),
  list(grp = G_DL, lab = "GSE1919",   kind = "del", drop = "GSE1919"),
  list(grp = G_DL, lab = "GSE206848", kind = "del", drop = "GSE206848"),
  list(grp = G_DL, lab = "GSE89408",  kind = "del", drop = "GSE89408")
)

# 三处格式必须与稿件 Table 3 逐字一致 —— 读者会拿图对着表看，位数不同就是"两个值"。
# Table 3 的舍入约定（不是猜的：下面本图与 Table 3 重叠的 4 个 log2FC 列、4 个 P 列、
# 2 个 I2 列共 80 格，都要与从稿件里**机器解析出来的** Table 3 单元格逐字相等，猜错就报错）：
#   log2FC  : 3 位小数、带正负号
#   P       : >=0.01 取 3 位有效数字（0.0601 / 0.0140 / 0.790）；<0.01 取 2 位
#             有效数字（0.0034 / 0.0028）；<0.001 写科学计数法（3.9e-04）
#   I2      : 1 位小数，恰为 0 时写 "0"
fm <- function(x, d = 3) if (!is.finite(x)) "n/a" else sprintf(paste0("%+.", d, "f"), x)
fp <- function(x) {
  if (!is.finite(x)) return("n/a")
  if (x < 0.001) return(formatC(x, format = "e", digits = 1))
  d <- if (x >= 0.01) 3 else 2
  formatC(signif(x, d), format = "f", digits = max(0L, d - 1L - floor(log10(x))))
}
fi2   <- function(x) if (!is.finite(x)) "n/a" else if (round(x, 1) == 0) "0" else formatC(x, format = "f", digits = 1)
frank <- function(x) if (!is.finite(x)) "n/a" else sprintf("%d", round(x))

rows <- list()
for (ci in seq_along(COLS)) {
  C <- COLS[[ci]]
  for (gi in seq_along(GENES)) {
    g <- GENES[gi]
    if (C$kind == "del") {
      r <- loo_p[loo_p$gene == g & loo_p$dropped == C$drop, ]
      if (nrow(r) == 0) { txt <- "n/a"; fill <- NA_real_; note <- "" }
      else {
        stopifnot(nrow(r) == 1)
        v <- n2(r$logFC); pv <- n2(r$P)
        # 变号标记的参照量取**同一估计量**的全样本值：删除列用的就是discovery层 IVW-FE，
        # 故参照也是 disc_ivw_logFC（不是主估计的 disc_primary_logFC，换量比较会误标 *）。
        full <- n2(fam$disc_ivw_logFC[fam$gene == g])
        note <- if (is.finite(v) && is.finite(full) && sign(v) != sign(full)) " *" else ""
        txt <- paste0(fm(v), note, "\n", fp(pv))
        fill <- -log10(pv)   # 与 P 列同轴；用原始 P 会与 P 列的深浅含义相反
      }
      rows[[length(rows) + 1]] <- data.frame(gene = g, ci = ci, txt = txt, fill = fill,
                                             kind = "del", stringsAsFactors = FALSE)
      next
    }
    tb <- if (identical(C$tab, "ovw")) ovw else fam
    v <- n2(tb[[C$src]][tb$gene == g])
    if (C$kind == "lfc")      { txt <- fm(v); fill <- v }
    else if (C$kind == "p")   { txt <- fp(v); fill <- -log10(v) }
    else if (C$kind == "i2")  { txt <- fi2(v); fill <- v }
    else if (C$kind == "rank"){ txt <- frank(v); fill <- if (is.finite(v)) 1 - (v - 1) / (C$n - 1) else NA_real_ }
    # rank 的 fill 已是 0-1 百分位（按各层名次分母 57/52），下面不再列内归一
    rows[[length(rows) + 1]] <- data.frame(gene = g, ci = ci, txt = txt, fill = fill,
                                           kind = C$kind, stringsAsFactors = FALSE)
  }
}
D <- do.call(rbind, rows)
D$col_lab <- vapply(D$ci, function(i) COLS[[i]]$lab, "")
D$grp     <- vapply(D$ci, function(i) COLS[[i]]$grp, "")
# x 轴键用列号：列标签本身有重名（log2FC 与 P 各出现四次、I2 (%) 两次，分属两层与两组估计量），
# 直接拿标签当 factor 水平会因重复水平报错；键与标签分开，标签由 scale_x_discrete 映射。
D$col_key <- factor(sprintf("%02d", D$ci), levels = sprintf("%02d", seq_along(COLS)))
COL_LAB <- setNames(vapply(COLS, function(C) C$lab, ""), sprintf("%02d", seq_along(COLS)))
D$grp     <- factor(D$grp,     levels = unique(vapply(COLS, function(C) C$grp, "")))
D$gene    <- factor(D$gene, levels = rev(GENES))   # rev：第一水平画在最下
D$ylab    <- factor(ifelse(as.character(D$gene) == PANEL,
                           paste0(as.character(D$gene), " *"), as.character(D$gene)),
                    levels = rev(ifelse(GENES == PANEL, paste0(GENES, " *"), GENES)))

# ── 逐格回验：图上印出的每一格数值必须等于源表（防止 cbind/match 错位后"看着合理但取错列"）
# 删除列的目标串按同一套格式化函数重建，不靠解析已生成的 txt —— 自己解析自己等于没验。
n_checked <- 0L
for (g in GENES) for (ci in seq_along(COLS)) {
  C <- COLS[[ci]]; got <- D$txt[D$gene == g & D$ci == ci]
  if (C$kind == "del") {
    r <- loo_p[loo_p$gene == g & loo_p$dropped == C$drop, ]
    if (nrow(r) == 0) { want <- "n/a" } else {
      stopifnot(nrow(r) == 1)
      v <- n2(r$logFC)
      mark <- if (sign(v) != sign(n2(fam$disc_ivw_logFC[fam$gene == g]))) " *" else ""
      want <- paste0(fm(v), mark, "\n", fp(n2(r$P)))
    }
  } else {
    tb <- if (identical(C$tab, "ovw")) ovw else fam
    v <- n2(tb[[C$src]][tb$gene == g])
    want <- switch(C$kind, lfc = fm(v), p = fp(v), i2 = fi2(v), rank = frank(v))
  }
  if (got != want) stop("回验不符：", g, " ", C$lab, " 图=", got, " 源=", want)
  n_checked <- n_checked + 1L
}
stopifnot(n_checked == length(GENES) * length(COLS))

# ── 与稿件 Table 3 逐格对账。**Table 3 是机器从稿件里解析出来的，不是手抄常数** ——
#    手抄的数一旦与稿子脱钩就查不出来了（本项目有过手写 n 与脚本脱节的先例）。
#    本图有 10 列与 Table 3 重叠（另 4 列名次、3 列删除为 Table 3 所无），这 80 格必须逐字相同。
#    两表的重叠部分**相对列序完全一致**（discovery主 → discovery敏 → replication主 → replication敏），故按名映射即可；
#    Table 3 多出的是 Gene / In panel / 两列 cohorts measured 共 4 列。
MS_STEM <- "English_manuscript_共用清洁稿_2026-09-28i"
MS <- .ms
ms <- readLines(MS, encoding = "UTF-8", warn = FALSE)
i3 <- grep("^\\*\\*Table 3\\.", ms)
stopifnot(length(i3) == 1L)
tbl <- ms[i3:(i3 + 12)]; tbl <- tbl[grepl("^\\|", tbl)]
stopifnot(length(tbl) == 10L)                       # 表头 + 分隔 + 8 数据行
# 稿件表格行写作 `| a | b | ... | c | `（行尾仍有管道符），故 split 后首段为空、
# **末段是真单元格**而不是空串——按 [2:(n-1)] 取会静默丢掉最后一列（Replication I2）。
rowcells <- function(l) {
  x <- strsplit(l, "|", fixed = TRUE)[[1]]
  if (length(x) && !nzchar(trimws(x[1]))) x <- x[-1]
  if (length(x) && !nzchar(trimws(x[length(x)]))) x <- x[-length(x)]
  trimws(x)
}
hdr <- rowcells(tbl[1])
# 表头不手抄就用不了（下面按列名取值），故先把**期望表头**写成常量再与解析结果逐字比：
# 少于或多于一列、或改了任何一列的估计量措辞，这里立刻失败，而不是在 D3 的取值处
# 静默取到 NULL（NULL 会让整列 `got != want` 报长度为零的错，看着像别的问题）。
HDR_WANT <- c("Gene", "In panel",
              "Discovery REML+HK log2FC", "Discovery HK P", "Discovery I² (%)",
              "Discovery IVW-FE log2FC", "Discovery IVW-FE P",
              "Cohorts measured (discovery)", "Cohorts measured (replication)",
              "Replication REML+HK log2FC", "Replication HK P", "Replication I² (%)",
              "Replication pooled log2FC", "Replication pooled P")
stopifnot(length(hdr) == length(HDR_WANT), !anyDuplicated(hdr), identical(hdr, HDR_WANT))
t3 <- as.data.frame(do.call(rbind, lapply(tbl[3:10], rowcells)), stringsAsFactors = FALSE)
names(t3) <- hdr
stopifnot(nrow(t3) == 8L, identical(t3[[1]], GENES))
# Table 3 的科学计数法写作 "3.9×10⁻⁴"（上标），图上只能用 ASCII，故先归一到 "3.9e-04" 再比
SUP <- c("⁰"="0","¹"="1","²"="2","³"="3","⁴"="4","⁵"="5","⁶"="6","⁷"="7","⁸"="8","⁹"="9","⁻"="-")
norm_t3 <- function(s) {
  if (grepl("^—", s)) return("n/a")
  # 稿件里的负号是 U+2212 MINUS SIGN，不是 ASCII 连字符；不归一就会把 -0.226 判成不符
  s <- gsub("−", "-", s, fixed = TRUE)
  if (grepl("×10", s, fixed = TRUE)) {
    mant <- sub("×10.*$", "", s); ex <- sub("^.*×10", "", s)
    for (k in names(SUP)) ex <- gsub(k, SUP[[k]], ex, fixed = TRUE)
    neg <- substr(ex, 1, 1) == "-"; num <- if (neg) substring(ex, 2) else ex
    return(sprintf("%se-%02d", mant, as.integer(num)))
  }
  s
}
D3 <- data.frame(
  ci = 1:10,                        # 本图前 10 列，顺序即图注顺序，与 Table 3 的 10 个数值列一一对应
  col = c("Discovery REML+HK log2FC", "Discovery HK P", "Discovery I² (%)",
          "Discovery IVW-FE log2FC", "Discovery IVW-FE P",
          "Replication REML+HK log2FC", "Replication HK P", "Replication I² (%)",
          "Replication pooled log2FC", "Replication pooled P"),
  stringsAsFactors = FALSE)
stopifnot(all(D3$col %in% hdr))     # 表头变了就在此报错，不在下面的取值处静默取 NULL
n_t3 <- 0L
for (k in seq_len(nrow(D3))) for (g in GENES) {
  want <- norm_t3(t3[[D3$col[k]]][match(g, t3[[1]])])
  got  <- D$txt[D$ci == D3$ci[k] & as.character(D$gene) == g]
  if (got != want) stop("与 Table 3 不符：", g, " | ", D3$col[k], " | 图=", got, " 表=", want)
  n_t3 <- n_t3 + 1L
}
stopifnot(n_t3 == length(GENES) * nrow(D3))
cat(sprintf("与 Table 3 逐格对账：%d 格全部相同（%d 列 x %d 行）\n", n_t3, nrow(D3), length(GENES)))

# 正文 §3.3 与 §4 点名的字面值，单独钉一遍（这些是读者会去正文对照的）。
# 注意每一条都绑到**该层该估计量**的列上：换列会让数值看着合理而归属错误。
stopifnot(fp(n2(fam$disc_primary_P[fam$gene == "LOXL1"])) == "0.0037")   # §3.3 第三点·主估计
stopifnot(fp(n2(fam$disc_ivw_P[fam$gene == "LOXL1"]))     == "3.9e-04")  # §3.3 第三点·敏感性
stopifnot(fp(n2(fam$disc_ivw_P[fam$gene == "LOXL4"]))     == "3.0e-04")  # 比 LOXL1 小
stopifnot(n2(fam$disc_ivw_logFC[fam$gene == "LOXL4"]) < 0)               # 方向相反
stopifnot(fm(n2(fam$rep_primary_logFC[fam$gene == "LOXL1"])) == "+0.773")  # §4 主估计
stopifnot(fm(n2(fam$rep_pooled_logFC[fam$gene == "LOXL1"]))  == "+0.884")  # §4 策略比较
stopifnot(fp(n2(fam$rep_primary_P[fam$gene == "LOXL1"])) == "0.196")     # §4 主估计
stopifnot(fp(n2(fam$rep_pooled_P[fam$gene == "LOXL1"]))  == "0.0028")    # §4 敏感性
stopifnot(fp(n2(fam$rep_primary_P[fam$gene == "PLOD1"])) == "0.0032")    # §3.3 收束
# §3.3 收束句把discovery层 PLOD1 的同一个数写成 P=0.093（三位有效数字的散文化写法），
# 表里与图上都是 0.0929 —— 断言钉表值，散文化的 0.093 不作为判据。
stopifnot(fp(n2(fam$disc_primary_P[fam$gene == "PLOD1"])) == "0.0929")

# ── 图注/图内文字一律 ASCII：基础 pdf 设备经由 mbcsToSbcs 转 PostScript 编码，
#    非 ASCII 会触发 conversion failure 并**静默丢字**，而 ggsave 仍成功、PNG 却正常。
# 删除列的单元格是**两行**（log2FC 与 P），故换行符要放行；只拦可打印 ASCII 之外的字符。
# 写成「不属于 [空格..~] 且不是换行」比写白名单稳（后者会漏掉全角括号、U+2212 之类）。
lab_all <- unique(c(D$txt, as.character(D$col_lab), as.character(D$grp), as.character(D$ylab)))
bad <- lab_all[grepl("[^ -~\n]", lab_all)]
if (length(bad)) stop("含非 ASCII 字符，基础 pdf 设备会静默丢字：", paste(bad, collapse = " | "))

# ── 着色：四个 log2FC 列共用发散色标；其余按列内归一（见文件头）
KEY <- paste(
  "Shading, not a composite score: log2FC cells sit on one diverging scale shared by all four",
  "log2FC columns (both layers, both estimators; blue negative, red positive, limits +/-1.25);",
  "other cells are shaded within their own column -- P and deletion cells by -log10 P, I2 cells by",
  "0-100%, rank cells by position in the set (rank 1 darkest). Every value is printed; the shading",
  "carries no ordering across columns. * reported panel gene.  Deletion cells are labelled with the",
  "reduced-set log2FC (line 1) and its P (line 2); * marks a reduced-set effect opposite in sign to",
  "the full discovery estimate.  n/a = not computable (evaluable in too few cohorts), not zero.", sep = "\n")

LFC_MAX <- 1.25
# **只允许两个填充图层**，且顺序固定为 Dl（log2FC）在前、Do（其余）在后。
# 原因（2026-09-25 实测，ggplot2 4.0.3）：facet_grid(free_x) 下离散 scale 的面板断点序
# 由**图层出现顺序**决定，不是 factor 水平序。原先插入的第三个图层（n/a 格）在名次组里
# 只含 rep 两列 11、12，于是名次面板被排成 11 12 09 10 —— 列标签与数据整体错位，
# 而图仍照常生成、退出码 0。Dl 只含四个 log2FC 列（01、04、06、09），**每个都是其组内第一列**，
# 故不会扰动顺序；Do 的数据按列号升序，故其余各组的序天然正确。
# n/a 的灰由两个 scale 的 na.value 各自给出（两处都写死 grey93，避免同一行两种灰）。
Dl  <- D[D$kind == "lfc", ]
Do  <- D[D$kind != "lfc", ]
Do$fill01 <- NA_real_
w <- Do$kind %in% c("p", "del") & is.finite(Do$fill)   # 按列内归一（n/a 不参与）
Do$fill01[w] <- ave(Do$fill[w], Do$col_key[w], FUN = function(v) {
  r <- range(v, na.rm = TRUE)
  if (diff(r) < 1e-12) return(rep(0.55, length(v)))
  (v - r[1]) / diff(r)
})
Do$fill01[Do$kind == "i2"]   <- Do$fill[Do$kind == "i2"] / 100          # 自然量纲 0-100%
Do$fill01[Do$kind == "rank"] <- Do$fill[Do$kind == "rank"]              # 已是 0-1 百分位
ok <- is.finite(Do$fill01)
stopifnot(all(Do$fill01[ok] >= -1e-9 & Do$fill01[ok] <= 1 + 1e-9))
stopifnot(all(!ok == !is.finite(Do$fill)))   # 只有 n/a 格允许没有 fill01

p <- ggplot() +
  geom_tile(data = Dl, aes(x = col_key, y = ylab, fill = fill), colour = "grey80", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-LFC_MAX, LFC_MAX), oob = scales::squish,
                       na.value = "grey93", breaks = c(-1, 0, 1), name = "log2FC") +
  new_scale_fill() +
  geom_tile(data = Do, aes(x = col_key, y = ylab, fill = fill01), colour = "grey80", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "#33506B", na.value = "grey93",
                      breaks = c(0, 1), labels = c("low", "high"),
                      name = "within column") +
  geom_text(data = D, aes(x = col_key, y = ylab, label = txt), size = 1.95, lineheight = 0.86) +
  scale_x_discrete(labels = COL_LAB) +
  facet_grid(. ~ grp, scales = "free_x", space = "free_x", switch = "x") +
  labs(x = NULL, y = NULL, caption = KEY) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey93", colour = "grey70"),
        strip.text = element_text(face = "bold", size = 7.6),
        axis.text.x = element_text(size = 7.2, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.ticks = element_blank(),
        panel.spacing.x = unit(3.5, "pt"),
        legend.position = "right", legend.key.height = unit(0.62, "cm"),
        legend.key.width = unit(0.30, "cm"),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6.6),
        plot.caption = element_text(size = 6.1, hjust = 0, lineheight = 1.12, colour = "grey20"),
        plot.caption.position = "plot", plot.margin = margin(4, 8, 4, 4))

# ── 断言：渲染时的面板断点序必须就是 COLS 的列序。
#    单靠看图discovery不了：列整体错位后图仍"看着合理"（标签跟着数据一起动），退出码也仍是 0。
bd <- ggplot_build(p)
for (i in seq_len(nrow(bd$layout$layout))) {
  G <- as.character(bd$layout$layout$grp[i])
  want <- sprintf("%02d", which(vapply(COLS, function(C) C$grp, "") == G))
  got  <- as.character(bd$layout$panel_params[[i]]$x$get_breaks())
  if (!identical(got, want)) stop("面板 [", G, "] 的列序与 COLS 不符：\n  渲染=",
                                  paste(got, collapse = " "), "\n  应为=", paste(want, collapse = " "))
}
cat(sprintf("列序断言：%d 个面板的断点序均与 COLS 一致\n", nrow(bd$layout$layout)))

PDF <- file.path(OUT, "Fig_evidence_matrix.pdf"); PNG <- file.path(OUT, "Fig_evidence_matrix.png")
ggsave(PDF, p, width = 12.4, height = 5.6, device = pdf)
ggsave(PNG, p, width = 12.4, height = 5.6, dpi = 400)

for (f in c(PDF, PNG)) if (!file.exists(f)) stop("未产出：", f)
# 交付件名：E3 从 04_图表/ 取冻结件、以 Figure4_evidence_matrix.* 之名发出。先写冻结位、
# 再replication成交付件名并逐字节核对（各存一次必然漂移；见 E2/B3b/B7/D3 同款做法）。
for (pr in list(c(PDF, file.path(DELIV, "Figure4_evidence_matrix.pdf")),
                c(PNG, file.path(DELIV, "Figure4_evidence_matrix.png")))) {
  if (normalizePath(pr[1]) != normalizePath(pr[2], mustWork = FALSE)) {
    if (!file.copy(pr[1], pr[2], overwrite = TRUE)) stop("replication到交付件名失败：", pr[2])
    if (tools::md5sum(pr[1])[[1]] != tools::md5sum(pr[2])[[1]])
      stop("交付件与冻结件字节不同：", pr[2])
  }
}

raw <- readBin(PDF, "raw", file.size(PDF))
if (!identical(rawToChar(raw[1:5]), "%PDF-")) stop("不是 PDF：", PDF)
npage <- length(grepRaw("/Type /Page", raw, all = TRUE)) - 1L
if (npage < 1L) stop("PDF 无页面对象：", PDF)

# 审计 CSV：图上的每一格都落盘，事后可逐格核
AUD <- D[order(D$ci, match(as.character(D$gene), GENES)), c("grp", "col_lab", "gene", "txt")]
# 组名内嵌的换行只为图面排版，落盘时换成 " / " —— 否则 CSV 的字段里带裸换行，
# 下游按行读会把一条记录读成两条。
AUD$grp <- gsub("\n", " / ", as.character(AUD$grp))
names(AUD) <- c("column_group", "column", "gene", "cell_label")
write.csv(AUD, file.path(RES, "Fig4_values.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat(sprintf("已核验：%s（%d 页，%.1f KB）\n        %s（%.1f KB）\n",
            PDF, npage, file.size(PDF) / 1024, PNG, file.size(PNG) / 1024))
cat(sprintf("矩阵：%d 行 x %d 列（%d 组）；单元格 %d；n/a %d 格\n",
            length(GENES), length(COLS), nlevels(D$grp), nrow(D),
            sum(D$txt == "n/a")))
cat("n/a 明细：\n")
print(AUD[AUD$cell_label == "n/a", ], row.names = FALSE)
cat(sprintf("\n审计表：%s\n", file.path(RES, "Fig4_values.csv")))
print(AUD, row.names = FALSE)
