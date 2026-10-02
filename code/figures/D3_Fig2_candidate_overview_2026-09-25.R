# D3 图 2：57 个候选跨两层总览（稿件 Figure 2 / §3.1）
#
# 图注原文（现役主稿 English_manuscript_共用清洁稿_2026-09-28i.md L306）逐项落在数据列上。
# 本节把四项编码与图注**逐字绑定**，并把图注里的两处事实陈述（哪 5 个gene不可合并、
# 原因是无 GPL96 探针）回数据核验——图注是待画的规格书，不是可以目测放过的散文。
#
#   x    = discovery-layer REML+HK log2FC      → meta_effects_discovery.csv$primary_logFC（= hk_logFC）
#   y    = replication-layer REML+HK log2FC    → meta_effects_replication.csv$primary_logFC（= hk_logFC）
#          （两层同一估计量 = §2.3 声明的主估计；★ 2026-10-01 由 IVW-FE / pooled 改来，
#            又由 RE-DL 再改为 REML+HK）
#   色   = 两层主估计（REML+HK）是否**同号** → 本节现算，不读 B8 的 same_dir
#          （B8 那列按 primary_* 算，与本节的新坐标不是同一个估计量）
#          注意：这是**跨层**之比，不是层内逐队列投票一致性（后者另有 disc_unanimous/
#          rep_unanimous 两列，实测 24/33 与 24/28，与本通道不是同一个量）。
#          2026-09-29 回应 R3-m2：图例原文 `direction vs the layer's primary estimate`
#          会被读成层内比较，与主稿图注「concordant/opposite between the two layers'
#          primary estimates, not per-cohort vote consistency within a layer」冲突。
#          现图例标题改为 `Direction: discovery vs replication`、档名改 concordant/opposite，
#          与主稿图注用同一套词。
#   大小 = number of evaluable cohorts          → n_coh_total
#          （本脚本断言它等于逐队列表里有测值的队列数，而不是 6 或 3 这类常数）
#
# 纪律（本项目反复踩到的三条）：
#  ① **不可合并 ≠ 方向相反。** 图注要求二者在图上可区分。5 个无replication层值的gene画在
#     下方的独立条带里（"drawn on the discovery axis only"），与主面板的"方向相反"
#     不共用任何视觉通道；两者若同框同形，读图的人会把"测不了"读成"结论相反"。
#  ② **n/a 类缺失不得与"值 = 0"共用标记。** 条带上的点只承载 x，不画 y。
#  ③ **标签不得误导通道。** 图注只说"the six panel genes and the eight cross-linking
#     genes labelled"，未声明任何按gene类别区分的形状或字号，故本图**不**加该类编码——
#     加了就是图上有一条图注没交代的通道（见笔记 legend-title-shows-aes-variable-name）。
#     面孔gene与交联gene靠名字与正文对应，LOXL1 同时属于两类，只标一次。
#
# 设备：基础 pdf 设备（同 B7/D2）。行标签只有 ASCII —— 非 ASCII 会触发 mbcsToSbcs
# 的 conversion failure 并**静默丢字**，而 ggsave 仍返回成功。
#
# 输入：02_结果表/{57gene总览, 逐队列效应_力学gene, 合并效应_discovery层, 合并效应_replication层}.csv
#       01_审计/sample_manifest_per_sample.csv（平台列，用于核验图注的 GPL96 陈述）
# 输出：04_图表/Fig_candidate_overview.pdf / .png（冻结件位）
#       03_图表与源数据/Figure2_candidate_overview.pdf / .png（交付件名，逐字节replication自前者）
#       02_结果表/Fig2_values.csv

# ★ 2026-09-29 路径修正（与 B7 同款）。原 BASE 写死上一棵树
#   （…/滑膜基质刚度生成与刚度感知 0912/08_修订_2026-09-17）：本脚本跑出来的图落进旧树，
#   本树 04_图表/ 里的 Fig_candidate_overview.* 并不会被更新，交付目录自然也换不上新图。
#   实测形态：改完图例措辞后重跑，交付的 Figure2 仍是 09-27 的旧字节，而退出码是 0。
#   现改为从本脚本位置推树根（scripts/ → 01_冻结结果与代码），并把产物同时写进交付件位。
#
# ★ 2026-10-01 第十轮（回应审稿 P0-1）：**两轴改取 RE-DL。**
#   原图 x = discovery层 IVW-FE、y = replication层 pooled —— 两层不同估计量，与 §2.3 的新声明
#   （两层同一主估计 RE-DL）打架，主稿 Figure 2 图注（h.md:306）已先行改定。
#   实测两处硬编码**仍然成立、不必改**：RE-DL 口径下 x∈[-0.925, 0.792]、
#   y∈[-0.880, 0.774]，仍包在 XL=c(-1.02,1.02) 内（无点被裁）；两层方向相反的
#   候选仍是 16 个，故 `must(N_DIFF == 16L)` 与图内题注的 16 同时不变。
#   审计表加 estimator 列自陈口径（坐标列名不变，仍 x_disc_logFC/y_rep_logFC）。
#
# ★ 2026-10-01 第十一轮（第二轮审稿 P0-1 再改判）：**两轴改取 REML+HK（primary_*）。**
#   B2 已把 primary_* 指向 hk_*，故坐标列由 re_* 改为 primary_*；本文件内的数据列名
#   x_re/y_re 随之改名为 x_prim/y_prim（名字若留 re_ 会与实际取列脱钩）。
#   实测口径变动的后果**为零**：x∈[-0.9249, 0.7919]、y∈[-0.8798, 0.7728]，仍包在
#   XL=c(-1.02,1.02) 内（无点被裁，逐点已核）；两层方向相反的候选仍是 16 个，
#   故 `must(N_DIFF == 16L)` 与图内题注的 16 同时不变；52/5 的缺失形态也不变。
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
AUD <- .src
OUT <- .out
DELIV <- .out
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(OUT), dir.exists(RES), dir.exists(AUD))

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ggrepel)
})

must <- function(cond, msg) { if (!isTRUE(cond)) stop("✗ ", msg, call. = FALSE); cat("  OK  ", msg, "\n") }
say  <- function(...) cat(sprintf(...), "\n")

cat("D3 绘制图 2：57 个候选跨两层总览\n")

# ── 1. 读数据
ov  <- read.csv(file.path(RES, "candidate_overview_57.csv"), stringsAsFactors = FALSE)
per <- read.csv(file.path(RES, "per_cohort_effects.csv"), stringsAsFactors = FALSE)
dis <- read.csv(file.path(RES, "meta_effects_discovery.csv"), stringsAsFactors = FALSE)
rep <- read.csv(file.path(RES, "meta_effects_replication.csv"), stringsAsFactors = FALSE)

must(nrow(ov) == 57L, sprintf("候选 %d 个（期望 57）", nrow(ov)))
must(sum(ov$in_panel) == 6L, sprintf("面板gene %d 个（期望 6）", sum(ov$in_panel)))
must(sum(ov$is_crosslink) == 8L, sprintf("交联gene %d 个（期望 8）", sum(ov$is_crosslink)))

# ── 2. 把图注的四项编码回数据核验（不是核对 ov 自己，而是回到上游表）
must(all(is.finite(dis$primary_logFC) & is.finite(dis$primary_P)),
     sprintf("discovery层上游表 primary_logFC/primary_P 逐行有值（%d 行，= 图注所称 discovery-layer REML+HK）", nrow(dis)))
must(all(is.finite(rep$primary_logFC) & is.finite(rep$primary_P)),
     sprintf("replication层上游表 primary_logFC/primary_P 逐行有值（%d 行，= 图注所称 replication-layer REML+HK）", nrow(rep)))
# 主估计必须是 REML+HK，且 primary_* 确实取自 hk_*（挂了新名字却仍读旧列是最易犯的错）
must(all(dis$primary_method == "REML+HK") && all(rep$primary_method == "REML+HK"),
     "两层上游表 primary_method 均为 REML+HK")
must(all(dis$primary_logFC == dis$hk_logFC) && all(rep$primary_logFC == rep$hk_logFC),
     "上游表 primary_logFC 逐格等于 hk_logFC（主估计确实是 HK 那一列）")

# x/y 列：逐gene取自两层上游表的 primary_logFC（§2.3 声明的主估计）。
# ov$disc_logFC/ov$rep_logFC 由 B8 按 primary_* 生成，与本节同源；此处**只借它的缺失形态**
# （哪 5 个无replication层值），取值仍显式从上游表取 —— 免得将来 B8 换了列而这里静默跟着变。
i <- match(ov$gene, dis$gene)
must(!any(is.na(i)), "57 个候选全部出现在discovery层合并表中")
ov$x_prim <- dis$primary_logFC[i]
must(all(is.finite(ov$x_prim)),
     sprintf("x 轴取自discovery层 primary_logFC，%d 个候选全部有值", nrow(ov)))
j <- match(ov$gene, rep$gene)
has_y <- !is.na(ov$rep_logFC)
must(all(!is.na(j[has_y])), "凡有replication层值的gene都出现在replication层合并表中")
ov$y_prim <- rep$primary_logFC[j]
must(identical(is.na(ov$y_prim), is.na(ov$rep_logFC)),
     "不可合并的 5 个与 B8 的 rep_logFC 缺失形态逐一相同（不是另一套缺失）")
must(all(is.finite(ov$y_prim[has_y])), "凡有replication层值的gene，其 primary_logFC 亦有值")

# 大小列：n_coh_total 必须是"有测值的队列数"，故先断言逐队列表无 NA，再逐gene数行
must(!any(is.na(per$logFC)), "逐队列效应表无 NA 值（故行数 = 有测值的队列数）")
n_by_gene <- table(per$gene)[ov$gene]
must(all(as.integer(n_by_gene) == ov$n_coh_total),
     sprintf("大小列 n_coh_total 逐gene等于逐队列表行数（取值 %s）",
             paste(sort(unique(ov$n_coh_total)), collapse = "/")))

# 色列：按**新的主估计**（REML+HK）重算，不可评估为 NA（不是 FALSE）。
# B8 落盘的 same_dir 是按旧的 primary_* 算的，故不读它、只重算 ——
# 读旧列会让颜色与坐标来自两个估计量。
same_dir <- sign(ov$x_prim) == sign(ov$y_prim)
same_dir[is.na(ov$y_prim)] <- NA
must(identical(is.na(same_dir), is.na(ov$y_prim)),
     "色列不可评估者为 NA 而非 FALSE（与条带的缺失形态一致）")
must(!any(is.na(same_dir[!is.na(ov$y_prim)])), "凡两层都有值的gene，方向一致性都可判定")
say("  方向一致性：一致 %d / 相反 %d / 不可评估 %d（合计 %d）",
    sum(same_dir, na.rm = TRUE), sum(!same_dir, na.rm = TRUE), sum(is.na(same_dir)), nrow(ov))

# ── 3. 核验图注的两处事实陈述
LEGEND_NC <- c("LOXL4", "LATS2", "TEAD2", "LOXL3", "MOB1B")     # 图注点名
nc <- ov$gene[is.na(ov$rep_logFC)]
must(setequal(nc, LEGEND_NC),
     sprintf("无replication层值者恰为图注点名的 5 个（实得 %s）", paste(sort(nc), collapse = ", ")))
must(all(is.na(ov$rep_k[is.na(ov$rep_logFC)])),
     "这 5 个gene的 rep_k 亦为 NA（未进replication层合并表，而非进了表却无值）")

# 图注称原因是"no GPL96 probe"：回样本清单取平台，再查这 5 个gene在这些队列有无测值
sm <- read.csv(file.path(AUD, "sample_manifest_per_sample.csv"), stringsAsFactors = FALSE)
plat <- tapply(sm$platform, sm$dataset, function(x) unique(x)[1])
must(identical(sort(names(plat)[plat == "GPL96"]), c("GSE55235", "GSE55457")),
     sprintf("GPL96 队列 = %s（回样本清单实测）", paste(sort(names(plat)[plat == "GPL96"]), collapse = ", ")))
gpl96 <- names(plat)[plat == "GPL96"]
absent <- vapply(LEGEND_NC, function(g) sum(per$gene == g & per$dataset %in% gpl96) == 0L, TRUE)
must(all(absent),
     sprintf("这 5 个gene在两个 GPL96 队列均无测值（故replication层只剩 GSE82107 一个队列，无法合并）——图注的因果陈述成立"))

# ── 4. 准备绘图数据
XL <- c(-1.02, 1.02)                       # 两轴同域：位置读作效应大小而非名次
main  <- ov[!is.na(ov$rep_logFC), ]
main$dir_lab <- ifelse(same_dir[!is.na(ov$rep_logFC)], "concordant", "opposite")
strip <- ov[is.na(ov$rep_logFC), ]
# 受标者 = 面板 ∪ 交联，**按gene身份筛，不按面板位置筛**。5 个不可合并gene里有 2 个
# （LOXL4、LOXL3）是交联gene，故受标点分散在主面板 11 个 + 条带 2 个；条带另 3 个
# （LATS2/TEAD2/MOB1B）既非面板也非交联，只作"不可合并"的点出现、不标名。
# 原先写 `lab <- main[main$in_panel | main$is_crosslink, ]` 再加 `nrow(lab)==13`，
# 因为 main 已排除那 5 个，断言当场失败——这条断言正是用来拦这种"按位置筛"的。
ID   <- ov$in_panel | ov$is_crosslink
lab   <- main[main$in_panel | main$is_crosslink, ]          # 主面板内 11 个
labs  <- strip[strip$in_panel | strip$is_crosslink, ]        # 条带内 2 个（LOXL4、LOXL3）

must(nrow(main) == 52L && nrow(strip) == 5L,
     sprintf("主面板 %d 个（两层皆可评）、条带 %d 个（仅discovery层）", nrow(main), nrow(strip)))
must(nrow(lab) + nrow(labs) == 13L,
     sprintf("受标gene合计 %d 个 = 主面板 %d + 条带 %d（6 面板 + 8 交联，LOXL1 两类兼属故去重）",
             nrow(lab) + nrow(labs), nrow(lab), nrow(labs)))
must(nrow(labs) == 2L && setequal(labs$gene, c("LOXL4", "LOXL3")),
     sprintf("条带内受标者恰为交联gene LOXL4、LOXL3（实得 %s）", paste(sort(labs$gene), collapse = ", ")))
must(all(!grepl("[^ -~]", c(main$gene, lab$gene, labs$gene))),
     "受标gene名全部 ASCII（基础 pdf 设备非 ASCII 会静默丢字）")

# 大小档位取自**主面板实测值**（定义须在 p_main 之前）
SZ_BREAKS <- sort(unique(main$n_coh_total))

COLS <- c(concordant = "#2E6E9E", opposite = "#C1543F")
TH <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        legend.key.height = unit(9, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(4, 6, 2, 2))

# ── 5. 主面板
# 建成函数而非一次赋值：负面对照要「用改写后的数据重建同一张图」。
# 若只写 `main$dir_lab <- "concordant"` 再 build 原对象，R 的写时replication使图仍指向旧数据，
# 对照会静默空转、恒真通过（本轮已实测踩到）。
build_main <- function(d) ggplot(d, aes(x = x_prim, y = y_prim)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3, linetype = "22") +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3, linetype = "22") +
  geom_point(aes(colour = dir_lab, size = n_coh_total), alpha = 0.85) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 2.3,
                           seed = 20260925, min.segment.length = 0.2,
                           segment.size = 0.25, segment.colour = "grey45",
                           box.padding = 0.28, point.padding = 0.18,
                           max.overlaps = Inf, segment.alpha = 0.8) +
  # 图例标题必须自陈是**哪两者**之比：只说 "the layer's primary estimate" 会被读成层内比较
  scale_colour_manual(values = COLS, name = "Direction: discovery vs\nreplication",
                      labels = c(concordant = "concordant", opposite = "opposite")) +
  # 大小档位取自**主面板实测值**，不写死 c(3,5,6)：取值为 3 的恰是那 5 个不可合并
  # gene，它们不在主面板，写死会画出一个没有任何点的图例项（图例声称一档，图上无此档）。
  scale_size_continuous(range = c(1.7, 4.4), breaks = SZ_BREAKS,
                        name = "Evaluable cohorts") +
  coord_cartesian(xlim = XL, ylim = XL) +
  labs(x = "Discovery-layer log2FC (REML+HK)", y = "Replication-layer log2FC (REML+HK)") +
  guides(colour = guide_legend(order = 1), size = guide_legend(order = 2)) +
  TH
set.seed(20260925)
p_main <- build_main(main)

# ── 6. 不可合并条带（只承载 x，y 不承载任何量）
# **两行错位是版面装置，不是第二个量。** LATS2(-0.364) 与 TEAD2(-0.357) 的 x 只差 0.007
# （面板宽度的 0.3%），单行时两个三角完全重合——5 个gene只看得到 4 个点，而图注说这里
# 有 5 个。错位按 x 排序后轮流分两行，确定性、可复现；y 在此条带不承载任何数值。
strip <- strip[order(strip$x_prim), ]
strip$row_y <- ifelse(seq_len(nrow(strip)) %% 2L == 1L, 0.17, -0.17)
# labels 子集**必须在此之后重取**：labs 早先是 strip 的子集，那时还没有 row_y 列，
# 而本层继承 plot 的 aes(y = row_y)，在旧子集里找不到该列会报 "object 'row_y' not found"。
labs <- strip[strip$in_panel | strip$is_crosslink, ]
stopifnot(nrow(labs) == 2L, all(c("row_y") %in% names(labs)))
p_strip <- ggplot(strip, aes(x = x_prim, y = row_y)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3, linetype = "22") +
  geom_point(shape = 2, size = 2.5, colour = "#6B6B6B", stroke = 0.7) +
  ggrepel::geom_text_repel(data = labs, aes(label = gene), size = 2.3, seed = 20260925,
                           direction = "x", nudge_y = -0.34,
                           segment.size = 0.25, segment.colour = "grey45",
                           box.padding = 0.22, max.overlaps = Inf) +
  coord_cartesian(xlim = XL, ylim = c(-0.5, 0.5), expand = FALSE) +
  labs(x = NULL, y = NULL) +
  annotate("text", x = XL[1], y = -0.44, hjust = 0, vjust = 0, size = 2.2,
           colour = "grey25",
           label = "Not combinable in the replication layer (x only; rows separate coincident points)") +
  TH +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "grey60", linewidth = 0.3))

N_DIFF <- sum(!same_dir, na.rm = TRUE)
must(N_DIFF == 16L, sprintf("两层方向相反的候选 %d 个（图注与正文所称的 16）", N_DIFF))

CAP <- paste(
  "52 of the 57 candidates are evaluable in both layers; 5 have no replication-layer value.",
  "Both axes are the layer's declared primary estimator (REML with the Hartung-Knapp variance",
  "correction; section 2.3). Point colour encodes whether the two layers' primary",
  "estimates agree in sign (concordant or opposite), not per-cohort vote consistency within a",
  "layer; point size is the number of evaluable cohorts.",
  sprintf("The 5 candidates below carry no replication value and are not combinable; this is distinct from the %d candidates whose direction differs between layers.", N_DIFF),
  "No rank is implied by position.")
# 手工折行：plot_annotation 的 caption 不自动换行，单行会溢出画布右缘被截断（本轮已踩到）。
CAP <- paste(strwrap(CAP, width = 96), collapse = "\n")

# ── 6b. 图注自检：五句必须逐句在位
# 起因：上一轮用「源脚本里搜子串」的办法核图注，第 4 句报 ✗。该句由 sprintf 在运行时拼出，
# 源码里根本不存在那个字面串——是检索方法的假阳性，不是图注缺字。判据必须落在**运行时的
# 字符串**上，故在这里就地断言，此后每次运行都自证，不依赖事后对源码做子串搜索。
CAP_FLAT <- gsub("\n", " ", CAP, fixed = TRUE)
CAP_MUST <- c(
  "52 of the 57 candidates are evaluable in both layers; 5 have no replication-layer value.",
  "Both axes are the layer's declared primary estimator (REML with the Hartung-Knapp variance",
  "correction; section 2.3). Point colour encodes whether the two layers' primary",
  "estimates agree in sign (concordant or opposite), not per-cohort vote consistency within a",
  "layer; point size is the number of evaluable cohorts.",
  sprintf(paste("The 5 candidates below carry no replication value and are not combinable;",
                "this is distinct from the %d candidates whose direction differs between layers."),
          N_DIFF),
  "No rank is implied by position.")
for (s in CAP_MUST) {
  must(grepl(s, CAP_FLAT, fixed = TRUE),
       sprintf("图注在位（%d 字）：%s…", nchar(s), substr(s, 1, 52)))
}

p <- p_main / p_strip +
  plot_layout(heights = c(6.2, 1)) +
  plot_annotation(caption = CAP,
    theme = theme(plot.caption = element_text(size = 6.2, hjust = 0, colour = "grey25",
                                              lineheight = 1.15, margin = margin(t = 3))))

# ── 7. 断言：编码与图注一致（画完再验图上通道，而不是只看脚本写了什么）
# 按 geom 的**类名**定位散点层，不按下标猜，也不按列名猜——ggrepel 的文本层同样带
# colour/size/x/y 四列，按列名定位会命中 2 个图层（本轮即如此）。
pt_layer <- function(g) which(vapply(g$layers, function(L) inherits(L$geom, "GeomPoint"), TRUE))
bd <- ggplot_build(p_main)
k <- pt_layer(bd$plot)
must(length(k) == 1L, sprintf("主面板恰有 %d 个散点层（按 geom 类名定位）", length(k)))
pts <- bd$data[[k]]
must(identical(sort(unique(as.character(pts$colour))), sort(unname(COLS))),
     sprintf("主面板实际用到的颜色恰为图注所述的两类（%s）",
             paste(sort(unique(as.character(pts$colour))), collapse = ", ")))
must(length(unique(pts$size)) == length(SZ_BREAKS) && length(SZ_BREAKS) >= 2L,
     sprintf("点大小恰 %d 档、与图例所报档位一一对应（%s）—— 大小通道非恒定，也无空档",
             length(unique(pts$size)), paste(SZ_BREAKS, collapse = "/")))
must(identical(SZ_BREAKS, c(5L, 6L)),
     sprintf("主面板 n_coh_total 实测档位 %s（取值为 3 的 5 个gene全在条带，故主面板无此档）",
             paste(SZ_BREAKS, collapse = "/")))
# 反向守卫：若把 same_dir 全置为一致，颜色应退化为 1 类；据此证明上面那条不是恒真
# 负面对照：用「dir_lab 全置一致」的数据重建同一张图（走同一个 build_main），
# 颜色应退化为 1 类。改写作用在传给 build_main 的副本上，故必然生效。
main_neg <- main; main_neg$dir_lab <- "concordant"
bd_neg <- ggplot_build(build_main(main_neg))
neg_k <- pt_layer(bd_neg$plot)
must(length(neg_k) == 1L && length(unique(as.character(bd_neg$data[[neg_k]]$colour))) == 1L,
     sprintf("负面对照：dir_lab 全置一致后颜色退化为 %d 类（证明颜色断言非恒真）",
             length(unique(as.character(bd_neg$data[[neg_k]]$colour)))))

must(nrow(strip) == 5L && length(unique(strip$row_y)) == 2L,
     sprintf("条带 5 个点分 %d 行错位（避免 x 相同者完全重合）", length(unique(strip$row_y))))
dmin <- min(abs(dist(strip$x_prim)))
must(dmin < 0.05,
     sprintf("错位确有必要：未错位时最近两点 x 仅差 %.4f 个 log2FC（面板宽的 %.2f%%）",
             dmin, 100 * dmin / diff(XL)))

# ── 8. 导出
PDF <- file.path(OUT, "Fig_candidate_overview.pdf")
PNG <- file.path(OUT, "Fig_candidate_overview.png")
ggsave(PDF, p, width = 6.6, height = 7.4, device = pdf)
ggsave(PNG, p, width = 6.6, height = 7.4, dpi = 400)

# 交付件名：E3 从 04_图表/ 取冻结件、以 Figure2_candidate_overview.* 之名发出。先写冻结位、
# 再replication成交付件名并逐字节核对（各存一次必然漂移；见 E2/B3b/B7 同款做法）。
for (pr in list(c(PDF, file.path(DELIV, "Figure2_candidate_overview.pdf")),
                c(PNG, file.path(DELIV, "Figure2_candidate_overview.png")))) {
  if (normalizePath(pr[1]) != normalizePath(pr[2], mustWork = FALSE)) {
    if (!file.copy(pr[1], pr[2], overwrite = TRUE)) stop("replication到交付件名失败：", pr[2])
    if (tools::md5sum(pr[1])[[1]] != tools::md5sum(pr[2])[[1]])
      stop("交付件与冻结件字节不同：", pr[2])
  }
}

raw <- readBin(PDF, "raw", 5L)
must(identical(rawToChar(raw), "%PDF-"), "PDF 魔数 %PDF- 正确")
txt <- readBin(PDF, "raw", file.size(PDF))
must(sum(txt == as.raw(0x50)) > 0, "PDF 非空")
say("已核验：%s（%.1f KB）", PDF, file.size(PDF) / 1024)
say("        %s（%.1f KB）", PNG, file.size(PNG) / 1024)

# ── 9. 落盘审计表（图上的每个点都要能追回数值）
aud <- data.frame(gene = ov$gene, in_panel = ov$in_panel, is_crosslink = ov$is_crosslink,
                  estimator = "REML+HK",
                  x_disc_logFC = ov$x_prim, y_rep_logFC = ov$y_prim,
                  colour_class = ifelse(is.na(ov$rep_logFC), "not combinable",
                                        ifelse(same_dir, "concordant", "opposite")),
                  size_n_coh = ov$n_coh_total, panel = ifelse(is.na(ov$rep_logFC), "strip", "main"),
                  stringsAsFactors = FALSE)
AUDIT <- file.path(RES, "Fig2_values.csv")
# 不变量：源数据里的色类措辞 == 图例档名（并集的第三类 "not combinable" 是条带的独立标识，
# 图上另有一段文字说明，不在色标里）。两处措辞若各改各的，这里就炸 —— B7 同款纪律。
must(setequal(setdiff(unique(aud$colour_class), "not combinable"), names(COLS)),
     sprintf("源数据 colour_class 的措辞集合 == 图例档名（不含条带专有的 not combinable）：CSV %s / 图例 %s",
             paste(sort(unique(aud$colour_class)), collapse = ","),
             paste(sort(names(COLS)), collapse = ",")))
write.csv(aud, AUDIT, row.names = FALSE)
say("审计表 %s（%d 行）", basename(AUDIT), nrow(aud))
say("  主面板 %d 行 / 条带 %d 行；色类计数 %s", sum(aud$panel == "main"), sum(aud$panel == "strip"),
    paste(names(table(aud$colour_class)), table(aud$colour_class), sep = "=", collapse = " "))
