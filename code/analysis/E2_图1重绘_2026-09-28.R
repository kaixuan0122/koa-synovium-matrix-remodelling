#!/usr/bin/env Rscript
# =============================================================================
# E2 — 图 1「研究设计与四层证据」重绘（2026-09-28）
#
# 本脚本由 D4_图1研究设计_2026-09-27.R **逐字复制后做四处精确替换**而成，D4 保持原状。
# 重绘起因：D4 从**旧稿**（05_稿件/English_manuscript_MDPI_重定位2026-09-23.md）读图注，
# 而终稿（02_共用主稿/…2026-09-28b.md）的 Figure 1 图注已删去 functional enrichment，
# 图内嵌题注与终稿图注因此不再逐字相同；另有「for the first time」与四类口径两项已废内容。
#
# 四处替换：
# ① 图注来源改指终稿（并加断言：终稿图注不含 functional enrichment）
# ② 括注 NOTE 只承诺机器学习（富集模块未交付任何补充材料，工作令 §六）
# ③ 删「compared under one framework for the first time」中的首次声明（§四）
# ④ 「4 cell types」改「5 cell types」（§2.5 与图 5 的口径均为五类，含 mast cells）
#
# 输出：04_图表/Fig_study_design.pdf / .png（冻结件位，E3 从这里取图）
#       03_图表与源数据/Figure1_study_design.pdf / .png（交付件名，逐字节复制自前者）
#       02_结果表/图1_研究设计_布局.csv、图1_研究设计_箱体.csv（审计表）
#       02_结果表/图1_研究设计_图注.txt、图1_研究设计_括注.txt（D4b 的比对基准）
#
# ★ 2026-09-29 修正输出路径（本脚本此前从未成功运行，这是原因之一）：
#   本脚本由 D4 复制而来时，把 TAB 与 FIG 一并改到了 03_图表与源数据/，于是
#     ① 图写到 03_图表与源数据/Fig_study_design.*，而 E3 是从 04_图表/Fig_study_design.*
#        取图的 —— 即便本脚本跑通，交付包里换上的仍是 D4 的旧图；
#     ② 审计表写到交付目录，而 E3 的源数据是从 02_结果表/ 取的 —— 图与源数据会各说各话。
#   现改为：审计表回 02_结果表/（D4b 与 E3 的正位），图先写 04_图表/ 再复制成交付件名并
#   逐字节核对。D4 原产物保留不动，仅被本脚本覆盖 04_图表/ 下同名的一对。
#
# 版式（用户 2026-09-25 裁定）：**横向四层阶梯**。左起数据层，向右依次四个证据层，
# 层间箭头；底部一条通栏括注「补充材料：机器学习 · 富集」。
#
# 三条纪律，逐条落在断言上：
#
# ① **阶梯的纵向位移是一条通道，故必须在图注里被声明。** 稿件 Figure 1 图注末句即
#    "The four layers are ordered by how strongly each bounds the priority claim, not by the
#    order in which the analyses were run."，故本图把这句所在的**整条图注**逐字放进图内
#    （运行时从稿件读取，单一来源，不会漂移）。除此之外不引入任何未被图注声明的通道。
#
# ② **基线 pdf 设备遇非 ASCII 会静默丢字**（`mbcsToSMBcs` 转换失败而退出码仍为 0）。
#    图注含一个 EM DASH，须显式映射为 ASCII；映射须覆盖实测到的**全部**非 ASCII 字符，
#    漏一个就静默丢一个。所有上图的字符串都要过 `!grepl("[^ -~]", s)`。
#
# ③ **箱体尺寸由实测字宽决定，且宽度上限要独立于箱体存在。** 箱子等宽（版面约束），
#    文字必须装进去——若把箱宽当自变量去凑文字，那就成了恒真断言。故：箱宽由画布宽度
#    和列数算出，再逐行**实测**字宽并断言装得下；负面对照把箱宽减半，断言检查会失败。
#
# 数值来源（两处独立核对，非同一来源）：
#   样本清单 01_审计/样本清单_逐样本.csv  —— 逐样本重新计数
#   稿件 Figure 1 图注里的 "74 samples: 34 OA / 40 normal" 等 —— 解析后与上者比对
# 二者一致才继续；不一致即失败（本项目有过手写 n 与脚本脱节的先例）。
#
# 输出：04_图表/Fig_study_design.pdf / .png；02_结果表/图1_研究设计_布局.csv（审计表）
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(grid) })

ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "..", ".."))
REV  <- ROOT
AUD  <- file.path(REV, "01_审计")
TAB  <- file.path(REV, "02_结果表")      # 审计表的正位：D4b 与 E3 都从这里读
FIG  <- file.path(REV, "04_图表")        # 冻结件位：E3_图件交付组装 从这里取图
DELIV<- file.path(REV, "..", "03_图表与源数据")   # 交付件位（Figure1_study_design.*）
MS_STEM <- "English_manuscript_共用清洁稿_2026-09-28i"   # 与 E3/H1/F7/G1/H2/H3/D2 同一常量名
MS   <- file.path(REV, "..", "02_共用主稿", paste0(MS_STEM, ".md"))

FAIL <- 0L
must <- function(cond, msg) {
  if (!isTRUE(cond)) { FAIL <<- FAIL + 1L; cat(sprintf("  FAIL %s\n", msg)) }
  else cat(sprintf("  OK   %s\n", msg))
}
cat("E2 图 1 研究设计与四层证据（重绘版）\n\n")

# ── 0. 版面常数（英寸；与 ggsave 的 width/height 同一坐标系，靠 coord_fixed 1:1 映射）
WT      <- 7.09      # 180 mm，MDPI 双栏上限
MARGIN  <- 0.10
GAP     <- 0.26      # 列间距：层间箭头就活在 6.6 mm 的缝里，故不能太窄
NCOL    <- 5L        # 数据层 + 四个证据层
W_BOX   <- (WT - 2 * MARGIN - (NCOL - 1) * GAP) / NCOL
PAD     <- 0.062      # 须 ≥ 半行盒（7.0 pt 标题 x 1.2 / 2 = 0.058 in）；下有一条断言在守
STEP    <- 0.30      # 阶梯每级下沉量
LINE_H  <- 0.115
GAP_H   <- 0.55 * LINE_H
DEV_PS  <- 12
SZ      <- 6.2       # 正文磅值
SZ_HEAD <- 7.0
SZ_CAP  <- 6.0
PT2MM   <- 25.4 / 72.27   # geom_text 的 size 单位是 mm

# 实测字宽（英寸）。字宽只取决于磅值与字体，与设备长宽无关，故可在临时设备上量。
tmp <- tempfile(fileext = ".pdf")
pdf(tmp, width = WT, height = 4, pointsize = DEV_PS)
tw <- function(s, sz) as.numeric(strwidth(s, units = "in", cex = sz / DEV_PS))

# ── 1. 数据层：从逐样本清单重算，并与稿件图注里的自报数对账
sm <- read.csv(file.path(AUD, "样本清单_逐样本.csv"), stringsAsFactors = FALSE)
must(all(c("dataset", "role", "assigned_group", "included") %in% names(sm)),
     "样本清单含 dataset/role/assigned_group/included 四列")
inc <- sm[sm$included %in% c("TRUE", "True", "true"), ]
must(nrow(inc) == 131L, sprintf("清单在册样本 %d 例（全部为 OA/Normal，RA 已排除）", nrow(inc)))

cnt <- function(role) {
  s <- inc[inc$role == role, ]
  c(n = nrow(s), oa = sum(s$assigned_group == "OA"), no = sum(s$assigned_group == "Normal"))
}
D <- cnt("发现"); R <- cnt("复制")
must(D[["n"]] == 74L && D[["oa"]] == 34L && D[["no"]] == 40L,
     sprintf("发现层 n=%d（%d OA / %d normal）", D[["n"]], D[["oa"]], D[["no"]]))
must(R[["n"]] == 57L && R[["oa"]] == 30L && R[["no"]] == 27L,
     sprintf("复制层 n=%d（%d OA / %d normal）", R[["n"]], R[["oa"]], R[["no"]]))

# 逐队列（名次按清单出现序，不按字母序——图上的顺序要与分层一致）
cohorts <- function(role) {
  s <- inc[inc$role == role, ]
  u <- unique(s$dataset)
  data.frame(dataset = u, n = as.integer(table(factor(s$dataset, u))), stringsAsFactors = FALSE)
}
cd <- cohorts("发现"); cr <- cohorts("复制")
must(nrow(cd) == 3L && nrow(cr) == 3L,
     sprintf("发现 %d 个队列、复制 %d 个队列", nrow(cd), nrow(cr)))
plat <- tapply(sm$platform, sm$dataset, function(x) unique(x)[1])
want <- c(GSE89408 = "RNA-seq (count matrix)", GSE206848 = "GPL570", GSE1919 = "GPL91",
          GSE55235 = "GPL96", GSE55457 = "GPL96", GSE82107 = "GPL570")
# 注意 `tapply` 返回的是带 dim 的一维数组，`unname()` 不会去掉 dim，`identical()` 会因此
# 误判为不等——须用 as.vector() 取纯向量再比。
must(identical(as.vector(plat[names(want)]), as.vector(want)),
     "六个队列的平台与稿件 §2.1 逐个相符（GPL96 只落在 GSE55235/GSE55457 上）")
must(sum(plat == "GPL96") == 2L && setequal(names(plat)[plat == "GPL96"], c("GSE55235", "GSE55457")),
     "GPL96 恰两个队列——图 2 图注「no GPL96 probe」的因果陈述即据此")

# 稿件图注里的自报数（独立来源）
ms_lines <- readLines(MS, encoding = "utf-8", warn = FALSE)
leg <- ms_lines[grepl("^\\*\\*Figure 1:", ms_lines)]
must(length(leg) == 1L, "稿件中恰有一条 Figure 1 图注")
leg <- leg[[1]]
grab <- function(pat) {
  m <- regmatches(leg, regexec(pat, leg))[[1]]
  must(length(m) == 4L, sprintf("图注可解析出 `%s`", pat))
  as.integer(m[2:4])
}
d_leg <- grab("Discovery cohorts \\(three synovial cohorts, (\\d+) samples: (\\d+) OA / (\\d+) normal\\)")
r_leg <- grab("replication cohorts \\(three synovial cohorts, (\\d+) samples: (\\d+) OA / (\\d+) normal\\)")
must(identical(as.integer(d_leg), as.integer(c(D[["n"]], D[["oa"]], D[["no"]]))),
     sprintf("发现层：(n, OA, normal) = (%d, %d, %d)，清单重算与图注自报一致", d_leg[1], d_leg[2], d_leg[3]))
must(identical(as.integer(r_leg), as.integer(c(R[["n"]], R[["oa"]], R[["no"]]))),
     sprintf("复制层：(n, OA, normal) = (%d, %d, %d)，清单重算与图注自报一致", r_leg[1], r_leg[2], r_leg[3]))

# ── 2. 图注：非 ASCII 显式映射（EM DASH → ASCII），映射须覆盖实测到的全部非 ASCII
NA_MAP <- c("—" = "-")   # em dash
na_found <- unique(strsplit(leg, "")[[1]][grepl("[^ -~]", strsplit(leg, "")[[1]])])
must(all(na_found %in% names(NA_MAP)),
     sprintf("图注里的非 ASCII 字符 %d 种，全部在映射表内：%s",
             length(na_found), paste(sprintf("U+%04X", utf8ToInt(paste(na_found, collapse = ""))), collapse = " ")))
for (k in names(NA_MAP)) leg <- gsub(k, NA_MAP[[k]], leg, fixed = TRUE)
must(!grepl("[^ -~]", leg), "图注经映射后为纯 ASCII")

TITLE <- sub("^\\*\\*(.+?)\\*\\*.*$", "\\1", leg)
BODY  <- sub("^\\*\\*.+?\\*\\*[[:space:]]*", "", leg)
must(nzchar(TITLE) && nzchar(BODY) && !grepl("\\*\\*", BODY), "图注的题名与正文已拆开，且无残留 ** 标记")

# 换行按**实测字宽**做，不按字符数（比例字体下两者不等）
wrap_measured <- function(text, maxw, sz) {
  words <- strsplit(text, "[[:space:]]+")[[1]]
  out <- character(0); cur <- ""
  for (w in words) {
    cand <- if (nzchar(cur)) paste(cur, w) else w
    if (tw(cand, sz) <= maxw || !nzchar(cur)) cur <- cand else { out <- c(out, cur); cur <- w }
  }
  c(out, cur)
}
CAP_W  <- WT - 2 * MARGIN
BODY_LINES <- wrap_measured(BODY, CAP_W, SZ_CAP)
cap_lines  <- c(TITLE, BODY_LINES)   # 题名与正文分行，题名在图内加粗
CAP_H  <- length(cap_lines) * SZ_CAP * 1.15 / 72

# ── 3. 箱体内容。行类型：head 标题 / full 整行 / pair 左右两端 / gap 空档
rows_data <- c(list(list(k = "head", t = "DATA"),
                    list(k = "pair", l = "Discovery", v = sprintf("n = %d", D[["n"]])),
                    list(k = "full", t = sprintf("%d OA / %d normal", D[["oa"]], D[["no"]]))),
               lapply(seq_len(nrow(cd)), function(i)
                 list(k = "full", t = sprintf("%s %s", cd$dataset[i],
                      ifelse(plat[cd$dataset[i]] == "RNA-seq (count matrix)", "RNA-seq", "array")))),
               list(list(k = "pair", l = "Replication", v = sprintf("n = %d", R[["n"]])),
                    list(k = "full", t = sprintf("%d OA / %d normal", R[["oa"]], R[["no"]]))),
               lapply(seq_len(nrow(cr)), function(i)
                 list(k = "full", t = sprintf("%s %s", cr$dataset[i],
                      ifelse(plat[cr$dataset[i]] == "RNA-seq (count matrix)", "RNA-seq", "array")))),
               list(list(k = "gap"), list(k = "full", t = "DE per cohort"),
                    list(k = "full", t = "DESeq2 RNA-seq"), list(k = "full", t = "limma array")))

ROWS <- list(
  rows_data,
  list(list(k = "head",  t = "1  Candidate scope"),
       list(k = "full",  t = "57 candidates,"),
       list(k = "full",  t = "13 categories,"),
       list(k = "full",  t = "annotation-driven,"),
       list(k = "full",  t = "not data-driven"),
       list(k = "gap"),
       list(k = "full",  t = "7 genes tested:"),
       list(k = "full",  t = "LOXL1 ITGA5 SDC4 SAV1"),
       list(k = "full",  t = "VGLL4 CCN1 SERPINE1"),
       list(k = "full",  t = "CCN1/CYR61 withdrawn"),
       list(k = "full",  t = "-> 6 reported"),
       list(k = "full",  t = "alpha = 0.05/7 = 0.00714"),
       list(k = "full",  t = "8 cross-linkers;"),
       list(k = "full",  t = "only LOXL1 in panel")),
  list(list(k = "head",  t = "2  Repeatability"),
       list(k = "full",  t = "Per-cohort effects"),
       list(k = "full",  t = "with 95% CI"),
       list(k = "full",  t = "Fixed- and random-"),
       list(k = "full",  t = "effects meta-analysis"),
       list(k = "full",  t = "Leave-one-cohort-out"),
       list(k = "full",  t = "re-ranking"),
       list(k = "gap"),
       list(k = "full",  t = "5 of 57 candidates"),
       list(k = "full",  t = "not combinable in"),
       list(k = "full",  t = "replication (no"),
       list(k = "full",  t = "GPL96 probe)")),
  list(list(k = "head",  t = "3  Cross-linker family"),
       list(k = "full",  t = "8 cross-linking genes"),
       list(k = "full",  t = "compared under one"),
       list(k = "full",  t = "framework; ranked by"),
       list(k = "full",  t = "|log2FC| and by P"),
       list(k = "gap"),
       list(k = "full",  t = "No composite score")),
  list(list(k = "head",  t = "4  Cellular source"),
       list(k = "full",  t = "Donor-resolved"),
       list(k = "full",  t = "single-cell"),
       list(k = "full",  t = "localisation"),
       list(k = "full",  t = "(GSE152805; 3 OA"),
       list(k = "full",  t = "donors; 15 clusters;"),
       list(k = "full",  t = "5 cell types)"),
       list(k = "gap"),
       list(k = "full",  t = "Exploratory cis-eQTL"),
       list(k = "full",  t = "x FinnGen R12 KOA:"),
       list(k = "full",  t = "no transferable support"))
)

h_of <- function(rows) 2 * PAD + sum(vapply(rows, function(r)
  switch(r$k, head = LINE_H * 1.25, gap = GAP_H, LINE_H), 0))
H_BOX <- vapply(ROWS, h_of, 0)

# ── 4. 画布高度：由最低箱底 + 括注 + 图注算出（文字行数已知，故不必先量字宽）
Y_BOTTOM_NEED <- max((seq_along(ROWS) - 1L) * STEP + H_BOX)
H <- MARGIN + Y_BOTTOM_NEED + 0.13 + 0.06 + 0.10 + CAP_H + MARGIN
cat(sprintf("  版面：%.2f x %.2f in（%.0f x %.0f mm）；箱宽 %.3f in（%.1f mm）\n",
            WT, H, WT * 25.4, H * 25.4, W_BOX, W_BOX * 25.4))

# ── 5. 几何：箱子、行、箭头、括注、图注
xL <- MARGIN + (seq_len(NCOL) - 1L) * (W_BOX + GAP)
xR <- xL + W_BOX
yT <- H - MARGIN - (seq_len(NCOL) - 1L) * STEP
yB <- yT - H_BOX

rects <- data.frame(xmin = xL, xmax = xR, ymin = yB, ymax = yT,
                    fill = c("grey94", rep("white", 4)), layer = seq_len(NCOL))

txts <- do.call(rbind, lapply(seq_len(NCOL), function(j) {
  # 从箱顶**再退一个 PAD** 起排：H_BOX 里算了上下各一个 PAD，若这里从 yT 直接起排，
  # 顶部内边距就是 0（首行贴框线），底部却空出两个 PAD。本轮实测即如此，已修。
  rows <- ROWS[[j]]; y <- yT[j] - PAD; out <- list()
  for (r in rows) {
    h <- switch(r$k, head = LINE_H * 1.25, gap = GAP_H, LINE_H)
    yc <- y - h / 2; y <- y - h
    if (r$k == "head")
      out[[length(out) + 1L]] <- data.frame(x = xL[j] + PAD, y = yc, label = r$t, hjust = 0,
                                            fontface = "bold", sz = SZ_HEAD, layer = j)
    else if (r$k == "full")
      out[[length(out) + 1L]] <- data.frame(x = xL[j] + PAD, y = yc, label = r$t, hjust = 0,
                                            fontface = "plain", sz = SZ, layer = j)
    else if (r$k == "pair")
      out[[length(out) + 1L]] <- rbind(
        data.frame(x = xL[j] + PAD,     y = yc, label = r$l, hjust = 0,
                   fontface = "plain", sz = SZ, layer = j),
        data.frame(x = xR[j] - PAD,     y = yc, label = r$v, hjust = 1,
                   fontface = "plain", sz = SZ, layer = j))
  }
  do.call(rbind, out)
}))

# 层间箭头：第 j 箱右缘中点 → 第 j+1 箱左缘中点
segs <- data.frame(x = xR[-NCOL] + 0.02, xend = xL[-1] - 0.02,
                   y = (yT[-NCOL] + yB[-NCOL]) / 2, yend = (yT[-1] + yB[-1]) / 2)
must(all(segs$xend > segs$x), sprintf("四条层间箭头的落点都在起点右侧（缝宽 %.2f in）", GAP))

y_low  <- min(yB)
y_brk  <- y_low - 0.13
NOTE   <- "Supplementary Materials: machine learning (Table S2)"
NOTE_Y <- y_low - 0.06

# ── 6. 断言：逐行实测字宽装得进箱体（箱宽是版面约束，不是凑文字的自变量）
fit <- do.call(rbind, lapply(seq_len(NCOL), function(j) {
  rs <- ROWS[[j]]
  # 标题行按 SZ_HEAD 量、正文行按 SZ 量——用错磅值量会低估标题宽度
  rs <- lapply(rs, function(r) {
    if (r$k == "head") return(list(txt = r$t, sz = SZ_HEAD))
    if (r$k == "gap")  return(NULL)
    if (r$k == "pair") return(list(txt = c(r$l, r$v), sz = SZ))
    list(txt = r$t, sz = SZ)
  })
  rs <- Filter(Negate(is.null), rs)
  data.frame(layer = j,
             text = unlist(lapply(rs, `[[`, "txt")),
             w = unlist(Map(function(z) vapply(z$txt, tw, 0, sz = z$sz), rs)),
             stringsAsFactors = FALSE)
}))
# 纵向：行距须容得下该磅值的行盒，内边距须容得下半个行盒——否则相邻行会挤在一起、
# 或最后一行压住框线。这两条是**版面常数之间的关系**，与具体文字无关，故独立于字宽检查。
must(LINE_H >= SZ * 1.15 / 72,
     sprintf("正文行距 %.3f in ≥ 行盒 %.3f in（%.1f pt x 1.15）", LINE_H, SZ * 1.15 / 72, SZ))
must(LINE_H * 1.25 >= SZ_HEAD * 1.15 / 72,
     sprintf("标题行距 %.3f in ≥ 行盒 %.3f in（%.1f pt x 1.15）", LINE_H * 1.25, SZ_HEAD * 1.15 / 72, SZ_HEAD))
must(PAD >= 0.5 * SZ_HEAD * 1.2 / 72,
     sprintf("内边距 %.3f in ≥ 半行盒 %.3f in（%.1f pt x 1.2 / 2）", PAD, 0.5 * SZ_HEAD * 1.2 / 72, SZ_HEAD))
fit$limit <- W_BOX - 2 * PAD
must(all(fit$w <= fit$limit),
     sprintf("箱内 %d 行文字全部装得下（最宽 %.3f in / 上限 %.3f in）", nrow(fit), max(fit$w), fit$limit[1]))
must(max(fit$w) / fit$limit[1] > 0.80,
     sprintf("最宽行占上限的 %.0f%%——上限是紧的，不是随手放大后必然通过的空断言",
             100 * max(fit$w) / fit$limit[1]))

# 负面对照：箱宽减半后必须**有人放不下**，否则说明这条检查没有牙齿
neg <- fit$w <= (W_BOX / 2 - 2 * PAD)
must(!all(neg), sprintf("负面对照：箱宽减半后 %d/%d 行放不下（证明上面的检查会失败，不是恒真）",
                        sum(!neg), length(neg)))

# wrap_measured 是「超宽才断行」，故这条检查唯一能失败的路径是**单个词本身就比画布宽**
# （那时 `|| !nzchar(cur)` 会把它整词吐出来）——这正是它守的情形，不是空断言。
cap_w <- vapply(cap_lines, tw, 0, sz = SZ_CAP)
must(all(cap_w <= CAP_W),
     sprintf("图注 %d 行全部不超画布宽（最宽 %.3f / %.3f in）；无单词宽于画布",
             length(cap_lines), max(cap_w), CAP_W))

# 上图的每一个字符串都必须是 ASCII（基线 pdf 设备非 ASCII 会静默丢字）
allstr <- c(txts$label, NOTE, cap_lines)
bad <- allstr[grepl("[^ -~]", allstr)]
must(length(bad) == 0L, sprintf("%d 个待绘制字符串全部 ASCII", length(allstr)))
must(abs(STEP - (yT[1] - yT[2])) < 1e-12 && all(diff(yT) < 0),
     sprintf("五级台阶逐级下沉 %.3f in，序严格递减", STEP))

# ── 6b. E2 追加：交付件与终稿的一致性守卫（四处修正各有一条对应断言）
must(file.exists(MS), sprintf("图注来源为终稿：%s", basename(MS)))
must(!grepl("functional enrichment", leg, fixed = TRUE),
     "图注不含 functional enrichment（该模块未交付任何补充材料）")
must(grepl("Machine learning is reported in the Supplementary Materials", leg, fixed = TRUE),
     "图注采用终稿口径 Machine learning is reported ...")
flatten_rows <- function(rs) unlist(lapply(rs, function(r) {
  if (r$k == "gap") return(character(0))
  if (r$k == "pair") return(c(r$l, r$v))
  r$t
}))
box_txt <- unlist(lapply(ROWS, flatten_rows))
must(!any(grepl("first time", box_txt, fixed = TRUE)),
     sprintf("箱内 %d 段文字无 for the first time 声明（工作令 §四）", length(box_txt)))
must(!any(grepl("[^ -~]", box_txt)), "箱内文字全部 ASCII")
must("5 cell types)" %in% box_txt,
     "箱内细胞类型为五类（含 mast cells，与 §2.5 及图 5 一致）")
must(!any(grepl("4 cell types", box_txt, fixed = TRUE)), "箱内不残留四类口径")
must(grepl("machine learning", NOTE, ignore.case = TRUE) && !grepl("enrichment", NOTE),
     sprintf("括注只承诺机器学习：%s", NOTE))
must(!grepl("enrichment", NOTE, ignore.case = TRUE), "括注不含 enrichment")

# ── 7. 出图
# 三个绘制层分开写，避免把长度>1 的向量当常数传给非 aes 参数（ggplot 会静默取错长度）
head_txt <- txts[txts$fontface == "bold", ]
body_txt <- txts[txts$fontface != "bold", ]
p <- ggplot() +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = rects$fill, colour = "grey35", linewidth = 0.4) +
  geom_segment(data = segs, aes(x = x, xend = xend, y = y, yend = yend),
               colour = "grey30", linewidth = 0.35,
               arrow = arrow(length = unit(0.055, "in"), type = "closed")) +
  geom_text(data = body_txt, aes(x = x, y = y, label = label, hjust = hjust),
            size = SZ * PT2MM, fontface = "plain", colour = "grey15", family = "") +
  geom_text(data = head_txt, aes(x = x, y = y, label = label, hjust = hjust),
            size = SZ_HEAD * PT2MM, fontface = "bold", colour = "grey10", family = "") +
  annotate("segment", x = xL[1], xend = xR[NCOL], y = y_brk, yend = y_brk,
           colour = "grey45", linewidth = 0.35) +
  annotate("segment", x = xL[1], xend = xL[1], y = y_brk, yend = y_brk + 0.05,
           colour = "grey45", linewidth = 0.35) +
  annotate("segment", x = xR[NCOL], xend = xR[NCOL], y = y_brk, yend = y_brk + 0.05,
           colour = "grey45", linewidth = 0.35) +
  annotate("text", x = (xL[1] + xR[NCOL]) / 2, y = NOTE_Y, label = NOTE,
           size = SZ * PT2MM, colour = "grey30", hjust = 0.5) +
  # 图注题名单独一层：稿件图注里它是加粗的（**Figure 1: …**），图内也照加粗，避免两处观感不一
  annotate("text", x = MARGIN, y = y_brk - 0.12, label = TITLE,
           size = SZ_CAP * PT2MM, fontface = "bold", colour = "grey20", hjust = 0, vjust = 1) +
  annotate("text", x = MARGIN, y = y_brk - 0.12 - SZ_CAP * 1.15 / 72, label = paste(BODY_LINES, collapse = "\n"),
           size = SZ_CAP * PT2MM, colour = "grey25", hjust = 0, vjust = 1, lineheight = 1.15) +
  coord_fixed(ratio = 1, xlim = c(0, WT), ylim = c(0, H), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

PDF <- file.path(FIG, "Fig_study_design.pdf")
PNG <- file.path(FIG, "Fig_study_design.png")
ggsave(PDF, p, width = WT, height = H, device = pdf)
ggsave(PNG, p, width = WT, height = H, dpi = 400)
must(file.exists(PDF) && file.exists(PNG), "PDF 与 PNG 均已写出")
raw <- readBin(PDF, "raw", 8)
must(identical(rawToChar(raw[1:5]), "%PDF-"), "PDF 魔数 %PDF- 正确")

# 交付件名：E3 从 04_图表/ 取冻结件、以 Figure1_study_design.* 之名发出。两处必须是同一份
# 字节 —— 各存一次必然漂移（本项目已有「共享文字被两个生成器各存一份」的先例）。
# 故先写到冻结位，再 file.copy 到交付位，并逐字节核对；不做第二次 ggsave。
DELIV_PDF <- file.path(DELIV, "Figure1_study_design.pdf")
DELIV_PNG <- file.path(DELIV, "Figure1_study_design.png")
for (pr in list(c(PDF, DELIV_PDF), c(PNG, DELIV_PNG))) {
  if (normalizePath(pr[1]) != normalizePath(pr[2], mustWork = FALSE)) {
    ok <- file.copy(pr[1], pr[2], overwrite = TRUE)
    must(ok, sprintf("已复制到交付件名 %s", basename(pr[2])))
    must(identical(tools::md5sum(pr[1])[[1]], tools::md5sum(pr[2])[[1]]),
         sprintf("%s 与冻结件 %s 逐字节一致", basename(pr[2]), basename(pr[1])))
  }
}

# ── 8. 审计表：把版面与实测字宽落盘，事后可核（不止 print 到 stdout）
fit$limit_used <- fit$limit
write.csv(fit[, c("layer", "text", "w", "limit_used")],
          file.path(TAB, "图1_研究设计_布局.csv"), row.names = FALSE)
lay <- data.frame(layer = seq_len(NCOL),
                  x_left = round(xL, 4), x_right = round(xR, 4),
                  y_top = round(yT, 4), y_bottom = round(yB, 4), height = round(H_BOX, 4),
                  head = vapply(ROWS, function(rs) rs[[1]]$t, ""),
                  # 每箱的文本片段数（pair 行出两段）。D4b 拿它去数 PDF 文本层，两边必须相等：
                  # 若某行被裁掉或漏画，片段数会对不上。
                  n_text = as.integer(table(factor(txts$layer, levels = seq_len(NCOL)))))
write.csv(lay, file.path(TAB, "图1_研究设计_箱体.csv"), row.names = FALSE)
# 图内实际绘制的图注与括注逐字落盘：D4b 拿 PDF 文本层与它比对，
# 免得把 ASCII 映射表在 R 与 Python 两侧各写一份（两份必然漂移）。
writeLines(c(TITLE, BODY_LINES), file.path(TAB, "图1_研究设计_图注.txt"), useBytes = FALSE)
writeLines(enc2utf8(NOTE), file.path(TAB, "图1_研究设计_括注.txt"), useBytes = FALSE)
cat(sprintf("\n写出 %s（%.1f KB）\n     %s（%.1f KB）\n",
            basename(PDF), file.size(PDF) / 1024, basename(PNG), file.size(PNG) / 1024))
cat(sprintf("审计表 图1_研究设计_布局.csv（%d 行）/ 图1_研究设计_箱体.csv（%d 行）\n", nrow(fit), nrow(lay)))
cat(sprintf("  画布 %.2f x %.2f in；箱宽 %.3f in；台阶下沉 %.2f in；图注 %d 行\n",
            WT, H, W_BOX, STEP, length(cap_lines)))

dev.off()
cat(sprintf("\n结果：%s\n", if (FAIL == 0L) "全部通过" else sprintf("%d 项未通过", FAIL)))
if (FAIL > 0L) quit(status = 1)
