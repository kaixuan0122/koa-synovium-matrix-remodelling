# B7 森林图（§4.2）：逐队列效应 + 合并估计，发现层与复制层分面
#
# 两条纪律：
#  ① 区间必须与正文同一检验。两行的区间都直接取 B2 的表列，不在这里另算 ——
#     主估计行取 primary_CI_low/primary_CI_high（即 hk_*，REML+HK：metafor 在
#     test="knha" 下用 t(k−1) 分位，且 SE 本身也换成了 HK 的方差估计量）；
#     敏感性行按层取（发现层 ivw_*：IVW-FE 的正态近似；复制层 pooled_*：
#     ~0+group+cohort + HC3 稳健 SE 的 t 分布）。三者口径不同，
#     若在此用同一个公式重算，图与正文就会在边界上打架。
#  ② 发现层的逐队列点与 **IVW-FE 敏感性行**是**同一个估计量**（都用各队列自身的估计：GSE89408 为
#     DESeq2 系数、其余为 limma），故 IVW 合并是这些点的加权组合，图上自洽。
#     复制层三个队列同为微阵列，limma 的 logFC 与 ~0+group+cohort 的均值差同一个量。
#     发现层另有一个 pooled(~0+group+cohort) 估计（表内对照列），它估的是「对数均值差」，
#     对 GSE89408 与 DESeq2 系数是**不同的量**（LOXL1 可差到符号相反），故不上图。
#
#  ③ ★ 2026-10-01 第二轮改换口径后（见文件末 ★ 注）：两层的主估计都是 **REML+HK**
#     （`method="REML", test="knha"`，§2.3 声明），画成红色菱形置于各层队列行之上；
#     蓝色三角是**敏感性**行 —— 发现层 IVW-FE
#     （与主估计同为逐队列估计的加权组合，二者只差一个组间方差分量，是模型比较）、
#     复制层 pooled(~0+group+cohort)（作用于样本级取值，与主估计既换估计量又换方差处理，
#     是估计策略比较，§2.3(iv)）。两层敏感性行不是同一个估计量，故标签按层给；
#     主估计行两层同名 "REML+HK (main)"。
#
#  ④ **预测区间不上图（有意为之，非遗漏）。** 审稿意见要的是「报告预测区间」，本稿把 PI
#     落在 B2 的冻结汇总表（pi_low/pi_high）与正文里。若叠到本图上会毁掉图：PI 比 HK 的
#     CI 宽 5～10 倍（SDC4 发现层 PI=(-12.93,+11.41) 而 CI=(-2.90,+1.37)），而
#     facet_grid 只让 y 自由、**x 轴跨全部基因共用**，故 x 量程由最宽的 PI 决定，
#     全部队列点与 CI 会被压成中间一条竖线。若将来要画 PI，须先改成 x 也自由（各自成图）
#     或另开一张专画 PI 的图，不能直接叠加。
#
# 输入：02_结果表/逐队列效应_力学基因.csv（B1）、合并效应_{发现,复制}层.csv（B2）
# 输出：04_图表/Fig_forest_logFC.pdf / .png（冻结件位，E3 从这里取图）
#       03_图表与源数据/Figure3_forest_logFC.pdf / .png（交付件名，逐字节复制自前者）
#       02_结果表/森林图_合并行数值.csv

# ★ 2026-09-29 两处修正（回应 R3-M2）：
#   ① **行标签必须自陈估计量。** 原图把两层的合并行都标成 "Pooled"，而发现层的主估计是
#      IVW-FE —— 读者会以为两层同一估计量，正是 §2.3(iv) 反复声明的那个混淆。现分层标：
#      发现层 "IVW-FE (main)"、复制层 "Pooled (cohort-adjusted)"，与主稿 Figure 3 图注逐字对应。
#      落盘的 森林图_合并行数值.csv 同步把 kind 列换成 estimator 列（同图上一套措辞），
#      并加断言：CSV 里的措辞集合必须与图上标签集合相等。
#   ② **路径原本指向旧基线树**（`…/滑膜基质刚度生成与刚度感知 0912/08_修订_2026-09-17`）。
#      那样跑出来的图落在旧树里，交付目录里的 图 3 不会被更新 —— 上游改了而交付件照旧。
#      现改为从本脚本位置推树根，与 E2/B3b 同款。
#
# ★ 2026-10-01 第十轮（回应审稿 P0-1「两层各用不同估计量、pooled 居复制层主位」）：
#   ① **主估计改为 RE-DL，两层同一。** 原红菱形取 primary_*（发现层 IVW-FE、复制层
#      pooled），现改取 re_* 列；IVW-FE 与 pooled 双双降为敏感性（蓝三角）。
#      数据映射不变 —— df 里那两行本来就分别装着 primary_* 与 re_*，本轮只把
#      形状/颜色/大小的归属与 kind 名对调，使图与 §2.3 的新声明一致。
#   ② kind 由 Pooled/RE 改名为 Sens/Prim —— 旧名会让发现层的 IVW-FE 行在代码里自称
#      "Pooled"，正是本文件 ① 记的那个混淆。落盘的 CSV 仍只出 estimator 列。
#   ③ 断言同步：`kind == "Pooled"` → `"Sens"`、`"RE"` → `"Prim"`；行标签常量
#      POOL_LAB/RE_LAB → SENS_LAB/PRIM_LAB。E3 的 [L] 段按新常量重取。
#   ④ 主稿 Figure 3 图注（h.md:308）本轮已先行改定，是本次改图**依据的规格**，
#      不是事后追认：红菱形 = 两层主估计 RE-DL，蓝三角 = 各层敏感性估计。
#
# ★ 2026-10-01 第十一轮（第二轮审稿 P0-1 改判：主估计必须是 REML + Hartung–Knapp）：
#   ① **主估计由 RE-DL 再改为 REML+HK，两层同一。** 红菱形取 primary_* —— B2 已把
#      primary_* 指向 hk_*（HK 换的是**方差估计量**，不只是分位数：se_HK² = Σwᵢ(yᵢ−β̂)²/((k−1)Σwᵢ)，
#      wᵢ=1/(vᵢ+τ²)，再用 t(k−1) 分位），故此处**不再需要 re_* 列**。
#   ② **敏感性行改为按层取列**：发现层 ivw_*、复制层 pooled_*。上一版两层的敏感性行都取
#      primary_*（当时发现层 = IVW-FE、复制层 = pooled，恰好对），primary_* 换成 hk 后
#      若照旧取 primary_*，两个敏感性行会双双变成主估计的复制品、图上少一个估计量。
#      行标签 SENS_LAB 不变（两层敏感性的估计量本来就没变）。
#   ③ 图注现写「REML+HK」。E3 的 [L] 段 EST_DESC 表与 L6c 的 "DerSimonian" 判据须同步。
suppressPackageStartupMessages({library(ggplot2)})

ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "..", ".."))
RES  <- file.path(ROOT, "02_结果表")
OUT  <- file.path(ROOT, "04_图表")                    # 冻结件位
DELIV<- file.path(ROOT, "..", "03_图表与源数据")        # 交付件位
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PANEL <- c("LOXL1", "VGLL4", "SDC4", "SERPINE1", "ITGA5", "SAV1")
DISCOVERY <- c("GSE89408", "GSE206848", "GSE1919")
REPLICATION <- c("GSE55235", "GSE55457", "GSE82107")
# 行标签。**只能用 ASCII**：基础 pdf 设备走 mbcsToSbcs 转 PostScript 编码，
# †(U+2020) / ‡(U+2021) 这类字符会触发 "conversion failure" 警告并**静默丢字**，
# 而 ggsave 仍返回成功、PNG 设备却画得出来 —— 只在 PDF 里坏，肉眼查 PNG 查不出来。
# 层间区分靠 facet 条（Discovery/Replication）+ 图注，不靠符号。
# ★ 2026-10-01：主估计两层同名；敏感性行的标签按层给（两层敏感性的估计量不同，见文件头 ③）。
# 行标签只写估计量名，角色（主/敏感性）由图注与文件末 ★ 注交代。
SENS_LAB <- c(Discovery = "IVW-FE", Replication = "Pooled (cohort-adjusted)")
PRIM_LAB <- c(Discovery = "REML+HK (main)", Replication = "REML+HK (main)")
# 敏感性行的**取列**规则：发现层 ivw_*、复制层 pooled_*。与上面的标签是两个变量 ——
# 标签是手写的字符串、取列是这里定的映射，二者错配后图上仍"看着合理"
# （见笔记 summary-table-label-must-bind-to-taken-coefficient）。故下面逐列断言取对了源列。
SENS_COL <- c(Discovery = "ivw", Replication = "pooled")
# 两层主估计是**同一个估计量**（§2.3 声明），故标签同名。但 factor(levels=) 不接受
# 重复水平 —— R >= 4.0 直接报 "factor level [k] is duplicated"，不是警告。两层主估计行
# 各在本人面版占一行、且要落在同一高度，共用一个水平正是想要的画法，故在此去重。
LV <- unique(c(unname(SENS_LAB), DISCOVERY, REPLICATION, unname(PRIM_LAB)))
stopifnot(!anyDuplicated(LV))

per <- read.csv(file.path(RES, "逐队列效应_力学基因.csv"), stringsAsFactors = FALSE)
per <- per[per$gene %in% PANEL, ]
per$layer <- ifelse(per$dataset %in% DISCOVERY, "Discovery", "Replication")
stopifnot(setequal(unique(per$layer), c("Discovery", "Replication")))

pf <- read.csv(file.path(RES, "合并效应_发现层.csv"), stringsAsFactors = FALSE)
pr <- read.csv(file.path(RES, "合并效应_复制层.csv"), stringsAsFactors = FALSE)
pool <- rbind(pf, pr)
pool$layer <- ifelse(pool$layer == "发现", "Discovery", "Replication")
pool <- pool[pool$gene %in% PANEL, ]
# 合并行的区间必须带出它用的方法，图注据此说明。
# 本轮起两层的主估计量同名（REML+HK），故这里是一个值而不是按层分支。
stopifnot(all(pool$primary_method == "REML+HK"))
# primary_* 必须确实取自 hk_*：挂了新名字却仍读旧列，是这类改口径最容易出的错。
stopifnot(all(pool$primary_logFC == pool$hk_logFC),
          all(pool$primary_P == pool$hk_P))
cat(sprintf("主估计量（两层）：%s\n", unique(pool$primary_method)))
cat(sprintf("敏感性行：发现层 %s（取 %s_* 列）/ 复制层 %s（取 %s_* 列）\n",
            SENS_LAB[["Discovery"]], SENS_COL[["Discovery"]],
            SENS_LAB[["Replication"]], SENS_COL[["Replication"]]))

# 队列水平的 factor 取 rev(LV)：ggplot 把第一个水平画在最下，故 rev 后 RE 行在最上、
# 合并行在其下、队列行在最下。不再用 scale_y_discrete(limits=) 钉住全部水平 —— 那会把
# 不属于本层的队列也留成空行，读者会当成这些队列没数据。改为 free_y + drop=TRUE。
LV_REV <- rev(LV)

# 敏感性行按层取列（见 SENS_COL）。先落到 pool 上再进 df —— 这样下面 chk 的回比可以
# 直接 merge 同一批列，不必把取列的规则写第二遍（写两遍必然漂移）。
sens_pick <- function(stem) {
  v <- rep(NA_real_, nrow(pool))
  for (L in names(SENS_COL)) {
    i <- pool$layer == L
    v[i] <- pool[[paste0(SENS_COL[[L]], "_", stem)]][i]
  }
  v
}
pool$sens_logFC    <- sens_pick("logFC")
pool$sens_CI_low   <- sens_pick("CI_low")
pool$sens_CI_high  <- sens_pick("CI_high")
pool$sens_P        <- sens_pick("P")
# 取列规则的自证：发现层那一半必须逐位等于 ivw_*，复制层那一半必须逐位等于 pooled_*。
# 若有人把 SENS_COL 写反或写重，上面两个和式的哪一个都会立刻不成立。
for (L in names(SENS_COL)) {
  i <- pool$layer == L
  stopifnot(all(pool$sens_logFC[i] == pool[[paste0(SENS_COL[[L]], "_logFC")]][i]))
}
stopifnot(all(is.finite(pool$sens_logFC)), all(is.finite(pool$sens_P)))

df <- rbind(
  data.frame(gene = per$gene, layer = per$layer, est = per$logFC, lo = per$CI_low,
             hi = per$CI_high, ylab = factor(per$dataset, levels = LV_REV),
             kind = "Cohort", n = NA_real_, I2 = NA_real_, P = per$P,
             stringsAsFactors = FALSE),
  data.frame(gene = pool$gene, layer = pool$layer, est = pool$sens_logFC,
             lo = pool$sens_CI_low, hi = pool$sens_CI_high,
             ylab = factor(unname(SENS_LAB[pool$layer]), levels = LV_REV), kind = "Sens",
             n = pool$k, I2 = pool$I2, P = pool$sens_P, stringsAsFactors = FALSE)
)

# 主估计行（REML+HK）只在 B2 真的算出了 primary_logFC 时再加。用它自己的区间列，不在这里重算 ——
# 与 ① 同一条纪律：图上的区间必须就是表里那一列。
n_prim <- sum(is.finite(pool$primary_logFC))
if (n_prim) {
  ok <- is.finite(pool$primary_logFC)
  df <- rbind(df, data.frame(
    gene = pool$gene[ok], layer = pool$layer[ok],
    est = pool$primary_logFC[ok], lo = pool$primary_CI_low[ok], hi = pool$primary_CI_high[ok],
    ylab = factor(unname(PRIM_LAB[pool$layer[ok]]), levels = LV_REV),
    kind = "Prim", n = pool$k[ok], I2 = pool$I2[ok], P = pool$primary_P[ok],
    stringsAsFactors = FALSE))
}
cat(sprintf("叠加主估计行（REML+HK）：%d（发现层 %d / 复制层 %d）\n", n_prim,
            sum(is.finite(pool$primary_logFC) & pool$layer == "Discovery"),
            sum(is.finite(pool$primary_logFC) & pool$layer == "Replication")))

# 断言各层只出现本层队列、合并行与 RE 行：上面那个空行缺陷就是从这一步漏掉的
for (L in c("Discovery", "Replication")) {
  got <- unique(as.character(df$ylab[df$layer == L]))
  want <- c(unname(SENS_LAB[[L]]), if (L == "Discovery") DISCOVERY else REPLICATION, PRIM_LAB[[L]])
  if (!setequal(got, want)) stop("层 ", L, " 的水平不符：", paste(got, collapse = ","))
}
# 结构断言：合并行在每层每基因恰好一行；图中不应出现无数据的层×基因组合
stopifnot(nrow(df[df$kind == "Sens", ]) == 12L)
stopifnot(all(is.finite(df$est)), all(df$lo < df$hi))
# 主估计行若出现，必须两层每基因都有 —— 只画一层会被读成「另一层没有异质性」
if (n_prim) stopifnot(nrow(df[df$kind == "Prim", ]) == 12L)
# 敏感性行与主估计行的数值必须逐格来自 B2 的表：把 df 里的值倒回去比对，
# 防止 cbind/rbind 列错位后图上出现「看着合理但取错列」的点。
chk <- merge(df[df$kind == "Sens", c("layer", "gene", "est")],
             pool[, c("layer", "gene", "sens_logFC")], by = c("layer", "gene"))
stopifnot(all(abs(chk$est - chk$sens_logFC) < 1e-12))
if (n_prim) {
  chk2 <- merge(df[df$kind == "Prim", c("layer", "gene", "est")],
                pool[, c("layer", "gene", "primary_logFC")], by = c("layer", "gene"))
  stopifnot(all(abs(chk2$est - chk2$primary_logFC) < 1e-12))
}

df$gene <- factor(df$gene, levels = PANEL)
df$kind <- factor(df$kind, levels = c("Cohort", "Sens", "Prim"))

# 硬拦非 ASCII 标签。上面那次 †/‡ 的教训：PNG 看着正常、PDF 静默丢字、退出码 0，
# 事后只能靠肉眼比对两种设备才发现。这里让它直接失败。
# 判据是「可打印 ASCII 之外还有没有字符」，不写成逐个字符白名单（后者会漏掉全角括号等）。
bad_lab <- grep("^[ -~]*$", unique(as.character(df$ylab)), value = TRUE, invert = TRUE)
if (length(bad_lab)) stop("y 轴标签含非 ASCII 字符，基础 pdf 设备会静默丢字：",
                          paste(bad_lab, collapse = " | "))
# 主估计行的 k 必须取自同一批逐队列数据（B2 里 k 与 k_meta 已验证相同，这里再钉一次，
# 防止将来 B2 改了 meta_ivw 的入组规则而 k 与 k_meta 分叉后图上印出错误的队列数）
if (n_prim) stopifnot(all(abs(pool$k[is.finite(pool$primary_logFC)] -
                              pool$k_meta[is.finite(pool$primary_logFC)]) < 1e-12))

# 图注要写出每个敏感性行/主估计行的 k 与 I²，故把它们落盘备查。
# 第一列写**估计量**（与图上标签同一套措辞），不写内部的 kind 名：kind 只是画图用的
# 因子（决定形状/颜色），把它当标签交出去会让发现层的 IVW-FE 行在源数据里自称 "Sens"。
out <- df[df$kind != "Cohort", c("kind", "layer", "gene", "est", "lo", "hi", "n", "I2", "P")]
out$estimator <- ifelse(out$kind == "Sens", unname(SENS_LAB[out$layer]),
                        unname(PRIM_LAB[out$layer]))
# 不变量：源数据里的估计量措辞集合 == 图上 y 轴标签集合。两处措辞若各改各的，这里就炸。
stopifnot(setequal(unique(out$estimator), setdiff(unique(as.character(df$ylab)), c(DISCOVERY, REPLICATION))))
write.csv(out[, c("estimator", "layer", "gene", "est", "lo", "hi", "n", "I2", "P")],
          file.path(RES, "森林图_合并行数值.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 面板基因两层「敏感性估计 vs 主估计（REML+HK）」的方向/显著性是否一致，落一份文字记录供正文引用
cat("\n面板基因 敏感性 vs 主估计（REML+HK）：\n")
cat(sprintf("  %-9s %-10s %10s %9s %10s %9s %s\n",
            "gene", "层", "敏感性", "敏感P", "主估计", "主P", "状态"))
for (L in c("Discovery", "Replication")) for (g in PANEL) {
  r <- df[df$gene == g & df$layer == L & df$kind != "Cohort", ]
  if (nrow(r) < 2) { cat(sprintf("  %-9s %-10s  （主估计不可用）\n", g, L)); next }
  p <- r[r$kind == "Sens", ]; e <- r[r$kind == "Prim", ]
  st <- if (sign(p$est) != sign(e$est)) "方向改变"
        else if ((p$P < 0.05) != (e$P < 0.05)) "显著性改变" else "不变"
  cat(sprintf("  %-9s %-10s %+10.3f %9.3g %+10.3f %9.3g %s\n",
              g, L, p$est, p$P, e$est, e$P, st))
}

p <- ggplot(df, aes(x = est, y = ylab, colour = kind, shape = kind, size = kind)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.18, linewidth = 0.45) +
  geom_point() +
  facet_grid(layer ~ gene, scales = "free_y", space = "free_y", switch = "y") +
  scale_shape_manual(values = c(Cohort = 16, Sens = 17, Prim = 18), guide = "none") +
  scale_colour_manual(values = c(Cohort = "grey30", Sens = "#2166AC", Prim = "#B2182B"),
                      guide = "none") +
  scale_size_manual(values = c(Cohort = 1.7, Sens = 2.3, Prim = 2.9), guide = "none") +
  scale_y_discrete(drop = TRUE) +
  labs(x = expression(log[2]~"fold change (OA vs Normal)"), y = NULL) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = "grey93", colour = "grey70"),
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 7.5),
        panel.spacing.x = unit(4, "pt"),
        panel.spacing.y = unit(6, "pt"))

# 本机 R 4.6 的 cairo 缺 X11（/opt/X11/lib/libXrender.1.dylib 不存在），cairo_pdf 会报
# "failed to load cairo DLL" 并**不产出文件**。所有标签均为拉丁字母与表达式，无 CJK，
# 故直接用基础 pdf 设备；若将来图上出现中文，须改 quartz 设备（见笔记 r-plot-cjk-needs-quartz-device）。
PDF <- file.path(OUT, "Fig_forest_logFC.pdf"); PNG <- file.path(OUT, "Fig_forest_logFC.png")
ggsave(PDF, p, width = 11.5, height = 5.2, device = pdf)
ggsave(PNG, p, width = 11.5, height = 5.2, dpi = 400)

# 交付件名：E3 从 04_图表/ 取冻结件、以 Figure3_forest_logFC.* 之名发出。先写冻结位、
# 再复制成交付件名并逐字节核对（各存一次必然漂移；见 E2/B3b 同款做法）。
for (pr in list(c(PDF, file.path(DELIV, "Figure3_forest_logFC.pdf")),
                c(PNG, file.path(DELIV, "Figure3_forest_logFC.png")))) {
  if (normalizePath(pr[1]) != normalizePath(pr[2], mustWork = FALSE)) {
    if (!file.copy(pr[1], pr[2], overwrite = TRUE)) stop("复制到交付件名失败：", pr[2])
    if (tools::md5sum(pr[1])[[1]] != tools::md5sum(pr[2])[[1]])
      stop("交付件与冻结件字节不同：", pr[2])
  }
}

# 断言而不是宣称：ggsave 在设备加载失败时会报错却可能已把错误吞掉。
# 判据用 PDF 的结构标记（页面对象数），不用文件大小 —— 大小阈值是拍的，矢量图 10 KB 上下都正常。
# 必须走原始字节：readLines/rawToChar 读二进制会在转宽字符时报 invalid，检查器自己先坏掉。
for (f in c(PDF, PNG)) if (!file.exists(f)) stop("未产出：", f)
raw <- readBin(PDF, "raw", file.size(PDF))
if (!identical(rawToChar(raw[1:5]), "%PDF-")) stop("不是 PDF：", PDF)
# "/Type /Page" 同时命中页面树节点 "/Type /Pages"，故页面数 = 命中数 - 1
npage <- length(grepRaw("/Type /Page", raw, all = TRUE)) - 1L
if (npage < 1L) stop("PDF 无页面对象：", PDF)
cat(sprintf("\n已核验：%s（%d 页，%.1f KB）\n        %s（%.1f KB）\n",
            PDF, npage, file.size(PDF) / 1024, PNG, file.size(PNG) / 1024))
cat("\n合并行（图上菱形 ● 与三角 ▲）：\n")
print(out[, c("estimator", "layer", "gene", "est", "lo", "hi", "n", "I2", "P")],
      row.names = FALSE)
