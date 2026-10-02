# B2 跨队列合并效应、留一法敏感性与异质性（§4.2）
#
# 主分析（按既定选择）：以 `~0 + group + cohort` 在 log 尺度矩阵上直接拟合，
#   队列作为协变量吸收平台差异。使用 HC3 稳健标准误——三个队列的平台、尺度、
#   残差方差都不同（发现 F6：线性尺度跨度 6.8 倍），经典 OLS 的等方差假设不成立，
#   稳健 SE 才能给出可信的区间。
#
# 对照分析：逆方差固定效应 Meta（metafor），它对每个队列按其自身估计精度加权。
#   为什么两个都做：主分析把三个队列的原始观测堆在一起拟合，每个队列贡献的是
#   「log 值均值差」；而逐队列估计里 GSE89408 来自 DESeq2 的计数 GLM 系数。
#   两者对同一队列可以是不同的量（实测 LOXL1：均值差 −0.29 vs GLM +0.74）。
#   两个方法都报，分歧处正是需要在正文交代的地方，不能只报一个。
#
# 输入：02_结果表/逐队列效应_力学基因.csv（B1）
#       04_分析数据/processed/{主分析,验证层}_genelevel.rds
#       03_分析/outputs/GSE89408_基因水平_双版本.rds（B0，size factor 归一化版）
# 输出：02_结果表/合并效应_发现层.csv、合并效应_复制层.csv、留一法_合并效应.csv
suppressPackageStartupMessages({library(sandwich); library(lmtest); library(metafor)})

BASE <- "/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/01.材料/01.博士材料/10.毕业/01.小论文/04.立项/滑膜基质刚度生成与刚度感知 0912"
P    <- file.path(BASE, "04_分析数据/processed")
# 2026-10-01：输出落点改由环境变量 B2_OUT 控制。输入（P 与 03_分析/outputs）一律只读
# 旧树，输出默认仍回旧树以保证历史可复现；本轮修订要求产物落在修订工作树内，
# 故交付/跑链时必须显式传 B2_OUT，否则会静默写回旧树。
OUT  <- Sys.getenv("B2_OUT", unset = file.path(BASE, "08_修订_2026-09-17/02_结果表"))

LAYER <- list(
  发现 = c("GSE206848", "GSE1919", "GSE89408"),
  复制 = c("GSE55235", "GSE55457", "GSE82107")
)
PANEL <- c("LOXL1", "VGLL4", "SDC4", "SERPINE1", "ITGA5", "SAV1")
CROSSLINK <- c("LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4", "PLOD1", "PLOD2", "TGM2")

## ---------- 1. 组装 log 尺度矩阵 ----------
main  <- readRDS(file.path(P, "主分析_genelevel.rds"))
valid <- readRDS(file.path(P, "验证层_genelevel.rds"))
g89   <- readRDS(file.path(BASE, "08_修订_2026-09-17/03_分析/outputs/GSE89408_基因水平_双版本.rds"))

MATS <- list(
  GSE206848 = main$GSE206848$mat, GSE1919 = main$GSE1919$mat,
  GSE55235  = valid$GSE55235$mat, GSE55457 = valid$GSE55457$mat, GSE82107 = valid$GSE82107$mat,
  GSE89408  = g89$size          # B0 的 size factor 归一化版，非 03_preprocess 的 colSums 版
)
GRPS <- list(
  GSE206848 = main$GSE206848$grp, GSE1919 = main$GSE1919$grp,
  GSE55235  = valid$GSE55235$grp, GSE55457 = valid$GSE55457$grp, GSE82107 = valid$GSE82107$grp,
  GSE89408  = g89$grp
)
cat("各队列矩阵（log 尺度）：\n")
for (nm in names(MATS)) cat(sprintf("  %-10s %5d 基因 x %2d 样本\n", nm, nrow(MATS[[nm]]), ncol(MATS[[nm]])))

eff <- read.csv(file.path(OUT, "逐队列效应_力学基因.csv"), stringsAsFactors = FALSE)
GENES <- unique(eff$gene)

# 基因符号逐队列可能不同：CTGF/CCN2、CYR61/CCN1、MKL1/MRTFA、MRTFB/MKL2。样本级矩阵的
# 行名是该队列**实际命中**的那个符号（B1 已逐队列记录在 symbol_used），故取基因必须按
# symbol_used 逐队列解析，不能拿规范名去硬查 rownames —— 硬查会静默 next，两种后果：
#   ① 整基因消失：CTGF/CYR61/MKL1 在复制层（三队列全叫 CCN2/CCN1/MRTFA）、
#      MRTFB 在发现层（GSE89408 叫 MKL2，只剩 GSE206848 一个队列）；
#   ② 更隐蔽的一种：CTGF/CYR61/MKL1 在发现层漏掉 GSE206848，只用 2 个队列算，
#      而报出的 k=2 与「本来就只有 2 个队列」的基因（GSE1919 无探针那批）撞巧一致，
#      数值看着完全正常。
SYM <- setNames(eff$symbol_used, paste(eff$gene, eff$dataset))

# 不变量：逐队列表里记录的每个 (基因, 队列)，其 symbol_used 都必须能在该队列矩阵里找到。
# 这条断言正是为上面那个缺陷设的——符号解析失败会静默 next，是整份产物里唯一
# 「不报错、不产出行、也不改变任何已有数值」的失败模式，靠肉眼看结果表发现不了。
.bad <- character(0)
for (i in seq_len(nrow(eff))) {
  s <- unname(SYM[paste(eff$gene[i], eff$dataset[i])])
  if (is.na(s) || !s %in% rownames(MATS[[eff$dataset[i]]]))
    .bad <- c(.bad, sprintf("%s / %s (symbol_used=%s)", eff$gene[i], eff$dataset[i], s))
}
if (length(.bad)) stop("符号解析失败，这些 (基因, 队列) 会被静默丢弃：\n  ", paste(.bad, collapse = "\n  "))
cat(sprintf("符号解析：%d 个 (基因, 队列) 全部命中矩阵行名\n", nrow(eff)))

## ---------- 2. 合并拟合 + 留一法 + 异质性 ----------
fit_pooled <- function(gene, cohorts) {
  # 返回 list(logFC, SE, CI_low, CI_high, P, n, k, method)，取不到返回 NULL
  ys <- list(); cs <- list(); gs <- list()
  for (nm in cohorts) {
    m <- MATS[[nm]]
    sym <- unname(SYM[paste(gene, nm)])   # 用 [ 而非 [[：键不存在时 [[ 会直接报 subscript out of bounds
    if (is.na(sym) || !sym %in% rownames(m)) next
    v <- m[sym, ]
    if (all(is.na(v)) || sd(v, na.rm = TRUE) == 0) next   # 恒定值不提供信息
    ys[[nm]] <- as.numeric(v); cs[[nm]] <- rep(nm, length(v))
    gs[[nm]] <- ifelse(GRPS[[nm]] == "OA", "OA", "Normal")
  }
  if (length(ys) < 2) return(NULL)
  d <- data.frame(y = unlist(ys), cohort = factor(unlist(cs)),
                  group = factor(unlist(gs), levels = c("Normal", "OA")))
  # 该基因在某队列内只有一组时无法提供组间信息，剔除该队列
  ok <- tapply(d$group, d$cohort, function(z) length(unique(z))) == 2
  if (sum(ok) < 2) return(NULL)
  d <- d[d$cohort %in% names(ok)[ok], ]; d$cohort <- droplevels(d$cohort)
  # ~0 + group + cohort：group 取全水平，cohort 吸收平台差异。
  # 参数化陷阱：group 取全水平后，groupNormal / groupOA 是**参照队列**两格的均值，
  # 不是组间差。组间效应须取对照 groupOA − groupNormal；只取 groupOA 得到的是表达
  # 水平（如 7.28）而不是差值（如 +0.65），数值看着仍"合理"，不会被任何断言拦住。
  fit <- lm(y ~ 0 + group + cohort, data = d)
  nm <- names(coef(fit))
  if (!all(c("groupNormal", "groupOA") %in% nm)) return(NULL)
  L <- setNames(rep(0, length(nm)), nm); L[["groupOA"]] <- 1; L[["groupNormal"]] <- -1
  V <- sandwich::vcovHC(fit, type = "HC3")
  b <- sum(L * coef(fit)); se <- sqrt(as.numeric(t(L) %*% V %*% L))
  df <- df.residual(fit); tq <- qt(0.975, df)
  # 该模型实际估的量：各队列内 log 值的均值差。落盘以便与 B1 的逐队列估计对照。
  coh <- sapply(split(d, d$cohort), function(z)
    mean(z$y[z$group == "OA"]) - mean(z$y[z$group == "Normal"]))
  list(logFC = unname(b), SE = unname(se),
       CI_low = unname(b - tq * se), CI_high = unname(b + tq * se),
       P = unname(2 * pt(-abs(b / se), df)),
       n = nrow(d), k = nlevels(d$cohort),
       used = names(ok)[ok],          # 实际进入拟合的队列，供下面的覆盖率断言比对
       cohort_logmeandiff = paste(sprintf("%s=%+.3f", sub("GSE", "", names(coh)), coh), collapse = "; "))
}

meta_ivw <- function(gene, cohorts) {
  sub <- eff[eff$gene == gene & eff$dataset %in% cohorts, ]
  sub <- sub[is.finite(sub$SE) & sub$SE > 0, ]
  if (nrow(sub) < 2) return(NULL)
  r <- try(metafor::rma(yi = sub$logFC, sei = sub$SE, method = "FE"), silent = TRUE)
  if (inherits(r, "try-error")) return(NULL)
  q <- metafor::rma(yi = sub$logFC, sei = sub$SE, method = "DL")
  # 随机效应（DerSimonian–Laird）作为**敏感性**估计，与固定效应并列落盘。
  # 本列不替换任何现有数值，只回答「I² 这么高时换成 RE 会不会翻」。
  # DL 的权重 1/(vi+tau2) 全为正，故 RE 估计必落在逐队列取值范围内——与上面的
  # 固定效应同一不变式，可被下面那条越界断言一并覆盖。

  # 2026-10-01（第二轮审稿意见第 1 项）：主估计量升为 REML + Hartung–Knapp。
  # 两件事必须分开记，否则会把机制的来源写错：
  #   ① REML 只换 tau² 的估计量（DL 的矩估计 → REML）。在 τ²=0 的格子上两者数值恒等；
  #      实测两层 109 行里 REML 与 DL 的显著/不显著判定只在复制层 2 格上有别。
  #   ② HK 换的是**方差估计量本身**，不只是分位数：
  #        se_HK² = (1/(k-1)) · Σ w_i (y_i − beta_hat)² / Σ w_i,   w_i = 1/(v_i + tau²)
  #      常规的 1/sqrt(Σw) 被换掉。队列彼此吻合时 se_HK 远**小于** 1/sqrt(Σw)
  #      （LOXL1 发现层 0.0387 vs 0.1800），队列分歧时远**大于**（VGLL4 复制层
  #      0.4579 vs 0.1359）。故 **HK 不是一致更保守**，正文必须按实测写。
  #      分位数同时由正态换成 t(df = k−1)：k=3 时 df=2、乘数 4.3027。
  h <- try(metafor::rma(yi = sub$logFC, sei = sub$SE, method = "REML", test = "knha"), silent = TRUE)
  if (inherits(h, "try-error")) return(NULL)
  # 无 HK 校正的 REML 拟合：只为取**常规** se 与 τ² 供预测区间使用，并留一列 rml_P
  # 作为「REML 单独 vs REML+HK」的对照（正文据此说明位移来自 HK 而非 REML）。
  m <- try(metafor::rma(yi = sub$logFC, sei = sub$SE, method = "REML"), silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  k <- q$k
  # 预测区间沿用正文已声明的 Higgins 式：beta_hat ± t(0.975, k−2)·sqrt(se^2 + tau²)，
  # se 取**常规** RE se（不是 HK 的 se_HK）——se_HK 已含队列间离散，再加 τ² 会重复计入。
  # k=2 时 k−2=0，无 t 分位，预测区间按定义不存在，落 NA（正文已如此声明）。
  pi_lb <- if (k >= 3) unname(h$beta[1] - qt(0.975, k - 2) * sqrt(m$se^2 + m$tau2)) else NA_real_
  pi_ub <- if (k >= 3) unname(h$beta[1] + qt(0.975, k - 2) * sqrt(m$se^2 + m$tau2)) else NA_real_
  list(ivw_logFC = unname(r$beta[1]), ivw_SE = unname(r$se), ivw_P = unname(r$pval),
       ivw_CI_low = unname(r$ci.lb), ivw_CI_high = unname(r$ci.ub),
       re_logFC = unname(q$beta[1]), re_SE = unname(q$se), re_P = unname(q$pval),
       re_CI_low = unname(q$ci.lb), re_CI_high = unname(q$ci.ub),
       hk_logFC = unname(h$beta[1]), hk_SE = unname(h$se), hk_P = unname(h$pval),
       hk_CI_low = unname(h$ci.lb), hk_CI_high = unname(h$ci.ub), hk_tau2 = unname(h$tau2),
       pi_low = pi_lb, pi_high = pi_ub,
       rml_logFC = unname(m$beta[1]), rml_P = unname(m$pval),
       Q = unname(q$QE), Q_p = unname(q$QEp),
       I2 = unname(100 * pmax(0, (q$QE - (q$k - 1)) / q$QE)), tau2 = unname(q$tau2), k_meta = q$k)
}

res <- list(); loo <- list()
for (layer in names(LAYER)) {
  cohorts <- LAYER[[layer]]
  cat(sprintf("\n========== %s层（%s）==========\n", layer, paste(cohorts, collapse = ", ")))
  for (gene in GENES) {
    fp <- fit_pooled(gene, cohorts)
    if (is.null(fp)) next
    mi <- meta_ivw(gene, cohorts)
    row <- data.frame(layer = layer, gene = gene, k = fp$k, n = fp$n,
                       pooled_logFC = fp$logFC, pooled_SE = fp$SE,
                       pooled_CI_low = fp$CI_low, pooled_CI_high = fp$CI_high,
                       pooled_P = fp$P, pooled_SE_type = "HC3",
                       cohorts_used = paste(fp$used, collapse = ";"),
                       cohort_logmeandiff = fp$cohort_logmeandiff,
                       stringsAsFactors = FALSE)
    if (!is.null(mi)) row <- cbind(row, as.data.frame(mi))
    res[[length(res) + 1]] <- row
    # 留一法：逐个剔除队列重拟合
    if (fp$k >= 3) {
      for (drop in cohorts) {
        f2 <- fit_pooled(gene, setdiff(cohorts, drop))
        if (is.null(f2)) next
        loo[[length(loo) + 1]] <- data.frame(
          layer = layer, gene = gene, dropped = drop, k = f2$k,
          logFC = f2$logFC, SE = f2$SE, CI_low = f2$CI_low, CI_high = f2$CI_high,
          P = f2$P, stringsAsFactors = FALSE)
      }
    }
  }
}

pooled <- do.call(rbind, res)

# 覆盖率对账：逐队列里有 ≥2 个队列的基因，必须都产出一行；产出行里记的 cohorts_used
# 必须等于该基因在该层的真实覆盖。上面符号解析那处缺陷的两种形态（整基因消失、
# 悄悄少用一个队列）都被这一条抓住——只查 k 不够，因为「真实覆盖 2 / 用了 2」与
# 「真实覆盖 3 / 用了 2」报出的 k 都是 2。
.exp <- unique(eff[eff$dataset %in% unlist(LAYER), c("gene", "dataset")])
for (layer in names(LAYER)) {
  e <- .exp[.exp$dataset %in% LAYER[[layer]], ]
  for (g in unique(e$gene)) {
    want <- sort(e$dataset[e$gene == g])
    if (length(want) < 2) next
    r <- pooled[pooled$layer == layer & pooled$gene == g, ]
    if (!nrow(r)) stop(sprintf("%s层 %s：逐队列覆盖 %d 个队列，却未产出合并行", layer, g, length(want)))
    got <- sort(strsplit(r$cohorts_used[1], ";")[[1]])
    if (!identical(got, want))
      stop(sprintf("%s层 %s：进入拟合的队列 [%s] ≠ 逐队列覆盖 [%s]",
                   layer, g, paste(got, collapse = ","), paste(want, collapse = ",")))
  }
}
cat(sprintf("覆盖率对账：两层共 %d 行全部与逐队列覆盖一致\n", nrow(pooled)))

# 不变量检查：合并估计是各队列效应的加权组合，必须落在**同一估计量**的逐队列取值范围内。
# 判据必须用模型实际估的量（对数均值差），不能用 B1 里 GSE89408 的 DESeq2 系数——
# 二者对 RNA-seq 是**不同的量**（均值差 vs 计数 GLM 系数，可差到符号相反）。
# 用错判据会把「估计量冲突」误报成「参数化写错」。
parse_coh <- function(s) as.numeric(sub(".*=", "", strsplit(s, "; ")[[1]]))
PER <- read.csv(file.path(OUT, "逐队列效应_力学基因.csv"), stringsAsFactors = FALSE)
bad <- 0L
for (i in seq_len(nrow(pooled))) {
  sub <- parse_coh(pooled$cohort_logmeandiff[i])
  if (!length(sub) || any(!is.finite(sub))) next
  if (pooled$pooled_logFC[i] < min(sub) - 1e-6 || pooled$pooled_logFC[i] > max(sub) + 1e-6) {
    bad <- bad + 1L
    cat(sprintf("  ⚠️ %s/%s 合并 %+.3f 越出同估计量逐队列范围 [%+.3f, %+.3f]\n",
        pooled$layer[i], pooled$gene[i], pooled$pooled_logFC[i], min(sub), max(sub)))
  }
}
cat(sprintf("\n合并估计越出「同估计量」逐队列范围的行数：%d / %d\n", bad, nrow(pooled)))
stopifnot(bad == 0L)   # 越界即参数化有问题，不得继续

# 同一不变式覆盖**全部四个** Meta 估计（IVW-FE、RE-DL、REML、REML+HK）：都是以 B1 逐队列
# logFC 为输入的加权平均，权重 1/vi 与 1/(vi+tau2) 全为正，故都必须落在**同一批 B1 取值**
# 的范围内。HK 换的是方差估计量而非点估计（beta_hat 仍为 Σw·y/Σw），故这条对它同样成立
# ——把它一并纳入判据，正是为了在点估计被写错时立刻响。
# 这里换用 B1 的取值作判据，与上一段用对数均值差判据是两件事：上一段查的是 pooled 那一列，
# 本段查的是 Meta 各列。判据取错会一边漏报一边误报。
bad_meta <- 0L
for (i in seq_len(nrow(pooled))) {
  sub <- PER$logFC[PER$gene == pooled$gene[i] & PER$dataset %in% LAYER[[pooled$layer[i]]]]
  sub <- sub[is.finite(sub)]
  if (length(sub) < 2) next
  for (col in c("ivw_logFC", "re_logFC", "rml_logFC", "hk_logFC")) {
    v <- pooled[[col]][i]
    if (!is.finite(v)) next
    if (v < min(sub) - 1e-6 || v > max(sub) + 1e-6) {
      bad_meta <- bad_meta + 1L
      cat(sprintf("  ⚠️ %s/%s %s=%+.3f 越出逐队列范围 [%+.3f, %+.3f]\n",
          pooled$layer[i], pooled$gene[i], col, v, min(sub), max(sub)))
    }
  }
}
cat(sprintf("Meta 四估计（IVW-FE / RE-DL / REML / REML+HK）越出逐队列范围的行数：%d / %d\n",
            bad_meta, 4 * sum(is.finite(pooled$ivw_logFC))))
stopifnot(bad_meta == 0L)

# 估计量冲突标记：合并模型（对数均值差）与逐队列汇报值（DESeq2/limma 系数）方向相反。
# 这只可能出现在发现层——复制层三个队列都是微阵列，两者本就是同一个量。
pooled$estimand_conflict <- FALSE
for (i in seq_len(nrow(pooled))) {
  sub <- PER$logFC[PER$gene == pooled$gene[i] & PER$dataset %in% LAYER[[pooled$layer[i]]]]
  if (!length(sub)) next
  if (all(sign(sub) == sign(sub[1])) && sign(pooled$pooled_logFC[i]) != sign(sub[1]))
    pooled$estimand_conflict[i] <- TRUE
}
cat(sprintf("估计量冲突（合并方向与各队列汇报值方向相反）：%d 行\n", sum(pooled$estimand_conflict)))
if (any(pooled$estimand_conflict))
  print(pooled[pooled$estimand_conflict,
               c("layer", "gene", "pooled_logFC", "ivw_logFC", "cohort_logmeandiff")], row.names = FALSE)

# 主估计列（2026-10-01 改；依据第二轮审稿意见第 1 项）。
# **两层同一估计量：REML + Hartung–Knapp**，作用在 B1 的逐队列效应上。
# 为什么两层都能用它、而不再一层用 IVW-FE、一层用 ~0+group+cohort：
#   分层原本是「一层 RNA-seq+微阵列混平台、一层纯微阵列」的产物。逐队列 Meta 让每个队列
#   贡献**自己的**估计量，混平台不再是障碍——发现层当初不用 Meta 的唯一理由（~0+group+cohort
#   把 GSE89408 的 DESeq2 系数换成对数均值差，把 LOXL1 拉到 +0.010）在 Meta 下不存在。
#   于是两层可以用同一个估计量，估计量本身不再随层变化，这是本轮最实质的简化。
# 被降级的两个（**均不再参与主估计，只作描述性敏感性**）：
#   IVW-FE  → ivw_* 列。权重 1/vi，不给 τ² 留位置；本轮降级。
#   RE-DL   → re_* 列。矩估计 τ²；与 REML 数值相近但非同一估计量，留作对照。
# pooled(~0+group+cohort) → pooled_* 列。样本级 HC3，改变分析单位，自始就是
#   「估计策略比较」而非「模型比较」（§2.3(iv) 已如此限定），本轮起明确降为敏感性。
pooled$primary_method <- "REML+HK"
pooled$primary_logFC   <- pooled$hk_logFC
pooled$primary_CI_low  <- pooled$hk_CI_low
pooled$primary_CI_high <- pooled$hk_CI_high
pooled$primary_P       <- pooled$hk_P

# 敏感性列的标注。三列各自对照的对象必须写清，否则会被读成三个主估计。
pooled$re_method <- "DL"          # re_* 列确实是 DL 的数值，标签不得写成 REML
pooled$re_layer_note <- paste0(
  "主估计为 REML+HK（hk_* 列，两层同一估计量）；",
  "ivw_* 为固定效应、re_* 为 DL、pooled_* 为样本级 HC3，三者均降为描述性敏感性；",
  "rml_* 为不加 HK 校正的 REML，仅供核对「位移来自 HK 而非 REML」")
# 预测区间的可算性逐行标注：k=2 时 k−2=0，t 分位不存在，区间按定义不存在。
pooled$pi_status <- ifelse(is.na(pooled$pi_low), "not defined (k=2, k-2=0)", "reported")

pooled_f <- pooled[pooled$layer == "发现", ]; pooled_r <- pooled[pooled$layer == "复制", ]
write.csv(pooled_f, file.path(OUT, "合并效应_发现层.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(pooled_r, file.path(OUT, "合并效应_复制层.csv"), row.names = FALSE, fileEncoding = "UTF-8")
LOO <- do.call(rbind, loo)
write.csv(LOO, file.path(OUT, "留一法_合并效应.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- 3. 控制台：面板基因 + 交联基因 ----------
show <- function(layer, genes) {
  cat(sprintf("\n---------- %s层（主估计：REML+HK，两层同一估计量）----------\n", layer))
  # 首位两列必须印 primary_*，与标题声明的估计方法一致；印 pooled_* 会让读者以为
  # 主估计就是 pooled 那一行，而发现层两者方向可差到符号相反（如 LOXL1）。
  cat(sprintf("%-9s %9s %8s %18s %9s | %9s %9s %7s %6s\n",
              "gene", "主估计logFC", "SE", "95%CI", "P",
              "对照logFC", "对照 P", "I²%", "k"))
  for (g in genes) {
    r <- pooled[pooled$layer == layer & pooled$gene == g, ]
    if (!nrow(r)) { cat(sprintf("%-9s  （可用队列不足，无法合并）\n", g)); next }
    cat(sprintf("%-9s %+9.3f %8.3f [%+7.3f,%+7.3f] %9.3g | %+9.3f %9.3g %7.1f %6d\n",
                g, r$primary_logFC, r$pooled_SE, r$primary_CI_low, r$primary_CI_high, r$primary_P,
                r$pooled_logFC, r$pooled_P,
                if ("I2" %in% names(r)) r$I2 else NA, r$k))
  }
}
for (layer in names(LAYER)) { show(layer, PANEL); show(layer, setdiff(CROSSLINK, PANEL)) }

## ---------- 4. 一致性检查：主分析与 IVW 对照 ----------
cat("\n========== 主分析(~0+group+cohort) 与 IVW 对照 ==========\n")
cmp <- pooled[!is.na(pooled$ivw_logFC), ]
sgn <- sign(cmp$pooled_logFC) == sign(cmp$ivw_logFC)
cat(sprintf("  可比基因数 %d；方向一致 %d（%.1f%%）\n", nrow(cmp), sum(sgn), 100 * mean(sgn)))
cat(sprintf("  |ΔlogFC| 中位 %.3f  最大 %.3f（%s）\n",
            median(abs(cmp$pooled_logFC - cmp$ivw_logFC)),
            max(abs(cmp$pooled_logFC - cmp$ivw_logFC)),
            cmp$gene[which.max(abs(cmp$pooled_logFC - cmp$ivw_logFC))]))
cat("  方向不一致的基因：\n")
if (any(!sgn)) print(cmp[!sgn, c("layer", "gene", "pooled_logFC", "ivw_logFC")], row.names = FALSE) else cat("    （无）\n")

## ---------- 4b. 随机效应敏感性：换成 RE-DL 会不会翻结论 ----------
# 这一步只回答一个问题：I² 偏高时主估计的结论是否依赖「固定效应」这个选择。
# 判据分两级：方向是否变（严重）、显著性是否跨 0.05（次要，且不必然算「翻转」——
# 见笔记 p-value-shift-is-not-verdict-flip）。
cat("\n========== 随机效应敏感性（RE-DL vs IVW-FE）==========\n")
re_cmp <- pooled[is.finite(pooled$ivw_logFC) & is.finite(pooled$re_logFC), ]
re_sgn <- sign(re_cmp$ivw_logFC) == sign(re_cmp$re_logFC)
re_flip <- sum(!re_sgn)
cat(sprintf("  可比行 %d（发现层 %d / 复制层 %d）\n", nrow(re_cmp),
            sum(re_cmp$layer == "发现"), sum(re_cmp$layer == "复制")))
cat(sprintf("  方向改变的行：%d\n", re_flip))
if (re_flip) print(re_cmp[!re_sgn, c("layer", "gene", "ivw_logFC", "re_logFC", "I2", "tau2")],
                   row.names = FALSE)
sig_ivw <- re_cmp$ivw_P < 0.05; sig_re <- re_cmp$re_P < 0.05
cat(sprintf("  显著性跨 0.05 的行：%d（IVW 显著而 RE 不显著 %d；RE 显著而 IVW 不显著 %d）\n",
            sum(sig_ivw != sig_re), sum(sig_ivw & !sig_re), sum(!sig_ivw & sig_re)))
if (any(sig_ivw != sig_re))
  print(re_cmp[sig_ivw != sig_re,
               c("layer", "gene", "ivw_logFC", "ivw_P", "re_logFC", "re_P", "I2")], row.names = FALSE)
cat("  面板基因逐条对照：\n")
cat(sprintf("  %-9s %-4s %10s %10s %10s %10s %8s\n",
            "gene", "层", "IVW-FE", "IVW P", "RE-DL", "RE P", "I2%"))
for (L in c("发现", "复制")) for (g in PANEL) {
  r <- re_cmp[re_cmp$layer == L & re_cmp$gene == g, ]
  if (!nrow(r)) { cat(sprintf("  %-9s %-4s  （无可比行）\n", g, L)); next }
  cat(sprintf("  %-9s %-4s %+10.3f %10.3g %+10.3f %10.3g %8.1f\n",
              g, L, r$ivw_logFC[1], r$ivw_P[1], r$re_logFC[1], r$re_P[1], r$I2[1]))
}

## ---------- 5. 留一法：合并效应是否被单一队列支配 ----------
cat("\n========== 留一法（面板基因）==========\n")
for (layer in names(LAYER)) for (g in PANEL) {
  L <- LOO[LOO$layer == layer & LOO$gene == g, ]
  if (!nrow(L)) next
  full <- pooled$pooled_logFC[pooled$layer == layer & pooled$gene == g]
  cat(sprintf("  %-4s %-9s 全模型 %+7.3f | %s\n", layer, g, full,
              paste(sprintf("剔%s %+7.3f", sub("GSE", "", L$dropped), L$logFC), collapse = "  ")))
}

cat(sprintf("\n已写出：\n  02_结果表/合并效应_发现层.csv  %d 行\n  02_结果表/合并效应_复制层.csv  %d 行\n  02_结果表/留一法_合并效应.csv  %d 行\n",
            sum(pooled$layer == "发现"), sum(pooled$layer == "复制"), nrow(LOO)))
