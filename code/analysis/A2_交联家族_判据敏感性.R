# A-2.1 / A-2.2 交联家族 8 基因同框架比较表 + 判据敏感性汇总
#
# 本脚本**不做任何新的估计**。它只把既有产物（合并效应_发现层/复制层.csv、57基因总览.csv）
# 重新排列成可直接进主表的形状，并量化「优先级随模型而变」。所有数值均可在上述 CSV 里逐格查到。
#
# ---- 2026-10-01 口径变更：两层统一为 REML+HK 主估计 ----
# 本轮之前两层的主估计**不同**：发现层 = IVW-FE，复制层 = pooled(~0+group+cohort)（HC3）。
# 本轮按审稿意见把两层的主估计统一为 **REML + Hartung–Knapp**（B2 的 `hk_*` 列，
# 亦即 B2 写出的 `primary_*`）；IVW-FE（`ivw_*`）与 DL（`re_*`）降为描述性敏感性。
# 于是本表两层列集合**完全相同**（旧版两层不同，读表时极易错配），
# 而 §三 的跨线对照改为两对：「IVW-FE → REML+HK」（同一逐队列输入的模型比较）
# 与「样本级 pooled → 队列级 REML+HK」（同时改分析单位与方差估计量的估计策略比较）。
# REML（`rml_*`，不含 HK）一并保留，用来把位移拆成「REML 的份」与「HK 的份」。
#
# ---- 一个必须写进正文的性质（实测，非推断） ----
# 主估计与固定效应对照换的不只是 SE：随机效应会改变权重，**点估计也随之变化**。
# 故「按 |logFC| 的名次在固定效应/随机效应下不变」是错的——两类名次都要报。
# 另须写明：HK **不是一致更保守**，τ²≈0 的格子区间反而收窄（见 §一 与审计记录）。
#
# 输入：$A2_RES/{合并效应_发现层,合并效应_复制层,57基因总览}.csv（默认落在旧树，见下）
# 输出：02_结果表/A2_交联家族_判据敏感性.csv      核心表（可直接进主表）
#       02_结果表/A2_名次随模型变动.csv           逐模型名次
#       01_审计/A2_判据敏感性_结果.md             三条读数的**实测**版本
suppressPackageStartupMessages({ library(stats) })

# 注意：R 不支持相邻字符串字面量的隐式拼接（那是 Python 的写法），必须显式 paste0
BASE <- paste0("/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/",
               "Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/",
               "01.材料/01.博士材料/10.毕业/01.小论文/04.立项/",
               "滑膜基质刚度生成与刚度感知 0912/08_修订_2026-09-17")
# 2026-10-01：落点改由环境变量控制。默认仍是旧树（历史可复现），但本轮修订必须显式传
# A2_RES/A2_AUD 指向修订工作树，否则读到的还是旧口径的合并表（口径断言会先拦下）。
RES  <- Sys.getenv("A2_RES", unset = file.path(BASE, "02_结果表"))
AUD  <- Sys.getenv("A2_AUD", unset = file.path(BASE, "01_审计"))
dir.create(AUD, showWarnings = FALSE, recursive = TRUE)

CROSSLINK <- c("LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4", "PLOD1", "PLOD2", "TGM2")
PANEL     <- c("LOXL1", "VGLL4", "SDC4", "SERPINE1", "ITGA5", "SAV1")
ALPHA     <- 0.05

# 主估计口径的**字面量**（与 C0 口径回源核对一致；不从数据现算，否则断言恒成立）
# 2026-10-01：两层统一为 REML+HK，估计量不再随层变化。
EXPECT_METHOD <- c(发现 = "REML+HK", 复制 = "REML+HK")

dis <- read.csv(file.path(RES, "合并效应_发现层.csv"), stringsAsFactors = FALSE)
rep <- read.csv(file.path(RES, "合并效应_复制层.csv"), stringsAsFactors = FALSE)
ov  <- read.csv(file.path(RES, "57基因总览.csv"),       stringsAsFactors = FALSE)

## ---------- 0. 口径断言：对不上就不许出表 ----------
stopifnot(all(dis$primary_method == EXPECT_METHOD[["发现"]]),
          all(rep$primary_method == EXPECT_METHOD[["复制"]]))
NEED <- c("ivw_logFC","ivw_P","re_logFC","re_P","I2","k",
          "hk_logFC","hk_P","hk_CI_low","hk_CI_high","hk_tau2","pi_low","pi_high",
          "primary_logFC","primary_P","primary_CI_low","primary_CI_high",
          "pooled_logFC","pooled_P")
stopifnot(all(NEED %in% names(dis)), all(NEED %in% names(rep)))
stopifnot(all(dis$primary_logFC == dis$hk_logFC), all(rep$primary_logFC == rep$hk_logFC))
stopifnot(all(CROSSLINK %in% ov$gene))

say <- function(...) cat(sprintf(...), "\n")
ln  <- character(0)
rec <- function(s = "") { cat(s, "\n"); ln <<- c(ln, s) }

rec("# A-2 交联家族判据敏感性（实测结果）")
rec()
rec("> 由 `03_分析/scripts/A2_交联家族_判据敏感性.R` 生成。**本脚本不做新估计**，")
rec("> 全部数值来自 `合并效应_{发现,复制}层.csv` 与 `57基因总览.csv`，可逐格回查。")
rec()

## ---------- 1. 核心表：8 交联基因同框架比较 ----------
f <- function(v) ifelse(is.na(v), NA_real_, as.numeric(v))
pick <- function(tab, g, cols) {
  k <- match(g, tab$gene)
  if (is.na(k)) return(setNames(rep(NA_real_, length(cols)), cols))
  setNames(vapply(cols, function(c) f(tab[[c]][k]), numeric(1)), cols)
}

# 两层列集合完全相同（旧版两层不同：发现层给 ivw_*、复制层给 pooled_*，读表时极易错配）。
# 主估计 = `primary_*`（= `hk_*`，REML+HK，由 B2 写死并已在上方断言）；`ivw_*` 与 `re_*`
# 一律是被降级为描述性敏感性的对照。
LAYER_COLS <- c("k",
                "primary_logFC", "primary_P", "primary_CI_low", "primary_CI_high",
                "pi_low", "pi_high", "hk_tau2",
                "ivw_logFC", "ivw_P", "ivw_CI_low", "ivw_CI_high",
                "re_logFC", "re_P", "I2",
                "rml_logFC", "rml_P",   # REML 但不含 HK：用来把位移拆成「REML 的份」与「HK 的份」
                # ★ 2026-10-01：pooled_* 加入。改口径后**复制层的声明敏感性估计量**是
                # pooled(~0+group+cohort)+HC3（发现层是 IVW-FE），D2 的图 4 要按层各画
                # 「主估计 + 该层声明敏感性」两组列 —— 主估计改指 hk 之后，复制层那组若不从
                # 这里取，D2 只能读到 hk，会把敏感性列画成主估计的复制品（图仍照常生成）。
                # 发现层的 pooled_* 一并落盘作表内对照列，但**不上图**：它估的是「对数均值差」，
                # 对 GSE89408 的 DESeq2 系数是不同的量（B7 文件头 ② 已记）。
                "pooled_logFC", "pooled_P")

core <- do.call(rbind, lapply(CROSSLINK, function(g) {
  d <- pick(dis, g, LAYER_COLS); names(d) <- paste0("disc_", LAYER_COLS)
  r <- pick(rep, g, LAYER_COLS); names(r) <- paste0("rep_",  LAYER_COLS)
  data.frame(gene = g, in_panel = g %in% PANEL, as.list(d), as.list(r),
             stringsAsFactors = FALSE)
}))
# 复制层缺行时 pick() 全 NA；以 rep_k 为准显式标出，避免被读成「无差异」
core$rep_available <- !is.na(core$rep_k)
# 把口径字面量随表一起落盘：下游（D2/B9）据以断言「读到的这一列确实是该层主估计」
core$disc_primary_method <- EXPECT_METHOD[["发现"]]
core$rep_primary_method  <- EXPECT_METHOD[["复制"]]

fm  <- function(x, d = 3) ifelse(is.na(x), "—", formatC(x, format = "f", digits = d, flag = "+"))
fP  <- function(x) ifelse(is.na(x), "—", formatC(x, format = "g", digits = 3))
# I² 是百分数且常恰为 0：用 %g 会把 0 印成「   0」，读起来像缺失，故单列一个格式化器
fI2 <- function(x) ifelse(is.na(x), "—", formatC(x, format = "f", digits = 1))
fCI <- function(lo, hi) ifelse(is.na(lo) | is.na(hi), "—",
                               sprintf("(%s, %s)", fm(lo), fm(hi)))

rec("## 一、交联家族 8 基因同框架比较（核心表）")
rec()
rec(sprintf("两层主估计同为 **%s**（`method=\"REML\", test=\"knha\"`），逐格取自 `primary_*` 列。",
            EXPECT_METHOD[["发现"]]))
rec("`FE` = IVW 固定效应 Meta；`DL` = DerSimonian–Laird 随机效应 Meta。二者在本轮**均降为")
rec("描述性敏感性**，不再承担任何主结论。预测区间（PI）按 Higgins 口径 `β ± t(0.975,k−2)·√(se²+τ²)`，")
rec("k=2 时无定义（k−2=0）。")
rec()
rec("### 1a 主估计（REML+HK）")
rec()
rec("| 基因 | 面板 | 发现 k | 发现 logFC | 发现 95%CI | 发现 P | 发现 PI | 发现 I² | 复制 k | 复制 logFC | 复制 95%CI | 复制 P | 复制 PI | 复制 I² |")
rec("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(core))) {
  z <- core[i, ]
  rec(sprintf("| **%s** | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
              z$gene, ifelse(z$in_panel, "是", ""),
              ifelse(is.na(z$disc_k), "—", z$disc_k),
              fm(z$disc_primary_logFC), fCI(z$disc_primary_CI_low, z$disc_primary_CI_high),
              fP(z$disc_primary_P), fCI(z$disc_pi_low, z$disc_pi_high), fI2(z$disc_I2),
              ifelse(is.na(z$rep_k), "—（无行）", z$rep_k),
              fm(z$rep_primary_logFC), fCI(z$rep_primary_CI_low, z$rep_primary_CI_high),
              fP(z$rep_primary_P), fCI(z$rep_pi_low, z$rep_pi_high), fI2(z$rep_I2)))
}
rec()
rec("### 1b 被降级的两个对照估计量（只列 P，供核对「换估计量改变了哪几格判定」）")
rec()
rec("| 基因 | 发现 FE P | 发现 DL P | 发现 主 P | 复制 FE P | 复制 DL P | 复制 主 P |")
rec("|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(core))) {
  z <- core[i, ]
  rec(sprintf("| **%s** | %s | %s | %s | %s | %s | %s |",
              z$gene, fP(z$disc_ivw_P), fP(z$disc_re_P), fP(z$disc_primary_P),
              fP(z$rep_ivw_P), fP(z$rep_re_P), fP(z$rep_primary_P)))
}
rec()
# 一个只在两列同时存在时才能验的性质：REML 单独（不含 HK）与 DL 的差异远小于 HK 带来的差异。
# 若 REML≡DL 逐格成立，则全部位移都来自 HK 这一层校正 —— 这句话必须在正文里说准。
# 列名一律小写 re_/ivw_/rml_。写成 core$disc_RE_logFC（大写）时 `$` 既非精确匹配也非
# 前缀匹配，会**静默返回 NULL** 而不是报错，下游 max() 于是给出 -Inf（实测踩到）。
# 故先按字面量取出并断言长度与有限性 —— 不许让「取不到的列」伪装成一个数。
pull2 <- function(a, b) {
  stopifnot(is.numeric(a), is.numeric(b), length(a) == length(b), length(a) > 0)
  v <- c(a, b)
  stopifnot(any(is.finite(v)))
  v
}
d_rml_dl <- max(abs(pull2(core$disc_primary_logFC, core$rep_primary_logFC) -
                    pull2(core$disc_re_logFC,     core$rep_re_logFC)), na.rm = TRUE)
rec(sprintf("REML+HK 与 DL 的**点估计**最大绝对差（家族 8 基因，两层合看）：%.4f。", d_rml_dl))
rec("HK 只换方差估计量、不动点估计，故该差值完全来自 REML 与 DL 对 τ² 的估计不同；")
rec("而**显著性判定**的位移另由 HK 那一层贡献 —— 两者不得混为一句。")
rec()

## ---------- 2. 名次随模型而变 ----------
rank_tab <- function(tab, layer) {
  # 在「交联家族内」与「全部基因内」分别按 |logFC| 与按 P 排名，FE 与主估计（REML+HK）各一套。
  # 列名沿用 FE_rank / RE_rank（schema 不变，下游 D1/D2/F2/F3/F8 无需改列名），
  # 但 RE_rank 自 2026-10-01 起指的是**该层主估计 REML+HK**，不再是 DL；
  # 故另落一列 re_estimator 把这一点写进产物本身，避免下游按旧含义读。
  out <- list()
  for (setname in c("交联家族(8)", "全部基因")) {
    genes <- if (setname == "交联家族(8)") CROSSLINK else tab$gene
    x <- tab[tab$gene %in% genes, ]
    for (basis in c("abs_logFC", "P")) {
      v_fe <- if (basis == "abs_logFC") abs(x$ivw_logFC)     else x$ivw_P
      v_re <- if (basis == "abs_logFC") abs(x$primary_logFC) else x$primary_P
      # 名次含义统一为「1 = 最靠前」：按 |logFC| 时降序（取负），按 P 时升序。
      # 直接 rank(v) 会让 |logFC| 的名次 1 = 效应量最小，与读者预期相反（实测踩过）。
      k_fe <- if (basis == "abs_logFC") -v_fe else v_fe
      k_re <- if (basis == "abs_logFC") -v_re else v_re
      r_fe <- rank(ifelse(is.na(k_fe), Inf, k_fe), ties.method = "min")
      r_re <- rank(ifelse(is.na(k_re), Inf, k_re), ties.method = "min")
      out[[length(out) + 1]] <- data.frame(
        layer = layer, set = setname, basis = basis, gene = x$gene,
        FE_rank = r_fe, RE_rank = r_re, n_in_set = length(genes),
        delta = r_re - r_fe, re_estimator = EXPECT_METHOD[[layer]],
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

rk <- rbind(rank_tab(dis, "发现"), rank_tab(rep, "复制"))
rk_cross <- rk[rk$gene %in% CROSSLINK, ]

rec("## 二、名次随模型而变（A-2.2）")
rec()
rec("「名次」有两套互不相同的判据：按 |logFC| 与按 P。**两套都要报**——")
rec(sprintf("因为对照的两端是 %s（主估计）与 IVW 固定效应，随机效应换的不只是 SE：",
            EXPECT_METHOD[["发现"]]))
rec("权重随 τ² 改变，点估计随之变化，故按 |logFC| 的名次同样会动。")
rec("本表 `RE_rank` 列自 2026-10-01 起指**该层主估计 REML+HK**（不再是 DL）；产物内 `re_estimator`")
rec("列记录了该取值，下游不得按旧含义读。")
rec()
rec("名次一律「1 = 最靠前」（按 |logFC| 取降序、按 P 取升序）。复制层家族内只有 6 个基因")
rec("——LOXL3、LOXL4 在复制层无合并行（见 §三 与读数 3），故其复制层名次为空，不是排在末位。")
rec()
for (ly in c("发现", "复制")) {
  rec(sprintf("### %s层（交联家族 8 基因内）", ly))
  rec()
  rec("| 基因 | FE 名次(|logFC|) | RE 名次(|logFC|) | 变动 | FE 名次(P) | RE 名次(P) | 变动 |")
  rec("|---|---|---|---|---|---|---|")
  for (g in CROSSLINK) {
    a <- rk_cross[rk_cross$layer == ly & rk_cross$set == "交联家族(8)" & rk_cross$basis == "abs_logFC" & rk_cross$gene == g, ]
    p <- rk_cross[rk_cross$layer == ly & rk_cross$set == "交联家族(8)" & rk_cross$basis == "P" & rk_cross$gene == g, ]
    if (!nrow(a)) next
    sgn <- function(d) ifelse(is.na(d), "—", sprintf("%+d", d))
    rec(sprintf("| %s | %s | %s | %s | %s | %s | %s |", g,
                ifelse(is.infinite(a$FE_rank) | is.na(a$FE_rank), "—", a$FE_rank),
                ifelse(is.infinite(a$RE_rank) | is.na(a$RE_rank), "—", a$RE_rank),
                sgn(a$delta),
                ifelse(is.infinite(p$FE_rank) | is.na(p$FE_rank), "—", p$FE_rank),
                ifelse(is.infinite(p$RE_rank) | is.na(p$RE_rank), "—", p$RE_rank),
                sgn(p$delta)))
  }
  rec()
}

## ---------- 3. 显著性跨线计数（判据依赖性的量化） ----------
cross_count <- function(tab, layer, basis) {
  if (basis == "FE_vs_primary（Meta 对：IVW-FE → REML+HK）") { a <- tab$ivw_P; b <- tab$primary_P }
  else if (basis == "pooled_vs_primary（估计策略对：样本级 HC3 → 队列级 REML+HK）") {
    if (!"pooled_P" %in% names(tab)) stop("该层没有 pooled_* 列（只有复制层有样本级 pooled 模型）")
    a <- tab$pooled_P; b <- tab$primary_P
  }
  else stop("未知 basis")
  k <- tab$gene %in% CROSSLINK & !is.na(a) & !is.na(b) & is.finite(a) & is.finite(b)
  crossed <- k & ((a < ALPHA) != (b < ALPHA))
  data.frame(layer = layer, basis = basis, n_eval = sum(k), n_crossed = sum(crossed),
             genes = paste(tab$gene[crossed], collapse = "、"),
             P_min = min(c(a[crossed], b[crossed]), na.rm = TRUE),
             P_max = max(c(a[crossed], b[crossed]), na.rm = TRUE),
             stringsAsFactors = FALSE)
}
cc <- rbind(
  cross_count(dis, "发现", "FE_vs_primary（Meta 对：IVW-FE → REML+HK）"),
  cross_count(rep, "复制", "FE_vs_primary（Meta 对：IVW-FE → REML+HK）"),
  cross_count(rep, "复制", "pooled_vs_primary（估计策略对：样本级 HC3 → 队列级 REML+HK）"))

rec("## 三、显著性跨 0.05 的基因（判据依赖性的量化）")
rec()
rec("> 报范围取**跨线那批基因的 min/max**，不拿全体或子集冒充上下界（记忆 `p-value-shift-is-not-verdict-flip`）。")
rec()
rec("| 层 | 对照 | 可评基因 | 跨线基因数 | 是哪些 | 跨线值的 min–max |")
rec("|---|---|---|---|---|---|")
for (i in seq_len(nrow(cc))) {
  z <- cc[i, ]
  rec(sprintf("| %s | %s | %d | %d | %s | %.4f – %.4f |", z$layer, z$basis, z$n_eval,
              z$n_crossed, ifelse(z$genes == "", "（无）", z$genes), z$P_min, z$P_max))
}
rec()

## ---------- 3.5 位移归因：四估计量的显著数（全层与家族各一套） ----------
# 主估计换成 REML+HK 后，显著数从 29/27 掉到 3/4。读者必然要问「是 REML 换掉 DL 的功劳，
# 还是 HK 的功劳」。四列都在表里，故直接数出来 —— 这正是把 REML（rml_*，不含 HK）
# 一并留在产物里的理由。
sig_n <- function(tab, setname) {
  x <- if (setname == "全部") tab else tab[tab$gene %in% CROSSLINK, ]
  x <- x[!is.na(x$ivw_P), ]
  c(n = nrow(x),
    FE      = sum(x$ivw_P     < ALPHA),
    DL      = sum(x$re_P      < ALPHA, na.rm = TRUE),
    REML    = sum(x$rml_P     < ALPHA, na.rm = TRUE),
    REML_HK = sum(x$primary_P < ALPHA, na.rm = TRUE))
}
SIG <- rbind(发现_全部 = sig_n(dis, "全部"), 发现_家族 = sig_n(dis, "家族"),
             复制_全部 = sig_n(rep, "全部"), 复制_家族 = sig_n(rep, "家族"))

rec("## 三之二、位移归因：四个估计量各自的显著数（名义 P<0.05）")
rec()
rec("| 层/集合 | 可汇总 n | FE | DL | REML | REML+HK（主） |")
rec("|---|---|---|---|---|---|")
for (i in seq_len(nrow(SIG))) {
  rec(sprintf("| %s | %d | %d | %d | %d | **%d** |",
              rownames(SIG)[i], SIG[i, "n"], SIG[i, "FE"], SIG[i, "DL"],
              SIG[i, "REML"], SIG[i, "REML_HK"]))
}
rec()
# 自洽性：若 REML 列与 DL 列逐行同数，则位移全部来自 HK；否则必须在正文分开陈述。
rec(sprintf("REML（不含 HK）与 DL 的显著数逐行%s。",
            ifelse(all(SIG[, "REML"] == SIG[, "DL"]),
                   "**相同** —— 位移全部来自 HK 那一层校正，不是 REML 换 DL 换来的",
                   "**不同** —— 位移同时来自 REML 与 HK 两层，正文须分开陈述")))
rec()

## ---------- 4. 读数（实测版） ----------
n_rep_cross     <- cc$n_crossed[cc$layer == "复制" & grepl("^FE_vs_primary", cc$basis)]
rep_cross_genes <- cc$genes[cc$layer == "复制" & grepl("^FE_vs_primary", cc$basis)]
n_pool_cross     <- cc$n_crossed[cc$layer == "复制" & grepl("^pooled_vs_primary", cc$basis)]
pool_cross_genes <- cc$genes[cc$layer == "复制" & grepl("^pooled_vs_primary", cc$basis)]

# 家族内两层都通过**主估计**的是哪些（旧版的同名读数用的是 DL；此处按新主估计重算）
both_primary <- core$gene[!is.na(core$disc_primary_P) & !is.na(core$rep_primary_P) &
                          core$disc_primary_P < ALPHA & core$rep_primary_P < ALPHA]
# 家族内「只在一层通过主估计」的，逐层点名 —— 这是「无候选跨两层」这句话的实际内容
disc_only <- core$gene[!is.na(core$disc_primary_P) & core$disc_primary_P < ALPHA &
                       (is.na(core$rep_primary_P) | core$rep_primary_P >= ALPHA)]
rep_only  <- core$gene[!is.na(core$rep_primary_P) & core$rep_primary_P < ALPHA &
                       (is.na(core$disc_primary_P) | core$disc_primary_P >= ALPHA)]

pl <- core[core$gene == "PLOD1", ]
loxl4 <- core[core$gene == "LOXL4", ]
loxl1 <- core[core$gene == "LOXL1", ]

rec("## 四、读数（实测版）")
rec()
rec(sprintf("1. **复制层的显著性高度依赖估计量：样本级 pooled(HC3) → 队列级 REML+HK 跨 0.05 的有 %d 个**（%s）。",
            n_pool_cross, ifelse(pool_cross_genes == "", "（无）", pool_cross_genes)))
rec(sprintf("   跨线值落在 %.4f–%.4f。同层的 Meta 对内（IVW-FE → REML+HK）跨线 %d 个：%s。",
            cc$P_min[cc$layer == "复制" & grepl("^pooled_vs_primary", cc$basis)],
            cc$P_max[cc$layer == "复制" & grepl("^pooled_vs_primary", cc$basis)],
            n_rep_cross, ifelse(rep_cross_genes == "", "（无）", rep_cross_genes)))
rec("   两者含义不同，不得混写：前者同时改了**分析单位**（样本级→队列级）与方差估计量，")
rec("   后者只改**估计量**（同一逐队列输入）。")
rec()
rec(sprintf("2. **交联家族 8 基因内，两层都通过主估计（REML+HK）的是 %s**。",
            ifelse(length(both_primary) == 0, "**一个都没有**", paste(both_primary, collapse = "、"))))
rec(sprintf("   只在发现层通过：%s；只在复制层通过：%s。",
            ifelse(length(disc_only) == 0, "（无）", paste(disc_only, collapse = "、")),
            ifelse(length(rep_only)  == 0, "（无）", paste(rep_only,  collapse = "、"))))
rec(sprintf("   PLOD1 两层主估计 P：发现 %s、复制 %s —— 它是**复制层**的通过者，发现层未过。",
            fP(pl$disc_primary_P), fP(pl$rep_primary_P)))
rec(sprintf("   LOXL1 两层主估计 P：发现 %s、复制 %s —— 它是**发现层**的通过者，复制层未过。",
            fP(loxl1$disc_primary_P), fP(loxl1$rep_primary_P)))
rec("   （该结论的范围是交联家族 8 基因；面板基因另在 §一 表内，未并入本读数的计数。）")
rec()
rec("3. **LOXL4 的发现层效应量大于 LOXL1，但复制层完全不可合并**：")
rec(sprintf("   LOXL4 发现主估计 logFC=%s、P=%s（%s）；LOXL1 发现主估计 logFC=%s、P=%s（%s）。",
            fm(loxl4$disc_primary_logFC), fP(loxl4$disc_primary_P),
            ifelse(loxl4$disc_primary_logFC > 0, "上调", "下调"),
            fm(loxl1$disc_primary_logFC), fP(loxl1$disc_primary_P),
            ifelse(loxl1$disc_primary_logFC > 0, "上调", "下调")))
rec("   LOXL4 与 LOXL3 复制层无行：仅 1 个复制队列测到该基因，低于 B2 第 92 行「至少 2 队列」的合并门槛。")
rec("   ——不是「复制层无差异」，是**复制层不可合并**，两者不得混写。")
rec()

## ---------- 5. 落盘 ----------
write.csv(core, file.path(RES, "A2_交联家族_判据敏感性.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(rk, file.path(RES, "A2_名次随模型变动.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(cc, file.path(RES, "A2_显著性跨线计数.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
# 注意分母的范围：`n_eval` 是**交联家族内**的可评数（≤8），不是该层全部候选（57/52）。
# 写成「复制层可评 N 个基因」会让读者当成普查分母（记忆 `filtered-detail-table-is-not-a-census`）。
writeLines(c(ln, "", "---",
             sprintf("复制层**家族内**可评 %d 个基因；IVW-FE → REML+HK 跨线 %d 个。",
                     cc$n_eval[2], cc$n_crossed[2]),
             "**本脚本不做新估计；全部数值来自既有 CSV，可逐格回查。**"),
           file.path(AUD, "A2_判据敏感性_结果.md"))

cat("\n已写出：\n")
cat("  02_结果表/A2_交联家族_判据敏感性.csv\n")
cat("  02_结果表/A2_名次随模型变动.csv\n")
cat("  02_结果表/A2_显著性跨线计数.csv\n")
cat("  01_审计/A2_判据敏感性_结果.md\n")
