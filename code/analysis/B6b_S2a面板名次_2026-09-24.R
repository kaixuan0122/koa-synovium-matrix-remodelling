# B6b：补跑 S2a 缺格——GSE89408 面板基因按 |log2FC| 的队列内名次
#
# 为什么单独跑：第 5 步重建 S2 时发现，S2a 表「Within-cohort rank of panel by |log2FC|」
#   一格（旧值 4,535–16,154）在 B6 的 ml11 重跑里**没有被输出**——run_ml11 算了 de_ord
#   却没算 panel_rank（旧 11 号脚本 L154 有 `panel_rank <- match(panel, de_ord)`）。
#   按规则 4，未取回前该格不得沿用旧值。本脚本即为取回该格的真实运行。
#
# 方法同一性检验（关键）：本脚本**同时**在旧 colSums 矩阵上复算一遍。
#   若旧矩阵复现不出旧值 4,535–16,154，说明口径与旧脚本不同，则新矩阵上的值也不可用。
#   （旧脚本的 tt 与 de_ord 见 04_分析数据/scripts/11_machine_learning_by_cohort.R L151-154。）
#
# 输入：08_修订_2026-09-17/03_分析/outputs/GSE89408_校正矩阵_cpm尺度.rds（B6 mat 部分产出）
#       04_分析数据/processed/主分析_genelevel.rds（旧 colSums 矩阵）
#       08_修订_2026-09-17/03_分析/outputs/GSE89408_基因水平_双版本.rds（g89$grp 分组）
# 用法：Rscript B6b_S2a面板名次_2026-09-24.R
# 输出：03_分析/logs/B6b_S2a面板名次.log（由调用方重定向）
#       02_结果表/F11下游_ml11_面板名次.csv

suppressPackageStartupMessages({ library(limma) })

BASE <- "/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/01.材料/01.博士材料/10.毕业/01.小论文/04.立项/滑膜基质刚度生成与刚度感知 0912"
P    <- file.path(BASE, "04_分析数据/processed")
OUT  <- file.path(BASE, "08_修订_2026-09-17/03_分析/outputs")
RES  <- file.path(BASE, "08_修订_2026-09-17/02_结果表")

PANEL <- c("LOXL1", "ITGA5", "SDC4", "SERPINE1", "SAV1", "VGLL4")

## 分组：与 B6 run_ml11 同一来源
g89 <- readRDS(file.path(OUT, "GSE89408_基因水平_双版本.rds"))
y <- factor(ifelse(g89$grp == "OA", "OA", "Normal"), levels = c("Normal", "OA"))

## 校正矩阵：直接用 B6 mat 部分的落盘产物（已断言与 size 版同维且换算生效）
mat_new <- readRDS(file.path(OUT, "GSE89408_校正矩阵_cpm尺度.rds"))$mat
mat_old <- readRDS(file.path(P, "主分析_genelevel.rds"))$GSE89408$mat
stopifnot(identical(rownames(mat_old), rownames(mat_new)),
          identical(colnames(mat_old), colnames(mat_new)))

## 与 B6 run_ml11 逐字相同的 limma 设定
design <- model.matrix(~ 0 + y); colnames(design) <- c("Normal", "OA")
rank_of <- function(mat, tag) {
  tt <- topTable(eBayes(contrasts.fit(lmFit(mat, design), makeContrasts(OA - Normal, levels = design))),
                 coef = 1, number = Inf, sort.by = "none")
  de_ord <- rownames(tt)[order(-abs(tt$logFC))]
  pnl <- intersect(PANEL, rownames(mat))
  rk  <- match(pnl, de_ord)                       # 旧脚本 L154 的原式
  data.frame(matrix = tag, gene = pnl, rank = rk, logFC = tt[pnl, "logFC"],
             adjP = tt[pnl, "adj.P.Val"], n_genes_ranked = nrow(tt), row.names = NULL)
}
res <- rbind(rank_of(mat_new, "corrected"), rank_of(mat_old, "old_colSums"))
res$range <- ""

## 对照：旧矩阵是否复现旧值 4,535–16,154
o <- res[res$matrix == "old_colSums", ]
cat(sprintf("方法同一性检验（旧矩阵应复现旧值 4,535–16,154）：实得 %d–%d\n",
            min(o$rank), max(o$rank)))
cat(sprintf("  分母（参与排序的基因数）：校正矩阵 %d；旧矩阵 %d\n",
            res$n_genes_ranked[res$matrix == "corrected"][1],
            res$n_genes_ranked[res$matrix == "old_colSums"][1]))
n <- res[res$matrix == "corrected", ]
cat(sprintf("S2a 待填格（校正矩阵，%d 个基因中）：%d–%d\n", n$n_genes_ranked[1], min(n$rank), max(n$rank)))
cat("\n逐基因：\n"); print(res, row.names = FALSE)

write.csv(res, file.path(RES, "F11下游_ml11_面板名次.csv"), row.names = FALSE, fileEncoding = "UTF-8")
