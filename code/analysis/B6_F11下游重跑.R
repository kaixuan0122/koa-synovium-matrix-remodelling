# B6 F11 下游重跑（审计发现 F11 的处置：用修正后的 GSE89408 矩阵重跑受影响的产物）
#
# 背景（见 B0 与审计报告 F11）：03_preprocess.R 对 GSE89408 用 `colSums` 当文库大小。
#   该数据集 11 个 rRNA/结构 RNA 符号占全部 reads 33.9%，且占比与分组系统性共变
#   （OA 中位 19.06% vs 正常 44.65%）。总和归一化把这一构成差异读成了表达差异，
#   使每个基因的 logFC 整体向正方向偏移约 +0.32~+0.42（面板基因实测）。
#   改用 DESeq2 中位数比 size factor 后偏移消失（B0）。
#
# 影响面（已逐个核对，不是猜的）：
#   受影响：
#     * 04_merge_batch.R —— PCA 图与 主分析_合并校正.rds
#     * 07_cibersort.R   —— 免疫解卷积（§3.3）
#     * 10_machine_learning.R / 11_machine_learning_by_cohort.R / 12_panel_transfer.R（§3.7）
#   不受影响（无需重跑，此处仅记录理由）：
#     * 05_de_meta.R —— GSE89408 走 DESeq2 raw count，B0 已实测剔除 rRNA 前后
#       logFC 相关 r=0.99996、面板基因符号全一致。故 DEG→方向投票→Fisher→富集/PPI 全链不变。
#     * 13_replication_per_cohort.R —— 只涉 3 个验证队列（均为微阵列）。
#     * 11 号脚本里 GSE89408 的**面板显著性**（稿件"校正 P=0.18–0.83"）来自该脚本内嵌的
#       limma；但同一数字在 B1 中由 DESeq2 独立给出且完全一致（0.1836–0.8254），
#       故该结论按 DESeq2 口径成立，不构成本次重跑的待改项。
#
# 用法：Rscript B6_F11下游重跑.R [mat] [pca] [cibersort] [ml10] [ml11] [ml12]
#       不给参数则全跑。重的部分（ml10/ml11）建议后台单跑并各自留存日志。
#
# 输出：03_分析/logs/B6_<部分>.log（由调用方重定向）
#       03_分析/outputs/GSE89408_校正矩阵_cpm尺度.rds
#       03_分析/outputs/F11下游_<部分>_新旧对照.csv
suppressPackageStartupMessages({ library(Biobase); library(limma) })

BASE <- "/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/01.材料/01.博士材料/10.毕业/01.小论文/04.立项/滑膜基质刚度生成与刚度感知 0912"
P    <- file.path(BASE, "04_分析数据/processed")
RAW  <- file.path(BASE, "04_分析数据/raw")
OUT  <- file.path(BASE, "08_修订_2026-09-17/03_分析/outputs")
RES  <- file.path(BASE, "08_修订_2026-09-17/02_结果表")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

PANEL <- c("LOXL1", "ITGA5", "SDC4", "SERPINE1", "SAV1", "VGLL4")

ARGS  <- commandArgs(trailingOnly = TRUE)
PARTS <- if (length(ARGS)) ARGS else c("mat", "pca", "cibersort", "ml10", "ml11", "ml12")
cat(sprintf("本次运行部分：%s\n", paste(PARTS, collapse = ", ")))

## ============================================================
## 公共：把 size factor 版换算到 CPM 尺度
## ------------------------------------------------------------
## 为什么需要换算：下游脚本（10 号的 `CPM-log2 > 1` 预过滤、CIBERSORT 的 QN 分支）
## 都按 CPM 尺度设定阈值。size factor 版是 log2(count/sf + 1)，sf 以**全体基因的
## 几何均值**为基准，故它不是组成型归一化，数值落在"中位文库大小"的计数尺度上。
## 直接套 CPM 阈值会因量纲不同而几乎不删基因，"候选基因数"随之失去可比性。
## 换算：乘常数 k = 1e6 / median(列和)，使中位样本的取值与 CPM 同尺度。该常数是全局
## 的、不随样本变，故不会把校正过的构成差异又引回来。
## ============================================================
g89 <- readRDS(file.path(OUT, "GSE89408_基因水平_双版本.rds"))
size_mat <- g89$size

cnt <- read.delim(gzfile(file.path(RAW, "GSE89408/GSE89408_GEO_count_matrix_rename.txt.gz")),
                  check.names = FALSE, stringsAsFactors = FALSE)
SYM <- sub(":.*", "", cnt[, 1])
M   <- as.matrix(cnt[, -1]); colnames(M) <- sub("^X", "", colnames(M))
M   <- M[, c(grep("^OA_tissue_", colnames(M), value = TRUE),
             grep("^normal_tissue_", colnames(M), value = TRUE))]
Ms  <- aggregate(M, by = list(SYM), FUN = sum); rownames(Ms) <- Ms[, 1]
Ms  <- as.matrix(Ms[, -1, drop = FALSE])
K   <- 1e6 / median(colSums(Ms))
cat(sprintf("CPM 尺度换算常数 k = %.4f（= 1e6 / 中位列和 %.0f）\n", K, median(colSums(Ms))))

## 已归一化后又按符号取均值，故换算只能在聚合后的 log 矩阵上做：log2(count/sf*k+1)
## = log2((2^y - 1)*k + 1)。这一步与 B0 的 size 版一一对应，不引入新近似。
mat89408 <- log2((2^size_mat - 1) * K + 1)
stopifnot(identical(dim(mat89408), dim(size_mat)),
          max(abs(mat89408 - size_mat)) > 0)     # 换算确实生效
saveRDS(list(mat = mat89408, grp = g89$grp, K = K),
        file.path(OUT, "GSE89408_校正矩阵_cpm尺度.rds"))

## 旧矩阵（03_preprocess 的 colSums 版）——用于新旧对照
old_main <- readRDS(file.path(P, "主分析_genelevel.rds"))
old89408 <- old_main$GSE89408
stopifnot(identical(rownames(old89408$mat), rownames(mat89408)),
          identical(colnames(old89408$mat), colnames(mat89408)))

## ============================================================
## 部分 pca：重跑 04_merge_batch.R 的 PCA（只换 GSE89408 那一支）
## ============================================================
run_pca <- function() {
  mats <- list(GSE89408 = mat89408, GSE206848 = old_main$GSE206848$mat,
               GSE1919   = old_main$GSE1919$mat)
  grps <- list(GSE89408 = g89$grp, GSE206848 = old_main$GSE206848$grp,
               GSE1919   = old_main$GSE1919$grp)
  asgrp <- function(g) if (is.factor(g)) as.character(g) else g
  common <- Reduce(intersect, lapply(mats, rownames))
  merged <- do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE]))
  groups  <- unlist(lapply(grps, asgrp))
  batch   <- rep(names(mats), sapply(mats, ncol))

  pca_plot <- function(m, tag) {
    pc <- prcomp(t(m), center = TRUE, scale. = TRUE)
    pve <- round(100 * summary(pc)$importance[2, 1:2], 1)
    pdf(file.path(OUT, paste0("F11重跑_PCA_", tag, ".pdf")), width = 7, height = 5.5)
    plot(pc$x[, 1], pc$x[, 2], col = as.integer(factor(batch)),
         pch = ifelse(groups == "OA", 17, 1),
         xlab = paste0("PC1 (", pve[1], "%)"), ylab = paste0("PC2 (", pve[2], "%)"),
         main = paste0("PCA ", tag))
    legend("topright", legend = unique(batch), col = 1:length(unique(batch)),
           pch = 15, cex = 0.7)
    dev.off()
    cat(sprintf("  PCA %s: PC1=%.1f%% PC2=%.1f%%\n", tag, pve[1], pve[2]))
    pve
  }
  cat("PCA（GSE89408 用 B0 校正矩阵）:\n")
  a <- pca_plot(merged, "主分析_校正前")
  merged_cor <- removeBatchEffect(merged, batch = batch, design = model.matrix(~groups))
  b <- pca_plot(merged_cor, "主分析_校正后")
  saveRDS(list(expr = merged_cor, groups = groups, batch = batch, common_genes = common),
          file.path(OUT, "F11重跑_主分析_合并校正.rds"))

  ## 与旧 PCA 对照：旧矩阵下同一套流程的 PVE 与样本坐标
  old_merged <- do.call(cbind, lapply(
    list(GSE89408 = old89408$mat, GSE206848 = old_main$GSE206848$mat,
         GSE1919 = old_main$GSE1919$mat), function(m) m[common, , drop = FALSE]))
  pc_old <- prcomp(t(old_merged), center = TRUE, scale. = TRUE)
  pve_old <- round(100 * summary(pc_old)$importance[2, 1:2], 1)
  cat(sprintf("  对照（旧矩阵）: PC1=%.1f%% PC2=%.1f%%\n", pve_old[1], pve_old[2]))
  cat(sprintf("  → PCA 是质控图，不进正文主图；此处只记录重跑结果，不改稿件结论\n"))
  data.frame(版本 = c("旧(colSums)", "新(size factor)"),
             PC1_PVE = c(pve_old[1], a[1]), PC2_PVE = c(pve_old[2], a[2]))
}

## ============================================================
## 部分 cibersort：重跑 07_cibersort.R（三队列全跑，含 GSE89408 校正矩阵）
## ============================================================
run_cibersort <- function() {
  suppressPackageStartupMessages({ library(e1071); library(preprocessCore) })
  lm22 <- "/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/01.材料/01.博士材料/03.小论文选题/生信分析/rnaseq-workflow-agent/data/LM22.txt"
  cib  <- "/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/01.材料/01.博士材料/03.小论文选题/生信分析/rnaseq-workflow-agent/data/Cibersort.R"
  source(cib)
  ## 三队列全跑：GSE206848/GSE1919 的输入矩阵未变，全跑可顺带验证 CIBERSORT 的
  ## 分数是否可复现（若这两个队列的新旧分数不同，说明流程本身有随机性，
  ## 那么 GSE89408 的新旧差异就不能全部归因于归一化修正）。
  queues <- list(
    GSE89408  = list(mat = mat89408,           grp = g89$grp,               QN = FALSE),
    GSE206848 = list(mat = old_main$GSE206848$mat, grp = old_main$GSE206848$grp, QN = TRUE),
    GSE1919   = list(mat = old_main$GSE1919$mat,   grp = old_main$GSE1919$grp,   QN = TRUE)
  )
  tmp <- file.path(OUT, "B6_ciber_tmp"); dir.create(tmp, showWarnings = FALSE)
  all_res <- list()
  for (nm in names(queues)) {
    q <- queues[[nm]]
    mf <- file.path(tmp, paste0(nm, "_mixture.txt"))
    write.table(q$mat, mf, sep = "\t", quote = FALSE, col.names = NA)
    owd <- getwd(); setwd(tmp)
    res <- CIBERSORT(lm22, mf, perm = 100, QN = q$QN)
    setwd(owd)
    res <- as.data.frame(res); res$sample <- rownames(res)
    res$group  <- q$grp[match(rownames(res), colnames(q$mat))]
    res$dataset <- nm
    all_res[[nm]] <- res
    cat(sprintf("  %s CIBERSORT 完成: %d 样本\n", nm, nrow(res)))
  }
  cell_types <- setdiff(colnames(all_res[[1]]), c("sample", "group", "dataset",
                                                  "P-value", "Correlation", "RMSE"))
  ## 1) 三队列全跑：新旧逐细胞类型对照（未变队列应完全一致）
  old_tab <- list()
  for (nm in names(all_res)) {
    f <- file.path(P, paste0("CIBERSORT_", nm, "_results.txt"))
    if (!file.exists(f)) next
    o <- read.table(f, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    m <- intersect(rownames(o), all_res[[nm]]$sample)
    d <- abs(as.matrix(o[m, cell_types]) - as.matrix(all_res[[nm]][match(m, all_res[[nm]]$sample), cell_types]))
    old_tab[[nm]] <- data.frame(dataset = nm, n_共同样本 = length(m),
                                max_abs_差 = max(d), 中位_abs_差 = median(d))
  }
  OLD <- do.call(rbind, old_tab)
  cat("\n  【可复现性核对】新旧 CIBERSORT 分数逐样本最大绝对差：\n"); print(OLD, row.names = FALSE)
  cat("  → 未变队列的差若非 0，则该队列的新旧差异不能全归因于归一化修正\n")

  ## 2) 合并三队列后的 OA vs Normal（与 07 号脚本同口径：合并样本 Wilcoxon）
  merged <- do.call(rbind, lapply(all_res, function(d) d))
  newdiff <- do.call(rbind, lapply(cell_types, function(ct) {
    oa <- merged[[ct]][merged$group == "OA"]; no <- merged[[ct]][merged$group == "Normal"]
    w <- wilcox.test(oa, no)
    data.frame(cell = ct, mean_OA = mean(oa), mean_Normal = mean(no), p = w$p.value)
  }))
  newdiff$p.adj <- p.adjust(newdiff$p, method = "BH")
  olddiff <- read.csv(file.path(P, "CIBERSORT_免疫细胞差异.csv"), stringsAsFactors = FALSE)
  cmp <- merge(olddiff, newdiff, by = "cell", suffixes = c(".old", ".new"))
  cmp <- cmp[order(cmp$p.new), ]
  write.csv(cmp, file.path(RES, "F11下游_cibersort_新旧对照.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(merged, file.path(OUT, "F11重跑_CIBERSORT_三队列_合并.csv"), row.names = FALSE)

  cat("\n  【结论核对】免疫细胞 OA vs Normal（按新 P 排序，前 6）：\n")
  print(head(cmp[, c("cell", "p.old", "p.adj.old", "p.new", "p.adj.new")], 6), row.names = FALSE)
  cat(sprintf("\n  BH 后显著（<0.05）的细胞类型数：旧 %d 个，新 %d 个\n",
              sum(olddiff$p.adj < 0.05), sum(newdiff$p.adj < 0.05)))
  ## 稿件 §3.3 点名的两个趋势
  for (ct in c("Macrophages M2", "T cells CD4 memory resting")) {
    r <- cmp[cmp$cell == ct, ]
    cat(sprintf("  %-30s P %.4f -> %.4f ；BH 后 %.3f -> %.3f\n",
                ct, r$p.old, r$p.new, r$p.adj.old, r$p.adj.new))
  }
  cmp
}

## ============================================================
## 部分 ml10：重跑 10_machine_learning.R（只换 GSE89408 矩阵）
## ============================================================
run_ml10 <- function() {
  suppressPackageStartupMessages({
    library(glmnet); library(e1071); library(randomForest); library(ranger)
    library(gbm); library(rpart); library(pROC); library(xgboost)
  })
  set.seed(20260913)   # 与 10 号脚本同种子；该脚本只处理 GSE89408，随机流可直接对齐
  OUTML <- file.path(OUT, "F11重跑_ML10"); dir.create(OUTML, showWarnings = FALSE)
  TOPK  <- 20
  step  <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

  d <- list(mat = mat89408, grp = g89$grp)
  X <- t(d$mat)
  y <- factor(ifelse(d$grp == "OA", "OA", "Normal"), levels = c("Normal", "OA"))
  cat(sprintf("GSE89408(校正): %d 样本 (%d OA / %d Normal) × %d 基因\n",
              nrow(X), sum(y == "OA"), sum(y == "Normal"), ncol(X)))
  ## 标签无关过滤：照抄 10 号脚本原条件。因已换算到 CPM 尺度，阈值语义与原版一致。
  keep <- colMeans(X > 1) >= 0.5 & apply(X, 2, sd) > 0
  n_raw <- ncol(X)
  X <- X[, keep, drop = FALSE]
  cat(sprintf("表达过滤后: %d 基因（原版为 24,894 -> 14,510）；面板基因均在候选集: %s\n",
              ncol(X), all(PANEL %in% colnames(X))))
  Xs <- scale(X); Xs[is.na(Xs)] <- 0

  sel <- list()
  step("LASSO")
  cv1 <- cv.glmnet(Xs, y, family = "binomial", alpha = 1, nfolds = 5)
  cf1 <- as.numeric(coef(cv1, s = "lambda.1se"))[-1]
  sel$LASSO <- names(sort(abs(cf1[cf1 != 0]), decreasing = TRUE))[1:min(TOPK, sum(cf1 != 0))]
  step("SVM-RFE")
  Xr <- Xs; kept <- colnames(Xs)
  while (length(kept) > TOPK) {
    m <- svm(Xr[, kept, drop = FALSE], y, kernel = "linear", scale = FALSE, type = "C-classification")
    w <- as.numeric(t(m$coefs) %*% m$SV); names(w) <- kept
    kept <- setdiff(kept, names(sort(abs(w)))[1:max(1, floor(0.1 * length(kept)))])
  }
  sel$`SVM-RFE` <- kept
  step("Random forest")
  rf <- randomForest(Xs, y, ntree = 1000, importance = TRUE)
  imp_rf <- sort(rf$importance[, "MeanDecreaseGini"], decreasing = TRUE)
  sel$`Random forest` <- names(imp_rf)[1:TOPK]
  step("ExtraTrees")
  et <- ranger(x = Xs, y = y, splitrule = "extratrees", importance = "impurity", num.trees = 1000)
  imp_et <- sort(et$variable.importance, decreasing = TRUE)
  sel$ExtraTrees <- names(imp_et)[1:TOPK]
  step("Gradient boosting")
  gb <- gbm.fit(Xs, as.numeric(y) - 1, distribution = "bernoulli", n.trees = 2000,
                shrinkage = 0.01, interaction.depth = 3, n.minobsinnode = 2, verbose = FALSE)
  imp_gb <- sort(summary(gb, plotit = FALSE)[, "rel.inf"], decreasing = TRUE)
  sel$`Gradient boosting` <- names(imp_gb)[1:TOPK]
  step("AdaBoost")
  df_ab <- data.frame(Xs, check.names = FALSE)
  colnames(df_ab) <- paste0("g", seq_len(ncol(Xs))); df_ab$OA <- y
  w <- rep(1 / nrow(Xs), nrow(Xs)); imp_ab <- setNames(numeric(ncol(Xs)), colnames(Xs)); n_skip <- 0
  for (i in 1:100) {
    o <- local({
      tr <- rpart(OA ~ ., data = df_ab, method = "class", weights = w, model = FALSE,
                  control = rpart.control(maxdepth = 4, minsplit = 5, cp = 1e-6, xval = 0))
      list(pred = as.character(predict(tr, df_ab, type = "class")), imp = tr$variable.importance)
    })
    err <- sum(w * (o$pred != as.character(y))) / sum(w)
    if (err >= 0.5) { n_skip <- n_skip + 1; next }
    a <- 0.5 * log((1 - err) / max(err, 1e-10))
    w <- w * exp(a * (o$pred != as.character(y))); w <- w / sum(w)
    idx <- as.integer(sub("^g", "", names(o$imp))); imp_ab[idx] <- imp_ab[idx] + a * as.numeric(o$imp)
  }
  cat(sprintf("  AdaBoost 100 轮中 %d 轮被跳过\n", n_skip))
  sel$AdaBoost <- names(sort(imp_ab, decreasing = TRUE))[1:TOPK]
  step("L1 logistic")
  cv7 <- cv.glmnet(Xs, y, family = "binomial", alpha = 1, nfolds = 5)
  cf7 <- as.numeric(coef(cv7, s = "lambda.min"))[-1]
  sel$`L1 logistic` <- names(sort(abs(cf7[cf7 != 0]), decreasing = TRUE))[1:min(TOPK, sum(cf7 != 0))]

  allg <- unique(unlist(sel))
  cnt <- data.frame(Gene = allg,
                    n_methods = sapply(allg, function(g) sum(sapply(sel, function(s) g %in% s))))
  cnt <- cnt[order(-cnt$n_methods, cnt$Gene), ]
  cat(sprintf("\n面板基因被任何算法选中的个数: %d / 6\n",
              sum(PANEL %in% cnt$Gene)))

  ## SHAP
  step("SHAP")
  dm  <- xgb.DMatrix(Xs, label = as.numeric(y) - 1)
  bst <- xgb.train(params = list(objective = "binary:logistic", max_depth = 3, learning_rate = 0.05),
                   data = dm, nrounds = 200, verbose = 0)
  shap <- predict(bst, dm, predcontrib = TRUE)
  shap <- shap[, !colnames(shap) %in% c("BIAS", "(Intercept)"), drop = FALSE]
  imp_shap <- sort(colMeans(abs(shap)), decreasing = TRUE)
  nz <- names(imp_shap)[imp_shap > 0]
  rank_tab <- data.frame(Gene = colnames(Xs), stringsAsFactors = FALSE)
  rank_tab$RF_used <- rank_tab$Gene %in% names(imp_rf)[imp_rf > 0]
  rank_tab$SHAP_used <- rank_tab$Gene %in% nz
  rank_tab$mean_abs_SHAP <- round(as.numeric(imp_shap[rank_tab$Gene]), 6)

  ## AUC
  step("AUC")
  auc1 <- sapply(PANEL, function(g) as.numeric(auc(roc(y, X[, g], quiet = TRUE))))
  df <- data.frame(Xs[, PANEL, drop = FALSE], y = y)
  cvpred <- numeric(nrow(df)); folds <- sample(rep(1:5, length.out = nrow(df)))
  for (k in 1:5) {
    tr <- folds != k; te <- folds == k
    m <- glm(y ~ ., data = df[tr, ], family = "binomial")
    cvpred[te] <- predict(m, df[te, ], type = "response")
  }
  auc_cv <- as.numeric(auc(roc(y, cvpred, quiet = TRUE)))
  v <- readRDS(file.path(P, "验证层_合并校正.rds"))
  gv <- intersect(PANEL, rownames(v$expr))
  yv <- factor(ifelse(v$groups == "OA", "OA", "Normal"), levels = c("Normal", "OA"))
  zscore <- function(m) scale(m)
  tr <- data.frame(zscore(X[, gv, drop = FALSE]), y = y)
  mdl <- glm(y ~ ., data = tr, family = "binomial")
  te <- data.frame(zscore(t(v$expr[gv, , drop = FALSE])), y = yv)
  auc_ext <- as.numeric(auc(roc(yv, predict(mdl, te, type = "response"), quiet = TRUE)))

  ## 新旧对照：原版数值一律**从旧产物读取**，不手写。
  ## 手写数字会与产物脱节（本稿已有前例：稿件自报的样本量与预处理脚本对不上）。
  MV  <- file.path(BASE, "04_分析数据/合并验证")
  old_auc1 <- read.csv(file.path(MV, "ML_单基因AUC.csv"), stringsAsFactors = FALSE)
  old_j    <- read.csv(file.path(MV, "ML_联合模型AUC.csv"), stringsAsFactors = FALSE)
  old_rank <- read.csv(file.path(MV, "ML_重要性排名.csv"), stringsAsFactors = FALSE)
  old_cnt  <- read.csv(file.path(MV, "ML_基因方法数.csv"), stringsAsFactors = FALSE)
  ## 原版候选基因数：对**旧矩阵**施加同一过滤条件复算，而不是引用稿件里的手写数
  Xo <- t(old89408$mat)
  old_ncand <- sum(colMeans(Xo > 1) >= 0.5 & apply(Xo, 2, sd) > 0)
  old_cv  <- old_j$value[old_j$metric == "AUC_5fold_CV_GSE89408"]
  old_ext <- old_j$value[old_j$metric == "AUC_external_3cohorts"]
  cmp1 <- merge(old_auc1, data.frame(Gene = PANEL, AUC_new = round(auc1, 3)), by = "Gene")
  cat("\n=== F11 重跑 vs 原版（10 号脚本，GSE89408）===\n")
  cat(sprintf("  候选基因数            %d -> %d\n", old_ncand, ncol(X)))
  cat(sprintf("  面板被任何算法选中    %d/6 -> %d/6\n",
              sum(PANEL %in% old_cnt$Gene), sum(PANEL %in% cnt$Gene)))
  cat(sprintf("  单基因 AUC 范围       %.3f-%.3f -> %.3f-%.3f\n",
              min(old_auc1$AUC), max(old_auc1$AUC), min(auc1), max(auc1)))
  cat(sprintf("  联合模型内部 CV AUC   %.3f -> %.3f\n", old_cv, auc_cv))
  cat(sprintf("  外部（3 验证队列合并）%.3f -> %.3f\n", old_ext, auc_ext))
  old_rf   <- sum(old_rank$RF_used); old_shap <- sum(old_rank$SHAP_used)
  cat(sprintf("  树模型用到的特征数    RF %d -> %d ; SHAP %d -> %d\n",
              old_rf, sum(imp_rf > 0), old_shap, length(nz)))
  cat(sprintf("  面板基因 mean|SHAP| 全为 0：原版 %s -> 重跑 %s\n",
              all(old_rank$mean_abs_SHAP[old_rank$Gene %in% PANEL] == 0),
              all(rank_tab$mean_abs_SHAP[rank_tab$Gene %in% PANEL] == 0)))
  print(cmp1, row.names = FALSE)

  out <- data.frame(指标 = c("候选基因数", "面板被选中数", "单基因AUC下限", "单基因AUC上限",
                             "联合模型CV_AUC", "外部验证AUC", "RF使用特征数", "SHAP使用特征数"),
                    原版 = c(old_ncand, sum(PANEL %in% old_cnt$Gene), min(old_auc1$AUC),
                             max(old_auc1$AUC), old_cv, old_ext, old_rf, old_shap),
                    重跑 = c(ncol(X), sum(PANEL %in% cnt$Gene), round(min(auc1), 3),
                             round(max(auc1), 3), round(auc_cv, 3), round(auc_ext, 3),
                             sum(imp_rf > 0), length(nz)))
  write.csv(out, file.path(RES, "F11下游_ml10_新旧对照.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(cnt, file.path(OUTML, "ML_基因方法数.csv"), row.names = FALSE)
  write.csv(rank_tab, file.path(OUTML, "ML_重要性排名.csv"), row.names = FALSE)
  out
}

## ============================================================
## 部分 ml11：重跑 11_machine_learning_by_cohort.R 的 GSE89408 一行
##   该脚本的种子按 match(nm, ORD_ALL) 逐队列确定，故单跑 GSE89408 与全跑结果一致
## ============================================================
run_ml11 <- function() {
  suppressPackageStartupMessages({
    library(glmnet); library(e1071); library(randomForest); library(ranger)
    library(gbm); library(rpart); library(pROC); library(limma)
  })
  ORD_ALL <- c("GSE55235", "GSE55457", "GSE206848", "GSE82107", "GSE89408")
  nm <- "GSE89408"
  set.seed(20260913 + match(nm, ORD_ALL))
  TOPK <- 20
  step <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

  run_ml <- function(Xs, y) {
    sel <- list()
    cv1 <- cv.glmnet(Xs, y, family = "binomial", alpha = 1, nfolds = 5)
    cf1 <- as.numeric(coef(cv1, s = "lambda.1se"))[-1]
    sel$LASSO <- names(sort(abs(cf1[cf1 != 0]), decreasing = TRUE))[1:min(TOPK, sum(cf1 != 0))]
    Xr <- Xs; kept <- colnames(Xs)
    while (length(kept) > TOPK) {
      m <- svm(Xr[, kept, drop = FALSE], y, kernel = "linear", scale = FALSE, type = "C-classification")
      w <- as.numeric(t(m$coefs) %*% m$SV); names(w) <- kept
      kept <- setdiff(kept, names(sort(abs(w)))[1:max(1, floor(0.1 * length(kept)))])
    }
    sel$`SVM-RFE` <- kept
    rf <- randomForest(Xs, y, ntree = 1000, importance = TRUE)
    sel$`Random forest` <- names(sort(rf$importance[, "MeanDecreaseGini"], decreasing = TRUE))[1:TOPK]
    et <- ranger(x = Xs, y = y, splitrule = "extratrees", importance = "impurity", num.trees = 1000)
    sel$ExtraTrees <- names(sort(et$variable.importance, decreasing = TRUE))[1:TOPK]
    gb <- gbm.fit(Xs, as.numeric(y) - 1, distribution = "bernoulli", n.trees = 2000,
                  shrinkage = 0.01, interaction.depth = 3, n.minobsinnode = 2, verbose = FALSE)
    sel$`Gradient boosting` <- names(sort(summary(gb, plotit = FALSE)[, "rel.inf"], decreasing = TRUE))[1:TOPK]
    df_ab <- data.frame(Xs, check.names = FALSE)
    colnames(df_ab) <- paste0("g", seq_len(ncol(Xs))); df_ab$OA <- y
    w <- rep(1 / nrow(Xs), nrow(Xs)); imp_ab <- setNames(numeric(ncol(Xs)), colnames(Xs))
    for (i in 1:100) {
      o <- local({
        tr <- rpart(OA ~ ., data = df_ab, method = "class", weights = w, model = FALSE,
                    control = rpart.control(maxdepth = 4, minsplit = 5, cp = 1e-6, xval = 0))
        list(pred = as.character(predict(tr, df_ab, type = "class")), imp = tr$variable.importance)
      })
      err <- sum(w * (o$pred != as.character(y))) / sum(w)
      if (err >= 0.5) next
      a <- 0.5 * log((1 - err) / max(err, 1e-10))
      w <- w * exp(a * (o$pred != as.character(y))); w <- w / sum(w)
      idx <- as.integer(sub("^g", "", names(o$imp))); imp_ab[idx] <- imp_ab[idx] + a * as.numeric(o$imp)
    }
    sel$AdaBoost <- names(sort(imp_ab, decreasing = TRUE))[1:TOPK]
    cv7 <- cv.glmnet(Xs, y, family = "binomial", alpha = 1, nfolds = 5)
    cf7 <- as.numeric(coef(cv7, s = "lambda.min"))[-1]
    sel$`L1 logistic` <- names(sort(abs(cf7[cf7 != 0]), decreasing = TRUE))[1:min(TOPK, sum(cf7 != 0))]
    sel
  }

  step("===== ", nm, " (校正矩阵) =====")
  mat <- mat89408
  y <- factor(ifelse(g89$grp == "OA", "OA", "Normal"), levels = c("Normal", "OA"))
  n <- ncol(mat)
  ## 面板显著数：原脚本用本队列矩阵跑 limma。此处**同时**给出校正矩阵与 DESeq2 两个口径，
  ## 因稿件"校正 P=0.18-0.83"实际取自 DESeq2（B1），不是本脚本的 limma。
  design <- model.matrix(~ 0 + y); colnames(design) <- c("Normal", "OA")
  tt <- topTable(eBayes(contrasts.fit(lmFit(mat, design), makeContrasts(OA - Normal, levels = design))),
                 coef = 1, number = Inf, sort.by = "none")
  panel <- intersect(PANEL, rownames(mat))
  n_sig <- sum(tt[panel, "adj.P.Val"] < 0.05, na.rm = TRUE)
  cat(sprintf("  面板基因在本队列 limma adjP（校正矩阵）: %s\n",
              paste(sprintf("%s=%.3f", panel, tt[panel, "adj.P.Val"]), collapse = ", ")))
  old_tt <- topTable(eBayes(contrasts.fit(
    lmFit(old89408$mat, design), makeContrasts(OA - Normal, levels = design))),
    coef = 1, number = Inf, sort.by = "none")
  cat(sprintf("  面板基因在本队列 limma adjP（原 colSums 矩阵）: %s\n",
              paste(sprintf("%s=%.3f", panel, old_tt[panel, "adj.P.Val"]), collapse = ", ")))

  gm <- apply(mat, 1, median)
  keep <- gm > quantile(gm, 0.25) & apply(mat, 1, sd) > 0
  mat <- mat[keep, , drop = FALSE]
  X <- t(mat); panel <- intersect(panel, colnames(X))
  Xs <- scale(X); Xs[is.na(Xs)] <- 0
  old_row <- read.csv(file.path(BASE, "04_分析数据/合并验证/ML_分队列_汇总.csv"),
                      stringsAsFactors = FALSE)
  old_row <- old_row[old_row$队列 == nm, ]
  cat(sprintf("  n=%d；候选基因 %d 个（原版 %d）\n", n, ncol(Xs), old_row$候选基因))

  sel <- run_ml(Xs, y)
  allg <- unique(unlist(sel))
  cnt <- data.frame(Gene = allg,
                    n_methods = sapply(allg, function(g) sum(sapply(sel, function(s) g %in% s))))
  cnt <- cnt[order(-cnt$n_methods, cnt$Gene), ]
  hit <- cnt[cnt$Gene %in% panel, ]
  cat(sprintf("  面板基因被选中: %s\n",
              if (nrow(hit) == 0) "无 (0/6)" else paste(sprintf("%s(%d)", hit$Gene, hit$n_methods), collapse = ", ")))

  de_ord <- rownames(tt)[order(-abs(tt$logFC))]
  sig_de <- rownames(tt)[!is.na(tt$adj.P.Val) & tt$adj.P.Val < 0.05]
  top50 <- head(intersect(de_ord, sig_de), 50)
  ov50 <- length(intersect(unique(unlist(sel)), top50))
  cat(sprintf("  阳性对照: 显著 DE %d 个（原版 %d）；前 50 强被回收 %d 个（原版 %d）\n",
              length(sig_de), old_row$显著DE数, ov50, old_row$前50强DE被回收数))

  auc1 <- sapply(panel, function(g) as.numeric(auc(roc(y, X[, g], quiet = TRUE))))
  cvp <- numeric(n); folds <- sample(rep(1:5, length.out = n))
  Pf <- Xs[, panel, drop = FALSE]; colnames(Pf) <- paste0("p", seq_along(panel))
  for (k in 1:5) {
    tr <- folds != k; te <- folds == k
    mm <- glm(y ~ ., data = data.frame(Pf[tr, , drop = FALSE], y = y[tr]), family = "binomial")
    cvp[te] <- predict(mm, data.frame(Pf[te, , drop = FALSE]), type = "response")
  }
  auc_cv <- as.numeric(auc(roc(y, cvp, quiet = TRUE)))
  cat(sprintf("  单基因 AUC %.3f-%.3f（原版 %.3f-%.3f）；CV AUC=%.3f（原版 %.3f）\n",
              min(auc1), max(auc1), old_row$单基因AUC下限, old_row$单基因AUC上限,
              auc_cv, old_row$联合模型CV_AUC))

  out <- data.frame(
    队列 = nm, n = n, OA = sum(y == "OA"), 正常 = sum(y == "Normal"),
    候选基因 = ncol(Xs), 面板显著数 = sprintf("%d/%d", n_sig, length(panel)),
    面板被选中数 = nrow(hit),
    被选中面板基因 = if (nrow(hit) == 0) "无" else paste(hit$Gene, collapse = ","),
    显著DE数 = length(sig_de), 前50强DE被回收数 = ov50,
    单基因AUC下限 = round(min(auc1), 3), 单基因AUC上限 = round(max(auc1), 3),
    联合模型CV_AUC = round(auc_cv, 3))
  write.csv(out, file.path(RES, "F11下游_ml11_GSE89408行.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(cnt, file.path(OUT, sprintf("F11重跑_ML_分队列_%s.csv", nm)), row.names = FALSE)
  cat("\n=== F11 重跑 vs 原版（11 号脚本 GSE89408 行）===\n")
  print(old_row, row.names = FALSE)
  print(out, row.names = FALSE)
  out
}

## ============================================================
## 部分 ml12：重跑 12_panel_transfer.R（五队列全跑，含 GSE89408 校正矩阵）
## ============================================================
run_ml12 <- function() {
  suppressPackageStartupMessages({ library(pROC) })
  NREP <- 10; NPERM <- 200
  A <- list(GSE89408 = list(mat = mat89408, grp = g89$grp),
            GSE206848 = list(mat = old_main$GSE206848$mat, grp = old_main$GSE206848$grp),
            GSE1919   = list(mat = old_main$GSE1919$mat,   grp = old_main$GSE1919$grp))
  V <- readRDS(file.path(P, "验证层_genelevel.rds"))
  ALL <- c(A, V)
  ORD <- c("GSE55235", "GSE55457", "GSE206848", "GSE82107", "GSE89408")
  ORD <- ORD[ORD %in% names(ALL)]
  get <- function(nm) {
    mat <- ALL[[nm]]$mat
    y <- factor(ifelse(ALL[[nm]]$grp == "OA", "OA", "Normal"), levels = c("Normal", "OA"))
    list(mat = mat, y = y, genes = PANEL[PANEL %in% rownames(mat)])
  }
  cv_once <- function(X, y) {
    n <- nrow(X); folds <- sample(rep(1:5, length.out = n)); p <- numeric(n)
    for (k in 1:5) {
      tr <- folds != k; te <- folds == k
      mu <- colMeans(X[tr, , drop = FALSE]); sg <- apply(X[tr, , drop = FALSE], 2, sd); sg[sg == 0] <- 1
      a <- scale(X[tr, , drop = FALSE], center = mu, scale = sg)
      b <- scale(X[te, , drop = FALSE], center = mu, scale = sg)
      colnames(a) <- colnames(b) <- paste0("p", seq_len(ncol(X)))
      d <- as.data.frame(a); d$y <- y[tr]
      m <- suppressWarnings(glm(y ~ ., data = d, family = "binomial"))
      p[te] <- predict(m, newdata = as.data.frame(b), type = "response")
    }
    as.numeric(auc(roc(y, p, quiet = TRUE)))
  }
  set.seed(20260913)
  cal <- list()
  for (nm in ORD) {
    d <- get(nm); Xp <- t(d$mat[d$genes, , drop = FALSE])
    obs  <- replicate(NREP, cv_once(Xp, d$y))
    pool <- setdiff(rownames(d$mat), d$genes)
    null <- replicate(NPERM, {
      g <- sample(pool, length(d$genes)); cv_once(t(d$mat[g, , drop = FALSE]), d$y)
    })
    pval <- (sum(null >= mean(obs)) + 1) / (NPERM + 1)
    cat(sprintf("  %-11s n=%2d  面板 CV AUC=%.3f [%.3f-%.3f]  零分布中位=%.3f 95分位=%.3f  P=%.3f\n",
                nm, nrow(Xp), mean(obs), min(obs), max(obs),
                median(null), quantile(null, 0.95), pval))
    cal[[nm]] <- data.frame(队列 = nm, n = nrow(Xp), 面板CV_AUC = round(mean(obs), 3),
                            CV下限 = round(min(obs), 3), CV上限 = round(max(obs), 3),
                            零分布中位 = round(median(null), 3),
                            零分布95分位 = round(quantile(null, 0.95), 3), 经验P = round(pval, 3))
  }
  CAL <- do.call(rbind, cal)
  old_cal <- read.csv(file.path(BASE, "04_分析数据/合并验证/ML_判别力定标.csv"), stringsAsFactors = FALSE)
  cmp <- merge(old_cal, CAL, by = "队列", suffixes = c(".old", ".new"))
  write.csv(cmp, file.path(RES, "F11下游_ml12_定标新旧对照.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  cat("\n【定标新旧对照】\n"); print(cmp[, c("队列", "面板CV_AUC.old", "面板CV_AUC.new",
                                               "零分布中位.old", "零分布中位.new",
                                               "经验P.old", "经验P.new")], row.names = FALSE)

  M <- matrix(NA_real_, length(ORD), length(ORD), dimnames = list(ORD, ORD))
  for (tr_nm in ORD) {
    dtr <- get(tr_nm); g <- dtr$genes
    Xtr <- t(dtr$mat[g, , drop = FALSE])
    mu <- colMeans(Xtr); sg <- apply(Xtr, 2, sd); sg[sg == 0] <- 1
    a <- scale(Xtr, center = mu, scale = sg); colnames(a) <- paste0("p", seq_along(g))
    d <- as.data.frame(a); d$y <- dtr$y
    m <- suppressWarnings(glm(y ~ ., data = d, family = "binomial"))
    for (te_nm in ORD) {
      dte <- get(te_nm)
      Xte <- scale(t(dte$mat[g, , drop = FALSE]), center = mu, scale = sg)
      colnames(Xte) <- paste0("p", seq_along(g))
      M[tr_nm, te_nm] <- suppressMessages(as.numeric(
        auc(roc(dte$y, predict(m, newdata = as.data.frame(Xte), type = "response"), quiet = TRUE))))
    }
  }
  write.csv(round(M, 3), file.path(OUT, "F11重跑_ML_迁移矩阵_AUC.csv"))
  old_M <- as.matrix(read.csv(file.path(BASE, "04_分析数据/合并验证/ML_迁移矩阵_AUC.csv"),
                              row.names = 1, check.names = FALSE))
  cm <- intersect(rownames(M), rownames(old_M))
  cat(sprintf("\n【迁移矩阵新旧对照】GSE89408 行/列\n"))
  print(data.frame(训练 = cm, GSE89408列_旧 = old_M[cm, "GSE89408"],
                   GSE89408列_新 = M[cm, "GSE89408"]), row.names = FALSE)
  print(data.frame(测试 = cm, GSE89408行_旧 = old_M["GSE89408", cm],
                   GSE89408行_新 = M["GSE89408", cm]), row.names = FALSE)
  cat(sprintf("\n对角均值 %.3f（原 %.3f）；非对角均值 %.3f（原 %.3f）\n",
              mean(diag(M)), mean(diag(old_M)),
              mean(M[upper.tri(M) | lower.tri(M)]), mean(old_M[upper.tri(old_M) | lower.tri(old_M)])))
  M
}

## ============================================================
res <- list()
if ("pca"       %in% PARTS) res$pca       <- run_pca()
if ("cibersort" %in% PARTS) res$cibersort <- run_cibersort()
if ("ml10"      %in% PARTS) res$ml10      <- run_ml10()
if ("ml11"      %in% PARTS) res$ml11      <- run_ml11()
if ("ml12"      %in% PARTS) res$ml12      <- run_ml12()
cat(sprintf("\nB6 完成：%s\n", paste(names(res), collapse = ", ")))
