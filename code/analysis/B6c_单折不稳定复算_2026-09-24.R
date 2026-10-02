# B6c：复算 GSE89408 单次 5 折 CV 的折间波动（S2b 第 4 条 / S2c 括注的替换来源）
#
# 为什么要跑：S2 有两处引用「单次 5 折 CV 在 n=50 下不稳定」的具体数：
#   S2b 第 4 条：同一数据同一流程（全队列标准化、单次折划分）三次运行 0.515 / 0.503 / 0.468
#   S2c 括注  ：改用全队列标准化，内部 CV AUC 会降至 0.468–0.515
# 审计发现这三/两个数**没有可复现的出处**：
#   * 12_panel_transfer.R L6 记的是 0.515/0.503/**0.578**，与 S2 的 0.468 打架；
#   * 该脚本 `set.seed(20260913)` 在顶层，反复运行不会给出不同的第三个值；
#   * 且该脚本的 cv_once 用的是**只有训练折**标准化，与 S2 写的「全队列标准化」不符。
# 三个数还都取自旧的 colSums 矩阵。故按规则 3/4，重建 S2 时不能沿用，改由本脚本复算。
#
# 本脚本对两个矩阵、两种标准化口径，各做 3 次固定种子的单次 5 折 CV：
#   whole  = 标准化参数取自**全队列**（含测试折，即泄漏版，对应 S2 原文所写的口径）
#   train  = 标准化参数只取自训练折（12 号脚本 cv_once 的实际口径）
# 3 个折划分用固定种子 1/2/3 生成，与任何其他脚本的 RNG 流无关，故完全可复现。
#
# 输入：03_分析/outputs/GSE89408_校正矩阵_cpm尺度.rds；04_分析数据/processed/主分析_genelevel.rds
# 用法：Rscript B6c_单折不稳定复算_2026-09-24.R
# 输出：03_分析/logs/B6c_单折不稳定.log（由调用方重定向）
#       02_结果表/F11下游_ml12_单折波动.csv

suppressPackageStartupMessages({ library(pROC) })

BASE <- "/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/01.材料/01.博士材料/10.毕业/01.小论文/04.立项/滑膜基质刚度生成与刚度感知 0912"
P    <- file.path(BASE, "04_分析数据/processed")
OUT  <- file.path(BASE, "08_修订_2026-09-17/03_分析/outputs")
RES  <- file.path(BASE, "08_修订_2026-09-17/02_结果表")

PANEL <- c("LOXL1", "ITGA5", "SDC4", "SERPINE1", "SAV1", "VGLL4")
NSPLIT <- 3

g89 <- readRDS(file.path(OUT, "GSE89408_基因水平_双版本.rds"))
y <- factor(ifelse(g89$grp == "OA", "OA", "Normal"), levels = c("Normal", "OA"))
mats <- list(corrected = readRDS(file.path(OUT, "GSE89408_校正矩阵_cpm尺度.rds"))$mat,
             old_colSums = readRDS(file.path(P, "主分析_genelevel.rds"))$GSE89408$mat)

## 单次 5 折；norm = "whole" 全队列标准化 / "train" 仅训练折标准化（12 号脚本口径）
cv_once <- function(X, y, split, norm) {
  n <- nrow(X)
  set.seed(1000 + split)                          # 固定折种子，独立于任何其他脚本的 RNG 流
  folds <- sample(rep(1:5, length.out = n))
  p <- numeric(n); sep <- FALSE
  for (k in 1:5) {
    tr <- folds != k; te <- folds == k
    if (norm == "whole") {
      mu <- colMeans(X); sg <- apply(X, 2, sd)
    } else {
      mu <- colMeans(X[tr, , drop = FALSE]); sg <- apply(X[tr, , drop = FALSE], 2, sd)
    }
    sg[sg == 0] <- 1
    a <- scale(X[tr, , drop = FALSE], center = mu, scale = sg)
    b <- scale(X[te, , drop = FALSE], center = mu, scale = sg)
    colnames(a) <- colnames(b) <- paste0("p", seq_len(ncol(X)))
    d <- as.data.frame(a); d$y <- y[tr]
    m <- suppressWarnings(glm(y ~ ., data = d, family = "binomial"))
    sep <- sep || any(abs(coef(m)[-1]) > 20)
    p[te] <- predict(m, newdata = as.data.frame(b), type = "response")
  }
  c(auc = as.numeric(auc(roc(y, p, quiet = TRUE))), sep = as.numeric(sep))
}

rows <- list()
for (mn in names(mats)) {
  gn <- intersect(PANEL, rownames(mats[[mn]]))
  X <- t(mats[[mn]][gn, , drop = FALSE])
  for (nm in c("whole", "train")) {
    v <- sapply(1:NSPLIT, function(s) cv_once(X, y, s, nm))
    rows[[paste(mn, nm)]] <- data.frame(matrix = mn, norm = nm,
      split1 = round(v["auc", 1], 3), split2 = round(v["auc", 2], 3), split3 = round(v["auc", 3], 3),
      min = round(min(v["auc", ]), 3), max = round(max(v["auc", ]), 3),
      separation = paste(which(v["sep", ] == 1), collapse = "/"), n = nrow(X), row.names = NULL)
  }
}
res <- do.call(rbind, rows)
cat("GSE89408 面板 6 基因联合模型，单次 5 折 CV（固定折种子 1001/1002/1003）\n")
print(res, row.names = FALSE)
cat("\n对照：S2 原文记 0.515 / 0.503 / 0.468（12 号脚本注释另记 0.515/0.503/0.578）\n")

## ------------------------------------------------------------
## 附带发现须验证：两种标准化给出逐位相同的 AUC，说明折内重拟合时
## 标准化口径不影响结果（含截距的 logistic 对预测变量的仿射变换不变）。
## 这里把「预测概率逐元素相同」直接测出来，而不是只凭 AUC 相同推断。
## 注意：该不变性只适用于**折内重拟合**；跨队列迁移时系数固定，缩放口径会影响预测。
## ------------------------------------------------------------
gn <- intersect(PANEL, rownames(mats$corrected))
X <- t(mats$corrected[gn, , drop = FALSE])
set.seed(1001); folds <- sample(rep(1:5, length.out = nrow(X)))
pv <- function(norm) {
  p <- numeric(nrow(X))
  for (k in 1:5) {
    tr <- folds != k; te <- folds == k
    if (norm == "whole") { mu <- colMeans(X); sg <- apply(X, 2, sd) } else {
      mu <- colMeans(X[tr, , drop = FALSE]); sg <- apply(X[tr, , drop = FALSE], 2, sd) }
    sg[sg == 0] <- 1
    a <- scale(X[tr, , drop = FALSE], center = mu, scale = sg)
    b <- scale(X[te, , drop = FALSE], center = mu, scale = sg)
    colnames(a) <- colnames(b) <- paste0("p", seq_len(ncol(X)))
    d <- as.data.frame(a); d$y <- y[tr]
    m <- suppressWarnings(glm(y ~ ., data = d, family = "binomial"))
    p[te] <- predict(m, newdata = as.data.frame(b), type = "response")
  }
  p
}
pw <- pv("whole"); pt <- pv("train")
cat(sprintf("\n不变性检验（split 1，折内重拟合）：两种标准化的预测概率最大绝对差 = %.3e\n",
            max(abs(pw - pt))))
cat("（若为 ~1e-15 量级，则折内 CV AUC 与标准化口径无关；S2c 括注所述『全队列标准化会降低 AUC』不成立。）\n")

## 跨队列情形与之相反：系数固定，缩放口径会影响预测——用 GSE206848 作测试队列出示差异
V <- readRDS(file.path(P, "验证层_genelevel.rds"))
ALLt <- c(list(GSE206848 = list(mat = readRDS(file.path(P, "主分析_genelevel.rds"))$GSE206848$mat,
                                grp = readRDS(file.path(P, "主分析_genelevel.rds"))$GSE206848$grp)), V)
if (all(c("GSE206848") %in% names(ALLt)) && all(gn %in% rownames(ALLt$GSE206848$mat))) {
  mu <- colMeans(X); sg <- apply(X, 2, sd); sg[sg == 0] <- 1
  a <- scale(X, center = mu, scale = sg); colnames(a) <- paste0("p", seq_along(gn))
  dd <- as.data.frame(a); dd$y <- y
  mm <- suppressWarnings(glm(y ~ ., data = dd, family = "binomial"))
  Xt <- t(ALLt$GSE206848$mat[gn, , drop = FALSE])
  yt <- factor(ifelse(ALLt$GSE206848$grp == "OA", "OA", "Normal"), levels = c("Normal", "OA"))
  sc <- function(cen, scl) { z <- scale(Xt, center = cen, scale = scl)
    colnames(z) <- paste0("p", seq_along(gn))
    predict(mm, newdata = as.data.frame(z), type = "response") }
  p_tr <- sc(mu, sg)                                      # 训练队列参数（S2c 实际做法）
  p_te <- sc(colMeans(Xt), apply(Xt, 2, sd))              # 测试队列自身参数
  a_tr <- as.numeric(auc(roc(yt, p_tr, quiet = TRUE)))
  a_te <- as.numeric(auc(roc(yt, p_te, quiet = TRUE)))
  cat(sprintf("跨队列对照（GSE89408 训练 → GSE206848 测试，n=14）：\n"))
  cat(sprintf("  训练队列参数：AUC=%.4f；测试队列自身参数：AUC=%.4f（该队列只有 7/7，AUC 分辨率 1/49=0.0204，\n", a_tr, a_te))
  cat(sprintf("    两者落到同一格不构成『不变』的证据）\n"))
  cat(sprintf("  但预测概率最大绝对差 = %.3e（远大于折内情形的 1e-16）→ 迁移时缩放口径确实进入预测。\n",
              max(abs(p_tr - p_te))))
}
write.csv(res, file.path(RES, "F11下游_ml12_单折波动.csv"), row.names = FALSE, fileEncoding = "UTF-8")
