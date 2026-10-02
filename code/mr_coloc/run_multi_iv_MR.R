# 多工具变量 MR：≥3 SNP 做 IVW+Egger+Median；2 SNP 做 IVW；1 SNP 做 Wald
#
# 用法：Rscript run_multi_iv_MR.R [GTEx|Alasoo|Synovium] [后缀]
#      后缀缺省为空 → 输入 multi_iv_exposure_outcome_{tag}.csv
#                     输出 multi_iv_MR_results_{tag}.csv
#      后缀 "_pruned" → 读 LD 剪枝后的独立工具集（14_prune_iv.py 产出）
# 输入：multi_iv_exposure_outcome_{tag}{后缀}.csv
# 输出：multi_iv_MR_results_{tag}{后缀}.csv
#
# 与原版的差异：① 按 tag 取输入/输出，不再写死路径（原版指向已改名的旧目录，跑不起来）；
# ② 一并落盘 Cochran's Q —— 原版只把 Q 打进 stdout，正文引用的 P=0.932 因此没有
#    机器可读的出处，事后无法核验。
suppressPackageStartupMessages(library(MendelianRandomization))

args <- commandArgs(trailingOnly = TRUE)
tag  <- if (length(args) >= 1) args[1] else "GTEx"
suf  <- if (length(args) >= 2) args[2] else ""
dir  <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))

dat <- read.csv(file.path(dir, paste0("multi_iv_exposure_outcome_", tag, suf, ".csv")),
                stringsAsFactors = FALSE)
cat("暴露数据 tag =", tag, "后缀 =", if (nzchar(suf)) suf else "(无)", "，共", nrow(dat), "行\n")

res_all <- list()
for (g in c("VGLL4","SDC4","LOXL1","SERPINE1","SAV1","ITGA5")) {
  d <- dat[dat$gene == g, ]
  cat("\n=====", g, "(", nrow(d), "个工具变量 ) =====\n")
  mri <- mr_input(bx = d$beta_exp, bxse = d$se_exp, by = d$beta_out, byse = d$se_out,
                  exposure = "gene expression", outcome = "KOA", snps = d$SNP)
  if (nrow(d) >= 3) {
    ivw <- mr_ivw(mri); egger <- mr_egger(mri); med <- mr_median(mri)
    het <- tryCatch(as.numeric(ivw@Heter.Stat), error = function(e) c(NA_real_, NA_real_))
    cat(sprintf("IVW:    beta=%.4f se=%.4f p=%.4g (F=%.1f)\n", ivw@Estimate, ivw@StdError, ivw@Pvalue, ivw@Fstat))
    cat(sprintf("Egger:  beta=%.4f se=%.4f p=%.4g (多效性截距 p=%.3f)\n", egger@Estimate, egger@StdError.Est, egger@Pvalue.Est, egger@Pleio.pval))
    cat(sprintf("Median: beta=%.4f se=%.4f p=%.4g\n", med@Estimate, med@StdError, med@Pvalue))
    cat(sprintf("Cochran Q: %.3f, p=%.3f\n", het[1], het[2]))
    add <- function(method, b, se, p) data.frame(
      gene=g, n_iv=nrow(d), method=method, beta=b, se=se, p=p,
      pleio_p=egger@Pleio.pval, Fstat=ivw@Fstat, Q=het[1], Q_p=het[2])
    res_all[[g]] <- rbind(add("IVW",    ivw@Estimate,   ivw@StdError,    ivw@Pvalue),
                          add("Egger",  egger@Estimate, egger@StdError.Est, egger@Pvalue.Est),
                          add("Median", med@Estimate,   med@StdError,    med@Pvalue))
  } else if (nrow(d) == 2) {
    ivw <- mr_ivw(mri)
    het <- tryCatch(as.numeric(ivw@Heter.Stat), error = function(e) c(NA_real_, NA_real_))
    cat(sprintf("IVW:    beta=%.4f se=%.4f p=%.4g (F=%.1f)\n", ivw@Estimate, ivw@StdError, ivw@Pvalue, ivw@Fstat))
    res_all[[g]] <- data.frame(gene=g, n_iv=2, method="IVW", beta=ivw@Estimate,
                               se=ivw@StdError, p=ivw@Pvalue, pleio_p=NA,
                               Fstat=ivw@Fstat, Q=het[1], Q_p=het[2])
  } else {
    b <- d$beta_out[1] / d$beta_exp[1]; se <- abs(d$se_out[1] / d$beta_exp[1]); p <- 2*pnorm(-abs(b/se))
    cat(sprintf("Wald:   beta=%.4f se=%.4f p=%.4g\n", b, se, p))
    res_all[[g]] <- data.frame(gene=g, n_iv=1, method="Wald", beta=b, se=se, p=p,
                               pleio_p=NA, Fstat=NA, Q=NA, Q_p=NA)
  }
}

out <- do.call(rbind, res_all)
f <- file.path(dir, paste0("multi_iv_MR_results_", tag, suf, ".csv"))
write.csv(out, f, row.names = FALSE)
cat("\n已保存", basename(f), "\n")
