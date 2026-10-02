# MR 敏感性检查：逐工具 Wald 对比 + 更严格剪枝（r²≤0.01）
#
# 目的：回应「2 变异 IVW 的残余相关只使 SE 轻微反保守」这一未经计算的断言。
# 做法一：把每个多工具 IVW 拆成其组成单变异的 Wald 比，看 IVW 是否由相互
#         一致的估计构成 —— 若两个单变异估计相差悬殊，IVW 的 P 值即不能独立
#         于该对工具的相关结构来解读。
# 做法二：用 ld_instruments.csv 里实测的逐对 r²，把阈值收紧到 0.01 重做贪心
#         剪枝，看哪些集合退化为单工具、退化后估计如何变化。
# 估计量：MendelianRandomization::mr_ivw，与主分析 (run_multi_iv_MR.R) 完全一致。
suppressPackageStartupMessages(library(MendelianRandomization))
DIR <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
TAGS <- c(GTEx = "GTEx", Alasoo = "Alasoo", Synovium = "Synovium")
ld <- read.csv(file.path(DIR, "ld_instruments.csv"), stringsAsFactors = FALSE)

wald <- function(bx, bxse, by, byse) {
  b <- by / bx
  se <- abs(byse / bx)          # 与主脚本同口径（一阶近似）
  c(beta = b, se = se, p = 2 * pnorm(-abs(b / se)))
}

rows <- list()
for (tag in names(TAGS)) {
  dat <- read.csv(file.path(DIR, paste0("multi_iv_exposure_outcome_", tag, "_pruned.csv")),
                  stringsAsFactors = FALSE)
  for (g in unique(dat$gene)) {
    d <- dat[dat$gene == g, ]
    n <- nrow(d)
    # --- 当前 IVW ---
    if (n >= 2) {
      mri <- mr_input(bx = d$beta_exp, bxse = d$se_exp, by = d$beta_out, byse = d$se_out)
      ivw <- mr_ivw(mri)
      ivw_b <- ivw@Estimate; ivw_se <- ivw@StdError; ivw_p <- ivw@Pvalue
    } else {
      w <- wald(d$beta_exp, d$se_exp, d$beta_out, d$se_out)
      ivw_b <- w["beta"]; ivw_se <- w["se"]; ivw_p <- w["p"]
    }
    # --- 组成单变异 Wald ---
    for (i in seq_len(n)) {
      w <- wald(d$beta_exp[i], d$se_exp[i], d$beta_out[i], d$se_out[i])
      rows[[length(rows) + 1]] <- data.frame(
        exposure = tag, gene = g, n_iv = n, snp = d$SNP[i],
        eqtl_p = d$p_exp[i], wald_beta = w["beta"], wald_se = w["se"], wald_p = w["p"],
        ivw_beta = ivw_b, ivw_se = ivw_se, ivw_p = ivw_p, role = "component")
    }
    # --- 更严格剪枝 r²<=0.01 ---
    if (n >= 2) {
      ord <- order(d$p_exp)                       # 与 14_prune_iv.py 同序：按 eQTL P 升序
      keep <- character(0)
      for (i in ord) {
        s <- d$SNP[i]; drop <- FALSE
        for (k in keep) {
          r <- ld$r2[ld$gene == g & ld$exposure == tag &
                     ((ld$snp1 == s & ld$snp2 == k) | (ld$snp1 == k & ld$snp2 == s))]
          if (length(r) && !is.na(r) && r > 0.01) { drop <- TRUE; break }
          if (length(r) && is.na(r)) { drop <- TRUE; break }   # LD 不可核验同样剔除
        }
        if (!drop) keep <- c(keep, s)
      }
      ds <- d[d$SNP %in% keep, ]
      if (nrow(ds) >= 2) {
        mri2 <- mr_input(bx = ds$beta_exp, bxse = ds$se_exp, by = ds$beta_out, byse = ds$se_out)
        i2 <- mr_ivw(mri2)
        r2b <- i2@Estimate; r2se <- i2@StdError; r2p <- i2@Pvalue
      } else {
        w <- wald(ds$beta_exp, ds$se_exp, ds$beta_out, ds$se_out)
        r2b <- w["beta"]; r2se <- w["se"]; r2p <- w["p"]
      }
      rows[[length(rows) + 1]] <- data.frame(
        exposure = tag, gene = g, n_iv = nrow(ds), snp = paste(keep, collapse = ";"),
        eqtl_p = NA, wald_beta = r2b, wald_se = r2se, wald_p = r2p,
        ivw_beta = ivw_b, ivw_se = ivw_se, ivw_p = ivw_p,
        role = paste0("strict_r2_0.01 (n_iv ", n, "->", nrow(ds), ")"))
    }
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path(DIR, "mr_sensitivity_component.csv"), row.names = FALSE)

# --- 逐对实测 r² 汇总（多工具集合用）---
pairs <- data.frame()
for (tag in names(TAGS)) {
  dat <- read.csv(file.path(DIR, paste0("multi_iv_exposure_outcome_", tag, "_pruned.csv")),
                  stringsAsFactors = FALSE)
  for (g in unique(dat$gene)) {
    d <- dat[dat$gene == g, ]
    if (nrow(d) < 2) next
    for (i in 1:(nrow(d) - 1)) for (j in (i + 1):nrow(d)) {
      a <- d$SNP[i]; b <- d$SNP[j]
      r <- ld$r2[ld$gene == g & ld$exposure == tag &
                 ((ld$snp1 == a & ld$snp2 == b) | (ld$snp1 == b & ld$snp2 == a))]
      pairs <- rbind(pairs, data.frame(exposure = tag, gene = g, snp1 = a, snp2 = b,
                                       r2 = if (length(r)) r[1] else NA_real_))
    }
  }
}
write.csv(pairs, file.path(DIR, "mr_sensitivity_residual_r2.csv"), row.names = FALSE)
cat("已保存 mr_sensitivity_component.csv 与 mr_sensitivity_residual_r2.csv\n")
