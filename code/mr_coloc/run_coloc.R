# coloc 共定位：cis-eQTL vs FinnGen KOA，是否共享因果位点
#
# 用法：Rscript run_coloc.R [GTEx|Alasoo|Alasoo_orig]     默认 GTEx（主分析）
# 输入：coloc_input_{gene}_{in_tag}.json
# 输出：coloc_results_{run}.csv
#
# 与原版的差异：
#  ① 输入不再写死 /tmp/merge_data（重启即丢），改为读本目录下 extract_all.py 的产出；
#  ② 样本量不再写死。原版用 N=500（eQTL）与 N=400000, s=0.08（KOA）三个常数，都与
#     实际数据无关：eQTL 分别是 483（GTEx）与 84（Alasoo），KOA 是 376,471 例对照合计、
#     病例占比 0.163 而非 0.08。s 差一倍会直接改变 ABF，进而改变 PP.H4。
#  ③ 结果落盘成 CSV —— 原版只 print，事后无法核验；原 coloc_results.csv 里
#     CYR61/SDC4 两行写的是 "~1" "~0" 这种手填近似值，不是脚本产出。
#  ④ 新增 Alasoo_orig 一档：用原版的错误参数跑同一份输入，用来把「换暴露数据」和
#     「修参数」两个改动的影响分开量化，该档只作诊断，不进正文。
library(coloc)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
run  <- if (length(args) >= 1) args[1] else "GTEx"
dir  <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))

# FinnGen R12 M13_ARTHROSIS_KNEE：官方 manifest 载明 61,356 例 / 315,115 对照
N_KOA <- 61356 + 315115
S_KOA <- 61356 / N_KOA

RUNS <- list(
  Synovium    = list(in_tag = "Synovium", N_eqtl = 72),   # 主分析：人 OA 滑膜
  GTEx        = list(in_tag = "GTEx",     N_eqtl = 483),
  Alasoo      = list(in_tag = "Alasoo",   N_eqtl = 84),
  Alasoo_orig = list(in_tag = "Alasoo",   N_eqtl = 500)   # 诊断用：原版常数
)
cfg <- RUNS[[run]]
if (is.null(cfg)) stop("未知 run：", run)
if (run == "Alasoo_orig") { N_koa <- 400000; s_koa <- 0.08 } else { N_koa <- N_KOA; s_koa <- S_KOA }
cat(sprintf("run=%s  输入 tag=%s  N_eqtl=%d  N_koa=%d  s_koa=%.4f\n",
            run, cfg$in_tag, cfg$N_eqtl, N_koa, s_koa))

res_all <- list()
# 基因清单的唯一来源是 gene_coords.py 的 COLOC_GENES。此处不再各抄一份：本目录
# gene_coords.py 的文件头记的就是「两处各抄一份、只改了一处」导致 SDC4 窗口取错的事故。
COLOC_GENES <- c("VGLL4", "SDC4", "SERPINE1")     # 2026-09-17 补入 SERPINE1
gc_line <- grep("^COLOC_GENES", readLines(file.path(dir, "gene_coords.py")), value = TRUE)
stopifnot(length(gc_line) == 1L,
          identical(COLOC_GENES, gsub('"', "", regmatches(gc_line, gregexpr('"[A-Za-z0-9]+"', gc_line))[[1]])))
# 缺输入必须报错，不能 next 跳过：静默跳过会产出一份少了一个基因、看着完全正常的结果表。
have <- file.exists(file.path(dir, paste0("coloc_input_", COLOC_GENES, "_", cfg$in_tag, ".json")))
if (!all(have)) stop("缺 coloc 输入：", paste(COLOC_GENES[!have], collapse = ", "),
                     " —— 先跑 extract_all.py ", cfg$in_tag)

for (sym in COLOC_GENES) {
  f <- file.path(dir, paste0("coloc_input_", sym, "_", cfg$in_tag, ".json"))
  d <- fromJSON(f)
  if (nrow(d) < 3) { cat(sym, ": 交集 SNP 太少(<3)，跳过\n"); next }
  d1 <- list(snp = d$rsid, beta = d$beta_eqtl, varbeta = d$se_eqtl^2,
             MAF = d$maf_eqtl, N = cfg$N_eqtl, type = "quant")
  d2 <- list(snp = d$rsid, beta = d$beta_koa, varbeta = d$se_koa^2,
             MAF = d$maf_koa, N = N_koa, type = "cc", s = s_koa)
  s <- coloc.abf(d1, d2)$summary
  cat(sprintf("\n=== %s coloc (%d SNP) ===\n", sym, s["nsnps"]))
  for (h in c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf"))
    cat(sprintf("  %-10s %.3f\n", sub("\\.abf$", "", h), s[h]))
  res_all[[sym]] <- data.frame(gene = sym, n_snp = s["nsnps"],
                               PP.H0 = s["PP.H0.abf"], PP.H1 = s["PP.H1.abf"],
                               PP.H2 = s["PP.H2.abf"], PP.H3 = s["PP.H3.abf"],
                               PP.H4 = s["PP.H4.abf"], row.names = NULL)
}

if (length(res_all)) {
  out <- file.path(dir, paste0("coloc_results_", run, ".csv"))
  write.csv(do.call(rbind, res_all), out, row.names = FALSE)
  cat("\n已保存", basename(out), "\n")
}
