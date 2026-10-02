# -*- coding: utf-8 -*-
"""按实测 LD 剪枝工具变量，产出独立工具集。

为什么必须做：§2.5 原先用「物理距离 ≥100 kb」代表近似独立。实测（13_ld_check.py）
发现 6 对工具 r²>0.1，其中 Alasoo VGLL4 的 rs347609–rs7624357 达 r²=0.69 —— 那正是
文中 P=1.8×10⁻¹⁰ 所依赖的 3 个工具中的 2 个。相关工具会让 IVW 的标准误被低估、
P 值被高估。

策略（保守）：
  · 工具按 eQTL P 升序（最重要的在前）贪心保留；
  · 与任一已保留工具 r²>0.1 者丢弃；
  · LD 无法核实的（不在 1000G 参考里）从多工具集中丢弃——宁可少一个工具，
    不可在独立性未知的情况下当作独立工具求和。
剪枝后 1 个工具 → Wald；≥2 个 → IVW（此时已可视为近似独立）。

不改动原文件，结果写入 *_pruned.csv，并落盘剪枝日志。
"""
import csv
import os
import glob
import re

MR = os.path.dirname(os.path.abspath(__file__))
R2_MAX = 0.10

# 读实测 LD：(暴露, 基因, 无序SNP对) → r²
ld = {}
for r in csv.DictReader(open(f"{MR}/ld_instruments.csv", encoding="utf-8")):
    if r["r2"] in ("", "NA"):
        continue                      # 缺失 = 无法核实，调用处按缺失处理
    k = (r["exposure"], r["gene"], frozenset((r["snp1"], r["snp2"])))
    ld[k] = float(r["r2"])


def r2_of(exp, gene, a, b):
    return ld.get((exp, gene, frozenset((a, b))))   # None 表示无法核实


log, all_drops = [], []
for path in sorted(glob.glob(f"{MR}/multi_iv_exposure_outcome_*.csv")):
    if "_pruned" in path:
        continue
    tag = re.search(r"_outcome_(\w+)\.csv", path).group(1)
    rows = list(csv.DictReader(open(path, encoding="utf-8")))
    by = {}
    for r in rows:
        by.setdefault(r["gene"], []).append(r)

    kept_all, dropped_all = [], []
    for gene, rs in sorted(by.items()):
        # eQTL P 升序：最重要的工具优先保留
        rs = sorted(rs, key=lambda x: float(x["p_exp"]))
        kept = []
        for cand in rs:
            why = None
            for k in kept:
                m = r2_of(tag, gene, cand["SNP"], k["SNP"])
                if m is None:
                    why = f"与 {k['SNP']} 的 LD 无法核实（不在 1000G 参考中）"
                    break
                if m > R2_MAX:
                    why = f"与 {k['SNP']} 相关 r²={m:.3f}"
                    break
            if why:
                dropped_all.append((gene, cand["SNP"], why))
                all_drops.append((tag, gene, cand["SNP"], why))
            else:
                kept.append(cand)
        kept_all.extend(kept)
        log.append("%-9s %-9s 工具 %d → %d  保留 %s"
                   % (tag, gene, len(rs), len(kept), ",".join(x["SNP"] for x in kept)))

    fn = ["gene", "SNP", "pos", "beta_exp", "se_exp", "p_exp",
          "beta_out", "se_out", "effect_allele", "other_allele", "eaf_out"]
    out = f"{MR}/multi_iv_exposure_outcome_{tag}_pruned.csv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fn)
        w.writeheader()
        for r in kept_all:
            w.writerow({k: r[k] for k in fn})
    print(f"✓ {tag}: 保留 {len(kept_all)} 个工具（原 {len(rows)}），丢弃 {len(dropped_all)}")

print()
print("被丢弃的工具及理由：")
for tag, gene, snp, why in all_drops:
    print(f"   {tag:<9} {gene:<9} {snp:<14} {why}")

open(f"{MR}/prune_log.txt", "w", encoding="utf-8").write(
    "工具变量 LD 剪枝日志（阈值 r²>" + str(R2_MAX) + "）\n\n"
    + "\n".join(log) + "\n\n被丢弃的工具及理由：\n"
    + "\n".join(f"   {t:<9} {g:<9} {s:<14} {w}" for t, g, s, w in all_drops) + "\n")
print(f"\n✓ 已落盘 prune_log.txt 与 3 个 *_pruned.csv")
