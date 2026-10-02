# 从落盘的 CSV 生成 Table 2 / Table 3 的 markdown 行，避免手抄数字。
#
# 用法：python3 15_gen_tables.py [--cn] [--raw]
#   默认读 LD 剪枝后的结果（*_pruned.csv）——那是统计上有效的主分析；
#   --raw 读未剪枝的旧结果，仅供对照，不用于正文。
import csv, sys

SUF = "" if "--raw" in sys.argv else "_pruned"

CN = "--cn" in sys.argv
EXPO = [("Synovium", "滑膜 (n=72)" if CN else "Synovium (n=72)"),
        ("GTEx", "GTEx 成纤维 (n=483)" if CN else "GTEx fibroblast (n=483)"),
        ("Alasoo", "巨噬 (n=84)" if CN else "Macrophage (n=84)")]
GENES = ["VGLL4", "SDC4", "LOXL1", "SERPINE1", "SAV1", "ITGA5"]

def fmt(b, se):  return f"{b:+.3f} ({se:.4f})"

def fnum(p):
    p = float(p)
    if p < 1e-3:
        e = 0; m = p
        while m < 1: m *= 10; e -= 1
        return f"{m:.1f}×10⁻{abs(e)}"
    return f"{p:.3f}"

def clean(v):    return "—" if v in ("", "NA", None) else v

def main_row(tag, g):
    """按 method 显式取主估计，不按行号。"""
    rs = [r for r in csv.DictReader(open(f"multi_iv_MR_results_{tag}{SUF}.csv")) if r["gene"] == g]
    for m in ("IVW", "Wald"):
        hit = [r for r in rs if r["method"] == m]
        if hit: return hit[0], m
    raise SystemExit(f"✗ {tag} {g} 无 IVW/Wald 行：{[r['method'] for r in rs]}")

coloc = {t: {r["gene"]: float(r["PP.H4"]) for r in csv.DictReader(open(f"coloc_results_{t}.csv"))}
         for t in ("Synovium", "GTEx", "Alasoo")}

print("| 基因 | 暴露 | 工具变量 | 方法 | β (SE) | P | 多效性截距 P | 基因水平 F | PP.H4 |" if CN else
      "| Gene | Exposure | n IV | Method | β (SE) | P | Intercept P | F | PP.H4 |")
print("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
for g in GENES:
    for tag, label in EXPO:
        r, _ = main_row(tag, g)
        h4 = coloc[tag].get(g)
        F = clean(r["Fstat"])
        print(f"| {g} | {label} | {r['n_iv']} | {r['method']} "
              f"| {fmt(float(r['beta']), float(r['se']))} | {fnum(r['p'])} "
              f"| {fnum(r['pleio_p']) if clean(r['pleio_p'])!='—' else '—'} "
              f"| {f'{float(F):.1f}' if F!='—' else '—'} "
              f"| {f'{h4:.3f}' if h4 is not None else '—'} |")

print("\n<<<TABLE3>>>")
print("| 暴露 | 基因 | SNP | 效应/其他等位基因 | EAF | eQTL β (SE) | eQTL P | KOA β (SE) | 基因水平 F |" if CN else
      "| Exposure | Gene | SNP | Effect/other allele | EAF | eQTL β (SE) | eQTL P | KOA β (SE) | Gene-level F |")
print("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
for tag, label in EXPO:
    for g in GENES:
        _, _ = main_row(tag, g)
        F = clean(main_row(tag, g)[0]["Fstat"])
        for d in csv.DictReader(open(f"multi_iv_exposure_outcome_{tag}{SUF}.csv")):
            if d["gene"] != g: continue
            print(f"| {label} | {g} | {d['SNP']} | {d['effect_allele']}/{d['other_allele']} "
                  f"| {float(d['eaf_out']):.3f} | {fmt(float(d['beta_exp']), float(d['se_exp']))} "
                  f"| {fnum(d['p_exp'])} | {fmt(float(d['beta_out']), float(d['se_out']))} "
                  f"| {f'{float(F):.1f}' if F!='—' else '—'} |")
