# -*- coding: utf-8 -*-
"""一次抓取，同时产出多工具变量 MR 的输入与 coloc 的输入。

取代原先分开的 extract_multi_iv.py 与 extract_coloc.py。合并有两个理由：
① 后者那份基因坐标是 GRCh37 的旧值、且未与前者同步，SDC4 的共定位窗口因此取错
   （详见 gene_coords.py）；
② EBI 对该 sumstats 限速约 11 KB/s，一次 ±500 kb 抓取要一分钟以上，两次抓取里的
   重叠基因（VGLL4、SDC4）不值得各抓一遍。

用法：python3 extract_all.py [GTEx|Alasoo]     默认 GTEx（主分析）

产出（均落在本目录）：
  cache/eqtl_{tag}_{gene}.tsv        原始 cis-eQTL 行，落盘；重跑直接读缓存，不再联网
  cache/finngen_{tag}.tsv            缓存基因窗内命中的 FinnGen 行
  multi_iv_exposure_outcome_{tag}.csv
  coloc_input_{gene}_{tag}.json
"""
import csv, gzip, hashlib, json, os, sys
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gene_coords import EXPOSURES, GENES, MR_GENES, COLOC_GENES, FETCH, FINNGEN

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache")
os.makedirs(CACHE, exist_ok=True)

tag = sys.argv[1] if len(sys.argv) > 1 else "GTEx"
exp = EXPOSURES[tag]
GENES_NEEDED = sorted(set(MR_GENES) | set(COLOC_GENES))
print(f"暴露数据：{exp['label']}  n={exp['n']}")

# ── 1. 抓 cis-eQTL（逐基因缓存；同一 rsid 只留首次出现，与原脚本一致） ──────
tb = pysam.TabixFile(exp["url"])
eqtl = {}
for sym in GENES_NEEDED:
    path = os.path.join(CACHE, f"eqtl_{tag}_{sym}.tsv")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        rows = [l.rstrip("\n").split("\t") for l in open(path)]
        print(f"{sym}: 读缓存 {len(rows)} 个 cis-eQTL")
    else:
        ensg, chrom, pos = GENES[sym]
        lo, hi = pos - FETCH, pos + FETCH
        seen, rows = set(), []
        for r in tb.fetch(chrom, lo, hi):
            f = r.split("\t")
            if f[0] == ensg and f[18] not in seen:
                seen.add(f[18])
                rows.append([f[2], f[3], f[4], f[7], f[8], f[9], f[10], f[18]])
        with open(path, "w") as fh:
            fh.write("".join("\t".join(r) + "\n" for r in rows))
        print(f"{sym}: 已抓取 {len(rows)} 个 cis-eQTL → 缓存")
    if not rows:
        sys.exit(f"✗ {sym} 在 ±{FETCH//1000} kb 窗内没有 cis-eQTL，坐标可能有问题")
    # 列: pos ref alt maf p beta se rsid
    eqtl[sym] = [dict(pos=int(r[0]), ref=r[1], alt=r[2], maf=float(r[3]), p=float(r[4]),
                      beta=float(r[5]), se=float(r[6]), rsid=r[7]) for r in rows]
    lo, hi = GENES[sym][2] - FETCH, GENES[sym][2] + FETCH
    pmin, pmax = min(d["pos"] for d in eqtl[sym]), max(d["pos"] for d in eqtl[sym])
    print(f"     位置 {pmin:,}–{pmax:,}，距窗口两端 {pmin-lo:,} / {hi-pmax:,} bp"
          + ("   ⚠️ 撞边，窗口是限制因素" if min(pmin - lo, hi - pmax) < 50_000 else ""))

# ── 2. 选工具变量：按 P 排序取 top3，两两相距 >100 kb（近似 LD 独立） ────────
gene_iv = {}
for sym in MR_GENES:
    selected = []
    for d in sorted(eqtl[sym], key=lambda x: x["p"]):
        if all(abs(d["pos"] - s["pos"]) > 100_000 for s in selected):
            selected.append(d)
        if len(selected) >= 3:
            break
    gene_iv[sym] = selected
    print(f"{sym}: 选 {len(selected)} 个工具变量, 最小 P={selected[0]['p']:.2e}")

# ── 3. FinnGen：全表扫一遍，取出上一步所有候选 rsid 的行 ────────────────────
# 缓存以「候选 rsid 集合」的哈希为键，而不是以暴露数据 tag 为键：三种暴露查的是同一份
# 结局 GWAS，但候选 rsid 各不相同。按 tag 缓存会让换了暴露仍读到旧的（可能缺行的）结果；
# 按集合哈希缓存则集合一变必然落空，不会静默复用。
need = {d["rsid"] for sym in GENES_NEEDED for d in eqtl[sym]}
fg_key = hashlib.sha1("|".join(sorted(need)).encode()).hexdigest()[:12]
fg_path = os.path.join(CACHE, f"finngen_{fg_key}.tsv")
if os.path.exists(fg_path):
    fg = {}
    for l in open(fg_path):
        p = l.rstrip("\n").split("\t")
        fg[p[4]] = dict(chrom=p[0], pos=int(p[1]), ref=p[2], alt=p[3],
                        beta=float(p[8]), se=float(p[9]), maf=float(p[10]))
    print(f"FinnGen: 读缓存 {len(fg)} 行")
else:
    print(f"FinnGen: 扫描全表（约 819 MB，需数分钟）…")
    fg, dup = {}, 0
    with gzip.open(os.path.join(HERE, FINNGEN), "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 11:
                continue
            for rs in p[4].split(","):
                if rs in need:
                    if rs in fg:
                        dup += 1
                    else:
                        fg[rs] = dict(chrom=p[0], pos=int(p[1]), ref=p[2], alt=p[3],
                                      beta=float(p[8]), se=float(p[9]), maf=float(p[10]))
    if dup:
        print(f"  ⚠️ FinnGen 中有 {dup} 个重复 rsid，已保留首次出现")
    # 落盘时保持与 sumstats 相同的列序：rsid 在 4、beta/se/maf 在 8/9/10
    with open(fg_path, "w") as fh:
        for rs, v in fg.items():
            fh.write("\t".join([v["chrom"], str(v["pos"]), v["ref"], v["alt"], rs,
                                ".", ".", ".", str(v["beta"]), str(v["se"]),
                                str(v["maf"])]) + "\n")
    print(f"FinnGen: 命中 {len(fg)} / {len(need)} 个候选 rsid → 缓存 {fg_key}")


def harmonise(e, g):
    """对齐效应等位基因；不匹配返回 None。"""
    if e["ref"] == g["ref"] and e["alt"] == g["alt"]:
        return g["beta"], g["maf"]
    if e["ref"] == g["alt"] and e["alt"] == g["ref"]:
        return -g["beta"], 1 - g["maf"]
    return None


# ── 4. 输出 MR 输入 ───────────────────────────────────────────────────────
rows = []
for sym in MR_GENES:
    for d in gene_iv[sym]:
        g = fg.get(d["rsid"])
        if not g:
            print(f"  {sym} {d['rsid']}: FinnGen 无匹配, 跳过"); continue
        h = harmonise(d, g)
        if h is None:
            print(f"  {sym} {d['rsid']}: 等位基因不匹配, 跳过"); continue
        beta_out, eaf_out = h
        rows.append(dict(gene=sym, SNP=d["rsid"], pos=d["pos"],
                         beta_exp=d["beta"], se_exp=d["se"], p_exp=d["p"],
                         beta_out=beta_out, se_out=g["se"],
                         effect_allele=d["alt"], other_allele=d["ref"], eaf_out=eaf_out))
out = os.path.join(HERE, f"multi_iv_exposure_outcome_{tag}.csv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(f"\n→ {os.path.basename(out)}：{len(rows)} 行")
for r in rows:
    print(f"   {r['gene']:9s} {r['SNP']:14s} eQTL P={r['p_exp']:.2e} "
          f"beta_exp={r['beta_exp']:+.3f} beta_out={r['beta_out']:+.3f}")

# ── 5. 输出 coloc 输入 ────────────────────────────────────────────────────
print()
for sym in COLOC_GENES:
    common = []
    for d in sorted(eqtl[sym], key=lambda x: x["pos"]):
        g = fg.get(d["rsid"])
        if not g:
            continue
        h = harmonise(d, g)
        if h is None:
            continue
        beta_koa, _ = h
        common.append(dict(rsid=d["rsid"], pos=d["pos"],
                           beta_eqtl=d["beta"], se_eqtl=d["se"], maf_eqtl=d["maf"],
                           beta_koa=beta_koa, se_koa=g["se"], maf_koa=g["maf"]))
    p = os.path.join(HERE, f"coloc_input_{sym}_{tag}.json")
    json.dump(common, open(p, "w"))
    print(f"{sym}: eQTL {len(eqtl[sym])} ∩ FinnGen → {len(common)} 个 SNP → "
          f"{os.path.basename(p)}")
