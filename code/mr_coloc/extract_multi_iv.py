# 提取多工具变量 MR 输入：GTEx cis-eQTL(每基因 top3 近似独立 SNP) + FinnGen KOA
import pysam, gzip

genes = {
    "VGLL4":    ("ENSG00000144560", "3",  11556068),
    "SDC4":     ("ENSG00000124145", "20", 45325288),   # 修正原稿 43952000
    "LOXL1":    ("ENSG00000129038", "15", 73926461),
    "SERPINE1": ("ENSG00000106366", "7",  101127103),
    "SAV1":     ("ENSG00000151748", "14", 50633579),
    "ITGA5":    ("ENSG00000161638", "12", 54395260),
}

# ---- 1. GTEx cis-eQTL 提取 ----
url = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000001/QTD000001/QTD000001.all.tsv.gz"
tb = pysam.TabixFile(url)

gene_iv = {}
for sym, (ensg, chrom, pos) in genes.items():
    lo, hi = pos - 500000, pos + 500000
    eqtls = {}
    for r in tb.fetch(chrom, lo, hi):
        f = r.split("\t")
        if f[0] == ensg and f[18] not in eqtls:
            eqtls[f[18]] = {"pos": int(f[2]), "ref": f[3], "alt": f[4],
                            "maf": float(f[7]), "p": float(f[8]),
                            "beta": float(f[9]), "se": float(f[10])}
    # 按 p 排序, 取 top3, 相邻距离 > 100kb (近似 LD 独立)
    selected = []
    for rsid, info in sorted(eqtls.items(), key=lambda x: x[1]["p"]):
        if all(abs(info["pos"] - s["pos"]) > 100000 for s in selected):
            selected.append({"rsid": rsid, **info})
        if len(selected) >= 3:
            break
    gene_iv[sym] = selected
    print(f"{sym}: cis-eQTL {len(eqtls)} 个, 选 {len(selected)} 个工具变量")

# ---- 2. FinnGen GWAS 提取(遍历一次) ----
needed = {}
for sym, selected in gene_iv.items():
    for s in selected:
        needed[s["rsid"]] = (sym, s)

finngen = {}
with gzip.open("/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/01.材料/01.博士材料/10.毕业/01.小论文/04.立项/KOA 滑膜力学转导基因 0912/04_分析数据/MR_coloc/mr_data/finngen_KOA.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 11:
            continue
        rsids = p[4].split(",")
        for rs in rsids:
            if rs in needed:
                finngen[rs] = {"chrom": p[0], "pos": int(p[1]), "ref": p[2], "alt": p[3],
                               "beta": float(p[8]), "se": float(p[9]), "maf": float(p[10])}

# ---- 3. 等位基因对齐 + 输出 ----
rows = []
for sym, selected in gene_iv.items():
    for s in selected:
        rs = s["rsid"]
        if rs not in finngen:
            print(f"  {sym} {rs}: FinnGen 无匹配, 跳过")
            continue
        e, g = s, finngen[rs]
        # 对齐: eQTL ref/alt vs GWAS ref/alt
        if e["ref"] == g["ref"] and e["alt"] == g["alt"]:
            beta_out, eaf_out = g["beta"], g["maf"]
        elif e["ref"] == g["alt"] and e["alt"] == g["ref"]:
            beta_out, eaf_out = -g["beta"], 1 - g["maf"]
        else:
            print(f"  {sym} {rs}: 等位基因不匹配 ({e['ref']}/{e['alt']} vs {g['ref']}/{g['alt']}), 跳过")
            continue
        rows.append({
            "gene": sym, "SNP": rs, "pos": e["pos"],
            "beta_exp": e["beta"], "se_exp": e["se"], "p_exp": e["p"],
            "beta_out": beta_out, "se_out": g["se"],
            "effect_allele": e["alt"], "other_allele": e["ref"],
            "eaf_out": eaf_out,
        })

import csv
out = "/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/01.材料/01.博士材料/10.毕业/01.小论文/04.立项/KOA 滑膜力学转导基因 0912/04_分析数据/MR_coloc/multi_iv_exposure_outcome.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(f"\n已输出 {len(rows)} 行到 multi_iv_exposure_outcome.csv")
for r in rows:
    print(f"  {r['gene']:9s} {r['SNP']:14s} eQTL p={r['p_exp']:.2e} beta_exp={r['beta_exp']:.3f} beta_out={r['beta_out']:.3f}")
