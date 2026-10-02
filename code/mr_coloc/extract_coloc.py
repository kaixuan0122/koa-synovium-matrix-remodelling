# 提取 coloc 输入：eQTL(GTEx QTD000001) + KOA(FinnGen) 的基因 cis 区域 summary
import pysam, gzip, json

url = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000001/QTD000001/QTD000001.all.tsv.gz"
tb = pysam.TabixFile(url)

genes = [
  ("VGLL4", "ENSG00000144560", "3",  11560000),
  ("CYR61", "ENSG00000142871", "1",  86046000),
  ("SDC4",  "ENSG00000124145", "20", 43952000),
]

finngen = "/tmp/finngen_KOA.gz"

def get_finngen_region(chrom, lo, hi):
    rows = []
    with gzip.open(finngen, "rt") as f:
        for line in f:
            p = line.split("\t")
            if p[0] == chrom and lo <= int(p[1]) <= hi:
                # chrom pos ref alt rsids ... pval(beta) sebeta af_alt
                rows.append({"rsid": p[4], "pos": int(p[1]), "ref": p[2], "alt": p[3],
                             "beta": float(p[8]), "se": float(p[9]), "maf": float(p[10])})
    return rows

for sym, ensg, chrom, pos in genes:
    lo, hi = pos - 500000, pos + 500000
    # eQTL
    eqtl = {}
    for r in tb.fetch(chrom, lo, hi):
        f = r.split("\t")
        if f[0] == ensg:
            eqtl[f[18]] = {"pos": int(f[2]), "ref": f[3], "alt": f[4], "maf": float(f[7]),
                           "beta": float(f[9]), "se": float(f[10]), "p": float(f[8])}
    # KOA
    koa = get_finngen_region(chrom, lo, hi)
    koa_map = {k["rsid"]: k for k in koa}
    # 交集 + 对齐等位基因
    common = sorted(set(eqtl.keys()) & set(koa_map.keys()))
    out = []
    for rs in common:
        e = eqtl[rs]; k = koa_map[rs]
        # 等位基因对齐：若 ref/alt 相反，翻转 KOA beta
        if e["ref"] == k["ref"] and e["alt"] == k["alt"]:
            beta_koa = k["beta"]
        elif e["ref"] == k["alt"] and e["alt"] == k["ref"]:
            beta_koa = -k["beta"]
        else:
            continue
        out.append({"rsid": rs, "pos": e["pos"],
                    "beta_eqtl": e["beta"], "se_eqtl": e["se"], "maf_eqtl": e["maf"],
                    "beta_koa": beta_koa, "se_koa": k["se"], "maf_koa": k["maf"]})
    with open(f"/tmp/merge_data/coloc_{sym}.json", "w") as f:
        json.dump(out, f)
    print(f"{sym}: eQTL {len(eqtl)} 个, KOA {len(koa)} 个, 交集 {len(out)} 个 SNP")
