# -*- coding: utf-8 -*-
"""工具变量间的实际连锁不平衡（LD）核查。

背景：§2.5 用「物理距离 ≥100 kb」代表近似独立。物理距离不是 LD 的替代品——
欧洲人群存在长程 LD 区块，相距数百 kb 仍可能强相关。相关工具会使 IVW 的标准误
被低估、P 值被高估。本脚本从 1000 Genomes 的 GRCh38 相位数据算出工具两两的 r²。

参考面板：1000 Genomes phase 3（GRCh38 重定位版，n=2548，其中 EUR 503）。
取远程 tabix 随机访问，只下工具所在的小窗口。

输出：ld_instruments.csv（逐对 r²）、ld_summary.txt（人读摘要）。
"""
import csv
import glob
import itertools
import os
import re
import sys
import time

import pysam

BASE = ("/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/"
        "Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/"
        "01.材料/01.博士材料/10.毕业/01.小论文/04.立项/滑膜基质刚度生成与刚度感知 0912")
MR = f"{BASE}/04_分析数据/MR_coloc"
PANEL = os.path.expanduser("~/bio-data/koa-stiffness-ld/integrated_call_samples_v3.20130502.ALL.panel")
VCF_BASE = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
            "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/"
            "ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz")

# 基因 → 染色体（GRCh38），与 gene_coords.py 的 GENES 表一致
CHROM = {"VGLL4": "3", "SDC4": "20", "LOXL1": "15",
         "SERPINE1": "7", "SAV1": "14", "ITGA5": "12"}

# LD 的判定门槛。r²>0.1 已足以让「独立工具」的假设不成立（IVW 的 SE 会明显低估），
# 故以此为「需改用相关工具模型」的界线；0.01~0.1 记为「轻微相关，需在文中交代」。
R2_ACTION, R2_NOTE = 0.10, 0.01


def eur_samples():
    out = []
    with open(PANEL, encoding="utf-8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 3 and row[2] == "EUR":
                out.append(row[0])
    return out


def load_instruments():
    """从三套暴露的逐 SNP 文件读出 (暴露, 基因, rsID, 位置)。"""
    recs = []
    for path in sorted(glob.glob(f"{MR}/multi_iv_exposure_outcome_*.csv")):
        exp = re.search(r"_outcome_(\w+)\.csv", path).group(1)
        for r in csv.DictReader(open(path, encoding="utf-8")):
            recs.append(dict(exposure=exp, gene=r["gene"], rsid=r["SNP"],
                             pos=int(r["pos"])))
    return recs


def fetch_window(chrom, start, end, keep_eur, cache):
    """远程取一个窗口，返回 {位置: 剂量向量}；按窗口缓存避免重复下载。"""
    key = (chrom, start, end)
    if key in cache:
        return cache[key]
    url = VCF_BASE.format(chrom=chrom)
    vcf = pysam.VariantFile(url)
    idx = {s: i for i, s in enumerate(vcf.header.samples)}
    sel = [i for s, i in idx.items() if s in keep_eur]
    out = {}
    t0 = time.time()
    for rec in vcf.fetch(chrom, start, end):
        if len(rec.alleles) != 2 or any(len(a) != 1 for a in rec.alleles):
            continue                      # 只用二等位 SNP，与 eQTL 工具口径一致
        doses = []
        ok = True
        for i in sel:
            gt = rec.samples[i].get("GT")
            if gt is None or None in gt:
                ok = False
                break
            doses.append(float(sum(gt)))
        if ok:
            out[rec.pos] = doses
    vcf.close()
    sys.stderr.write(f"    取回 chr{chrom}:{start}-{end} {len(out)} 个 SNP（{time.time()-t0:.0f}s）\n")
    cache[key] = out
    return out


def r2(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    va = sum((x - ma) ** 2 for x in a)
    vb = sum((x - mb) ** 2 for x in b)
    if va == 0 or vb == 0:
        return float("nan")
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    return (cov ** 2) / (va * vb)


def main():
    eur = set(eur_samples())
    sys.stderr.write(f"EUR 样本 {len(eur)} 个\n")
    recs = load_instruments()

    # 每个基因取一次覆盖全部暴露工具的窗口
    by_gene = {}
    for r in recs:
        by_gene.setdefault(r["gene"], []).append(r)

    rows, cache = [], {}
    for gene, rs in sorted(by_gene.items()):
        chrom = CHROM[gene]
        lo = min(x["pos"] for x in rs) - 5000
        hi = max(x["pos"] for x in rs) + 5000
        sys.stderr.write(f"  {gene} (chr{chrom}, {len(rs)} 个工具记录)\n")
        win = fetch_window(chrom, lo, hi, eur, cache)

        by_exp = {}
        for x in rs:
            by_exp.setdefault(x["exposure"], []).append(x)

        for exp, items in sorted(by_exp.items()):
            if len(items) < 2:
                continue
            items = sorted(items, key=lambda x: x["pos"])
            for a, b in itertools.combinations(items, 2):
                ga, gb = win.get(a["pos"]), win.get(b["pos"])
                if ga is None or gb is None:
                    miss = [x["rsid"] for x, g in ((a, ga), (b, gb)) if g is None]
                    rows.append(dict(exposure=exp, gene=gene, snp1=a["rsid"], snp2=b["rsid"],
                                     dist_kb=round(abs(a["pos"] - b["pos"]) / 1000, 1),
                                     r2="NA", n_eur="0",
                                     note=f"1000G 中未取到: {','.join(miss)}"))
                    continue
                # 两变异在 VCF 里各自独立编码，需取交集样本（此处同窗口同样本，直接算）
                m = r2(ga, gb) if len(ga) == len(gb) else float("nan")
                if m != m:
                    note = "单态，无法计算"
                elif m > R2_ACTION:
                    note = f"⚠️ 相关工具（r²>{R2_ACTION}），须改用相关工具模型"
                elif m > R2_NOTE:
                    note = "轻微相关，需在文中交代"
                else:
                    note = "近似独立"
                rows.append(dict(exposure=exp, gene=gene, snp1=a["rsid"], snp2=b["rsid"],
                                 dist_kb=round(abs(a["pos"] - b["pos"]) / 1000, 1),
                                 r2=round(m, 4) if m == m else "NA",
                                 n_eur=len(ga), note=note))

    out_csv = f"{MR}/ld_instruments.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["exposure", "gene", "snp1", "snp2",
                                           "dist_kb", "r2", "n_eur", "note"])
        w.writeheader()
        w.writerows(rows)

    lines = ["工具变量间实际 LD（1000 Genomes phase 3 GRCh38，EUR n=503）", ""]
    flag = [r for r in rows if str(r["r2"]) != "NA" and float(r["r2"]) > R2_ACTION]
    for r in rows:
        lines.append("%-9s %-9s %-14s %-14s %8.1f kb  r²=%-7s %s"
                     % (r["exposure"], r["gene"], r["snp1"], r["snp2"],
                        r["dist_kb"], r["r2"], r["note"]))
    lines += ["", f"需改用相关工具模型的对数: {len(flag)}"]
    open(f"{MR}/ld_summary.txt", "w", encoding="utf-8").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n✓ 已落盘 {out_csv}、{MR}/ld_summary.txt")


if __name__ == "__main__":
    main()
