#!/usr/bin/env python3
# B3 单细胞逐供者定位（§4.3）
#
# 为什么从原始 10x 重建，而不是直接用 oa_synovium_v2.h5ad：
#   该 h5ad 经 h5py 直接检查后确认**只含 obs 组**，没有 X / var / layers / obsm / uns。
#   即它只有细胞级元数据，没有表达矩阵，也没有 cell_type 列（审计报告 F10 记的
#   「obs 含 cell_type、七套正负 marker 打分」与实际不符——实际只有 3 类细胞、6 列打分）。
#   故 §4.3 的所有数字无法从该文件复算，只能回到 NCBI 重新下载的三份 10x 矩阵重建。
#
# 本脚本做四件事，每一件都落盘：
#   ① 逐供者细胞数（QC 前 → QC 后），并与 h5ad 的 obs 做**独立数据路径**的硬对账；
#   ② 逐供者细胞类型比例；
#   ③ 逐供者、逐细胞类型的面板基因表达比例与平均表达；
#   ④ 与正文 §3.6 已写的每一个数字逐条对照（12.80% / 37.48% / 65.45% / 69.1% 等）。
#
# 用途：python3 B3_单细胞逐供者定位.py   （日志由调用方重定向到 03_分析/logs/B3_scRNA.log）

import os, sys, gzip, json, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
import scanpy as sc

warnings.simplefilter("ignore")
sc.settings.verbosity = 1

BASE = "/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/01.材料/01.博士材料/10.毕业/01.小论文/04.立项/滑膜基质刚度生成与刚度感知 0912"
REV  = os.path.join(BASE, "08_修订_2026-09-17")
RAW  = os.path.expanduser("~/bio-data/koa-stiffness-scrna/GSE152805_raw")
H5AD = os.path.expanduser("~/bio-data/koa-stiffness-scrna/oa_synovium_v2.h5ad")
OUT  = os.path.join(REV, "03_分析/outputs/B3_scRNA")
RES  = os.path.join(REV, "02_结果表")
FIG  = os.path.join(REV, "04_图表")
for d in (OUT, RES, FIG):
    os.makedirs(d, exist_ok=True)

PANEL = ["LOXL1", "SERPINE1", "SDC4", "VGLL4", "ITGA5", "SAV1"]
SAMPLES = [("GSM4626763_SY_113", "SY_113"),
           ("GSM4626764_SY_116", "SY_116"),
           ("GSM4626765_SY_118", "SY_118")]
# h5ad obs 的独立对账目标（由 barcodes 交集解码得到，见 B3 结果 md 第一节）
EXPECT_POST = {"SY_113": 5096, "SY_116": 3875, "SY_118": 1629}

# 面板基因在该 CellRanger 参考（32,738 符号 / 32,643 唯一）中均为唯一符号，无需聚合。
# 但参考中**没有 PECAM1 / CD31**，内皮阳性 marker 须改用 CDH5 / VWF / CLDN5 / ENG。
POS = {
    "Fibroblast":  ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "DCN", "LUM", "PDGFRA", "THY1", "FAP"],
    "Macrophage":  ["CD68", "CD14", "LYZ", "AIF1", "CD163", "MRC1", "ITGAM", "FCGR3A"],
    "T cell":      ["CD3D", "CD3E", "CD3G", "CD2", "IL7R", "CD8A", "CD4"],
    "Endothelial": ["CDH5", "VWF", "CLDN5", "ENG", "SELE", "PROX1"],
    "NK":          ["NKG7", "GNLY", "KLRD1", "NCAM1", "KLRF1"],
    "B cell":      ["MS4A1", "CD79A", "CD79B", "CD19"],
    "Mast":        ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
}
NEG = {
    "Fibroblast":  ["PTPRC", "CD68", "CD14", "LYZ", "CD3D", "CD3E", "CDH5", "VWF", "CD19"],
    "Macrophage":  ["CD3D", "CD3E", "COL1A1", "CDH5", "VWF"],
    "T cell":      ["COL1A1", "CD68", "LYZ", "CDH5", "VWF"],
    "Endothelial": ["COL1A1", "CD68", "LYZ", "CD3D", "CD19"],
}
# 正文 §2.6 / §3.6 只用了这 4 类；其余类型用于判断簇归属是否干净，不进入正文比例
MAIN4 = ["Fibroblast", "Macrophage", "T cell", "Endothelial"]
CN = {"Fibroblast": "成纤维细胞", "Macrophage": "巨噬细胞", "T cell": "T 细胞",
      "Endothelial": "内皮细胞", "NK": "NK 细胞", "B cell": "B 细胞", "Mast": "肥大细胞",
      "未定型": "未定型", "Other": "其他"}

def log(*a):
    print(*a, flush=True)

def read_one(stem):
    """读一份 10x：MatrixMarket（基因×细胞）转置成 细胞×基因。"""
    f = os.path.join(RAW, f"{stem}.matrix.mtx.gz")
    with gzip.open(f, "rt") as fh:
        M = mmread(fh).tocsr()                       # genes × cells
    bc = [l.strip() for l in gzip.open(os.path.join(RAW, f"{stem}.barcodes.tsv.gz"), "rt") if l.strip()]
    gn = []
    for l in gzip.open(os.path.join(RAW, f"{stem}.genes.tsv.gz"), "rt"):
        p = l.rstrip("\n").split("\t")
        if len(p) > 1:
            gn.append(p[1])
    assert M.shape == (len(gn), len(bc)), f"{stem} 维度不符: {M.shape} vs {len(gn)}×{len(bc)}"
    A = sc.AnnData(X=sp.csr_matrix(M.T.astype(np.float32)),
                   obs=pd.DataFrame(index=pd.Index(bc, name="barcode")),
                   var=pd.DataFrame(index=pd.Index(gn, name="symbol")))
    # 该参考有 91 个符号重复（涉及 186 行），var 索引非唯一会让 concat 的重索引失败。
    # 加后缀使唯一；面板基因与全部 marker 基因本身都不重复，名字不受影响。
    A.var_names_make_unique()
    return A

def qc(A):
    """§2.6 的质控：每细胞 ≥200 基因、线粒体 reads <20%。返回布尔掩码与指标。"""
    A.var["mt"] = A.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    keep = (A.obs["n_genes_by_counts"] >= 200) & (A.obs["pct_counts_mt"] < 20)
    return keep, A.var["mt"].sum()

# ─────────────────────────────────────────────────────────────────────────────
log("=" * 78)
log("B3 单细胞逐供者定位（§4.3）")
log(f"输入：{RAW}   （2026-09-17 从 NCBI 重新下载并校验的三份 synovium 10x）")
log("=" * 78)

adatas, qc_rows = [], []
for stem, name in SAMPLES:
    A = read_one(stem)
    n_before = A.n_obs
    keep, n_mt = qc(A)
    if name == SAMPLES[0][1]:
        log(f"\n线粒体基因（MT- 前缀）共 {n_mt} 个")
    ob = A.obs
    n_after = int(keep.sum())
    qc_rows.append(dict(
        供者=name, GSM=stem.split("_")[0], QC前细胞数=n_before, QC后细胞数=n_after,
        被剔除=n_before - n_after, 剔除比例=round((n_before - n_after) / n_before * 100, 2),
        基因数下界=int(ob["n_genes_by_counts"].min()), 基因数中位=float(ob["n_genes_by_counts"].median()),
        基因数上界=int(ob["n_genes_by_counts"].max()),
        线粒体比例下界=round(float(ob["pct_counts_mt"].min()), 3),
        线粒体比例中位=round(float(ob["pct_counts_mt"].median()), 3),
        线粒体比例上界=round(float(ob["pct_counts_mt"].max()), 3),
        中位总counts=float(ob["total_counts"].median()),
    ))
    A = A[keep].copy()
    A.obs["donor"] = name
    adatas.append((name, A, stem))

df_qc = pd.DataFrame(qc_rows)
log("\n=== ① 逐供者细胞数（QC 前 → QC 后） ===")
log(df_qc.to_string(index=False))
df_qc.to_csv(os.path.join(OUT, "B3_逐供者细胞数.csv"), index=False, encoding="utf-8-sig")

log("\n=== ①b 与 h5ad obs 的独立对账（硬判据） ===")
chk = []
for name, A, stem in adatas:
    exp = EXPECT_POST[name]
    got = A.n_obs
    chk.append(dict(供者=name, h5ad期望=exp, 本次重建=got, 差=got - exp, 一致=("✓" if got == exp else "✗")))
df_chk = pd.DataFrame(chk)
log(df_chk.to_string(index=False))
df_chk.to_csv(os.path.join(OUT, "B3_与h5ad对账.csv"), index=False, encoding="utf-8-sig")
if (df_chk["差"] != 0).any():
    log("\n!!! 警告：本次重建的 QC 后细胞数与 h5ad 不一致。")
    log("    若不一致，则本脚本的细胞类型划分与正文原分析不是同一套细胞，")
    log("    下方所有比例均为「本文口径」，与正文数字的差异须按口径差异解读，不能当作复现失败。")
    log("    （脚本不中断：后续结果仍有独立价值，但必须在结果 md 中如实标注。）")
else:
    log("\n✓ 逐供者 QC 后细胞数与 h5ad 完全一致 —— 两条独立数据路径互相印证。")

# ─────────────────────────────────────────────────────────────────────────────
log("\n=== ② 逐样本 Scrublet 双细胞检测（§2.6 声称已做，此处实测） ===")
try:
    import scrublet as _scr
    HAS_SCR = True
except Exception as e:                                     # noqa: BLE001
    HAS_SCR = False
    log(f"!! scrublet 不可用：{e}")

dbl_rows = []
for name, A, stem in adatas:
    if not HAS_SCR:
        dbl_rows.append(dict(供者=name, 状态="scrublet 不可用", 检出双细胞=None))
        continue
    B = A.copy()
    try:
        sc.pp.scrublet(B, random_state=0)
        n = int(B.obs["predicted_doublet"].sum())
        dbl_rows.append(dict(供者=name, 状态="已运行", 检出双细胞=n,
                             双细胞比例=round(n / B.n_obs * 100, 2),
                             中位doublet_score=round(float(B.obs["doublet_score"].median()), 4),
                             最大doublet_score=round(float(B.obs["doublet_score"].max()), 4)))
    except Exception as e:                                 # noqa: BLE001
        dbl_rows.append(dict(供者=name, 状态=f"运行失败: {e}", 检出双细胞=None))
df_dbl = pd.DataFrame(dbl_rows)
log(df_dbl.to_string(index=False))
df_dbl.to_csv(os.path.join(OUT, "B3_Scrublet逐供者.csv"), index=False, encoding="utf-8-sig")

# ─────────────────────────────────────────────────────────────────────────────
log("\n=== ③ 合并 → 归一化 → 降维聚类（照 §2.6 的设定） ===")
ad = sc.concat([a for _, a, _ in adatas], index_unique=None, join="outer", fill_value=0)
ad.obs_names_make_unique()
log(f"合并后：{ad.n_obs} 细胞 × {ad.n_vars} 基因")
log(f"逐供者：{dict(ad.obs['donor'].value_counts())}")

# 原始 counts 必须在 normalize_total 之前另存一层：正文的「平均表达」口径未知，
# 第七节要拿真 counts 反推它到底用的是哪一种均值。
ad.layers["counts"] = ad.X.copy()
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)
ad.layers["lognorm"] = ad.X.copy()        # 面板表达一律在这个未 scale 的层上算
_c = ad.layers["counts"].data
assert np.allclose(_c, np.round(_c)), \
    "counts 层出现非整数——已被 log1p 改写，第七节的均值口径排查将失效"
log(f"counts 层完好：非零值 {_c.min():.0f}–{_c.max():.0f}，全为整数 ✓")
full = ad.copy()                          # 保留全基因，供 marker / 面板统计

sc.pp.highly_variable_genes(ad, n_top_genes=2000)         # §2.6：2000 高变基因
log(f"高变基因：{int(ad.var['highly_variable'].sum())}")
ad = ad[:, ad.var["highly_variable"]].copy()
sc.pp.scale(ad, max_value=10)
sc.tl.pca(ad, n_comps=15, svd_solver="arpack")           # §2.6：15 个主成分
sc.pp.neighbors(ad, n_neighbors=15, n_pcs=15)
sc.tl.leiden(ad, resolution=0.5, key_added="leiden",
             flavor="igraph", n_iterations=2, directed=False)   # §2.6：Leiden res 0.5
n_clu = ad.obs["leiden"].nunique()
log(f"Leiden(res=0.5) 得到 {n_clu} 个簇  —— h5ad 记为 14 个")
# 用 reindex 对齐而不是直接 .values —— 细胞顺序一旦不同，靠位置赋值会静默错配
full.obs["leiden"] = ad.obs["leiden"].astype(str).reindex(full.obs_names)
assert full.obs["leiden"].notna().all(), "簇标签未能对齐到全部细胞"

# ─────────────────────────────────────────────────────────────────────────────
log("\n=== ④ 细胞类型注释证据（阳性 + 阴性，两个方向都出表） ===")
for ct, genes in POS.items():
    g = [x for x in genes if x in full.var_names]
    sc.tl.score_genes(full, g, score_name=f"score_{ct}", random_state=0)

# 簇 × 类型 打分均值
score_cols = [f"score_{ct}" for ct in POS]
df_score = full.obs.groupby("leiden", observed=True)[score_cols].mean()
df_score.columns = [c.replace("score_", "") for c in df_score.columns]
df_score = df_score.round(3)
log("\n-- 4.1 簇 × 细胞类型 打分矩阵（score_genes 背景校正后）--")
log(df_score.to_string())
df_score.to_csv(os.path.join(OUT, "B3_簇类型打分.csv"), encoding="utf-8-sig")

# 簇 × marker 基因 平均表达（阳性与阴性同表，便于直接看方向）
allmk = sorted(set(sum(POS.values(), []) + sum(NEG.values(), [])))
allmk = [g for g in allmk if g in full.var_names]
df_mk = pd.DataFrame(np.asarray(full[:, allmk].layers["lognorm"].todense()),
                     index=full.obs_names, columns=allmk)
df_mk["leiden"] = full.obs["leiden"].values
mk_mean = df_mk.groupby("leiden", observed=True).mean().round(3)
df_mk.to_csv(os.path.join(OUT, "B3_逐细胞marker表达.csv.gz"), compression="gzip")
mk_mean.to_csv(os.path.join(OUT, "B3_簇marker均值.csv"), encoding="utf-8-sig")
log(f"\n-- 4.2 簇 × marker 平均表达（{len(allmk)} 个基因，含阳性与阴性）--")
log(mk_mean.to_string())

# 归属：簇打分的 argmax，且获胜打分须 > 0，否则记为未定型
assign = {}
for cl in df_score.index:
    s = df_score.loc[cl]
    win = s.idxmax()
    assign[cl] = win if s[win] > 0 else "未定型"
full.obs["cell_type"] = full.obs["leiden"].map(assign).astype(str)
df_assign = pd.DataFrame({"簇": list(assign.keys()),
                          "归属": [CN.get(v, v) for v in assign.values()],
                          "获胜打分": [round(df_score.loc[k, v], 3) if v != "未定型" else None
                                       for k, v in assign.items()],
                          "簇细胞数": full.obs["leiden"].value_counts().reindex(list(assign.keys())).values})
df_assign.to_csv(os.path.join(OUT, "B3_簇归属.csv"), index=False, encoding="utf-8-sig")
log("\n-- 4.3 簇 → 细胞类型 归属 --")
log(df_assign.to_string(index=False))

# ─────────────────────────────────────────────────────────────────────────────
log("\n=== ⑤ 细胞类型组成（逐供者） ===")
tab = pd.crosstab(full.obs["donor"], full.obs["cell_type"])
tab = tab.reindex(columns=[c for c in ["Fibroblast", "Macrophage", "T cell", "Endothelial", "NK", "B cell", "Mast", "未定型"] if c in tab.columns])
log("\n-- 5.1 细胞数 --")
log(tab.to_string())
pct = (tab.div(tab.sum(axis=1), axis=0) * 100).round(2)
log("\n-- 5.2 比例（分母 = 该供者全部细胞，含未定型）--")
log(pct.to_string())
main4 = [c for c in MAIN4 if c in tab.columns]
pct4 = (tab[main4].div(tab[main4].sum(axis=1), axis=0) * 100).round(2)
log("\n-- 5.3 比例（分母 = 该供者四类细胞之和，即正文口径）--")
log(pct4.to_string())
tab.to_csv(os.path.join(OUT, "B3_逐供者细胞类型_细胞数.csv"), encoding="utf-8-sig")
pct.to_csv(os.path.join(OUT, "B3_逐供者细胞类型_比例.csv"), encoding="utf-8-sig")
pct4.to_csv(os.path.join(OUT, "B3_逐供者细胞类型_比例_四类口径.csv"), encoding="utf-8-sig")
log("\n（全部细胞口径：")
log((tab.sum(axis=0).to_frame("细胞数").T.to_string()) + "）")

# ─────────────────────────────────────────────────────────────────────────────
log("\n=== ⑥ 面板基因表达：逐细胞类型 + 逐供者 ===")
# 在 lognorm 层上取面板基因
pg = [g for g in PANEL if g in full.var_names]
missing = [g for g in PANEL if g not in full.var_names]
if missing:
    log(f"!! 面板基因在表达矩阵中缺失：{missing}")
Xp = pd.DataFrame(np.asarray(full[:, pg].layers["lognorm"].todense()),
                  index=full.obs_names, columns=pg)           # log1p(CP10K)
CNT = pd.DataFrame(np.asarray(full[:, pg].layers["counts"].todense()),
                   index=full.obs_names, columns=pg)          # 原始 counts
TOT = full.obs["total_counts"].astype(float)
ct_of = full.obs["cell_type"]
XCP = CNT.div(TOT, axis=0) * 1e4                              # CP10K（线性，未取对数）
XCPM = XCP * 100.0                                            # CPM
meta = full.obs[["donor", "cell_type", "leiden"]].copy()

def stats_expr(sub_meta, by):
    """按 by 分组。同时给出 log1p 与线性两种尺度：正文的「平均表达」用的是线性 CP10K，
    若只给 log1p 会把它误判成数值不一致（见 ⑥b 的口径排查）。"""
    out = []
    for key, g in sub_meta.groupby(by, observed=True):
        if not isinstance(key, tuple):
            key = (key,)
        for gene in pg:
            v_log = Xp.loc[g.index, gene].values
            v_cp = XCP.loc[g.index, gene].values
            v_cpm = XCPM.loc[g.index, gene].values
            nz = v_cp > 0
            out.append(dict(zip(by, key)) | {
                "基因": gene, "细胞数": len(g), "非零细胞数": int(nz.sum()),
                "表达比例": round(nz.mean() * 100, 2),
                "平均log1p": round(float(v_log.mean()), 4),
                "平均CP10K": round(float(v_cp.mean()), 4),
                "非零细胞平均CPM": round(float(v_cpm[nz].mean()), 2) if nz.any() else 0.0,
                "组内最大CPM": round(float(v_cpm.max()), 2),
            })
    return pd.DataFrame(out)

by_ct = stats_expr(meta, ["cell_type"])
by_dc = stats_expr(meta, ["donor", "cell_type"])
by_ct.to_csv(os.path.join(OUT, "B3_面板表达_逐细胞类型.csv"), index=False, encoding="utf-8-sig")
by_dc.to_csv(os.path.join(OUT, "B3_面板表达_逐供者逐细胞类型.csv"), index=False, encoding="utf-8-sig")

log("\n-- 6.1 逐细胞类型：表达比例（%） --")
log(by_ct.pivot(index="cell_type", columns="基因", values="表达比例").to_string())
log("\n-- 6.2a 逐细胞类型：平均 log1p(CP10K)（作图与可视化用）--")
log(by_ct.pivot(index="cell_type", columns="基因", values="平均log1p").to_string())
log("\n-- 6.2b 逐细胞类型：平均 CP10K（**正文用的就是这个线性口径**，见 ⑥b）--")
log(by_ct.pivot(index="cell_type", columns="基因", values="平均CP10K").to_string())
log("\n-- 6.3 逐细胞类型：组内最大 CPM（低于 20 者按「未达可靠检出」记，见方法说明）--")
mx = by_ct.pivot(index="cell_type", columns="基因", values="组内最大CPM")
log(mx.to_string())
log("\n-- 6.4 逐供者 × 细胞类型：表达比例（%） --")
p_dc = by_dc.pivot_table(index=["donor", "cell_type"], columns="基因", values="表达比例")
log(p_dc.to_string())

# ─────────────────────────────────────────────────────────────────────────────
log("\n=== ⑥b 「平均表达」口径排查 ===")
log("正文给了 4 个平均表达值（LOXL1 成纤维 0.15；SERPINE1 内皮 5.32／成纤维 0.90／巨噬 0.49），")
log("但没说用哪一种平均值。而「表达比例」逐条对得上、「平均值」逐条偏高 1.6–4.2 倍，")
log("说明差异出在口径而不是细胞划分。此处把 10 种常见口径全部实测，看正文数对应哪一种。")

def all_cands(idx, gene):
    c = CNT.loc[idx, gene].values.astype(float)
    t = TOT.loc[idx].values
    p10 = np.divide(c, t, out=np.zeros_like(c), where=t > 0) * 1e4
    pm = p10 * 100.0
    m = c > 0
    def mu(v):
        return float(np.mean(v)) if len(v) else np.nan
    return {
        "①全细胞 mean log1p(CP10K)": mu(np.log1p(p10)),
        "②非零细胞 mean log1p(CP10K)": mu(np.log1p(p10[m])),
        "③全细胞 mean log1p(原始count)": mu(np.log1p(c)),
        "④非零细胞 mean log1p(原始count)": mu(np.log1p(c[m])),
        "⑤全细胞 mean CP10K": mu(p10),
        "⑥非零细胞 mean CP10K": mu(p10[m]),
        "⑦全细胞 mean CPM": mu(pm),
        "⑧非零细胞 mean CPM": mu(pm[m]),
        "⑨全细胞 mean log2(CPM+1)": mu(np.log2(pm + 1)),
        "⑩非零细胞 mean log2(CPM+1)": mu(np.log2(pm[m] + 1)),
    }

CLAIM_MEAN = [("LOXL1", "Fibroblast", 0.15, "讨论：LOXL1 成纤维平均 0.15"),
              ("SERPINE1", "Endothelial", 5.32, "讨论：SERPINE1 内皮平均 5.32"),
              ("SERPINE1", "Fibroblast", 0.90, "讨论：SERPINE1 成纤维平均 0.90"),
              ("SERPINE1", "Macrophage", 0.49, "讨论：SERPINE1 巨噬平均 0.49")]
rows = []
for gene, ct, claim, src in CLAIM_MEAN:
    idx = full.obs_names[ct_of == ct]
    d = all_cands(idx, gene)
    best = max(d.items(), key=lambda kv: -abs(kv[1] - claim) if np.isfinite(kv[1]) else -1e9)
    rows.append({"基因": gene, "细胞类型": CN.get(ct, ct), "正文值": claim, **d,
                 "最接近的口径": best[0], "该口径值": round(best[1], 4), "出处": src})
df_cand = pd.DataFrame(rows)
pd.set_option("display.width", 250, "display.max_columns", 50)
log("\n" + df_cand.drop(columns=["最接近的口径", "该口径值"]).to_string(index=False))
log("\n-- 每个正文值最接近哪一种口径 --")
log(df_cand[["基因", "细胞类型", "正文值", "最接近的口径", "该口径值", "出处"]].to_string(index=False))
df_cand.to_csv(os.path.join(OUT, "B3_平均表达口径排查.csv"), index=False, encoding="utf-8-sig")

# 正文的 4 个值是否指向同一个口径
uniq = df_cand["最接近的口径"].nunique()
log(f"\n4 个正文值指向 {uniq} 种不同口径。")
if uniq == 1:
    log(f"→ 可以判定：正文的「平均表达」= {df_cand['最接近的口径'].iloc[0]}")
else:
    log("→ 4 个值并不指向同一种口径，说明正文的平均值不能由本数据以任何单一口径复现；")
    log("   须在结果 md 中记为「口径不明、无法复现」，不得据此改写正文数值。")

# ─────────────────────────────────────────────────────────────────────────────
log("\n=== ⑦ 与正文 §3.6 / 摘要 / 讨论 逐条对照 ===")
CMP = [
    ("LOXL1", "Fibroblast",  "表达比例", 12.80, "摘要与 §3.6：成纤维细胞 12.8%"),
    ("LOXL1", "Fibroblast",  "平均CP10K", 0.15, "讨论：LOXL1 平均 0.15"),
    ("SDC4",  "Fibroblast",  "表达比例", 37.48, "§3.6：成纤维细胞 37.48%"),
    ("SDC4",  "Macrophage",  "表达比例", 17.71, "§3.6：巨噬细胞 17.71%"),
    ("VGLL4", "Endothelial", "表达比例", 65.45, "§3.6：内皮 65.45%"),
    ("VGLL4", "Fibroblast",  "表达比例", 45.31, "§3.6：成纤维 45.31%"),
    ("SERPINE1", "Endothelial", "表达比例", 69.1, "讨论：内皮 69.1%"),
    ("SERPINE1", "Endothelial", "平均CP10K", 5.32, "讨论：内皮平均 5.32"),
    ("SERPINE1", "Fibroblast",  "表达比例", 39.6, "讨论：成纤维 39.6%"),
    ("SERPINE1", "Fibroblast",  "平均CP10K", 0.90, "讨论：成纤维平均 0.90"),
    ("SERPINE1", "Macrophage",  "表达比例", 21.4, "讨论：巨噬 21.4%"),
    ("SERPINE1", "Macrophage",  "平均CP10K", 0.49, "讨论：巨噬平均 0.49"),
    ("SERPINE1", "T cell",      "表达比例", 4.4,  "讨论：T 细胞 4.4%"),
]
rows = []
bt = by_ct.set_index(["cell_type", "基因"])
for gene, ct, metric, claim, src in CMP:
    key = (ct, gene)
    if key not in bt.index or metric not in bt.columns:
        rows.append(dict(基因=gene, 细胞类型=CN.get(ct, ct), 指标=metric,
                         正文值=claim, 本次重算=None, 差=None, 出处=src, 判定="无法计算"))
        continue
    got = float(bt.loc[key, metric])
    rows.append(dict(基因=gene, 细胞类型=CN.get(ct, ct), 指标=metric,
                     正文值=claim, 本次重算=round(got, 3) if abs(got) >= 0.001 else got,
                     差=round(got - claim, 3), 出处=src,
                     判定=("一致" if abs(got - claim) <= max(0.05 * abs(claim), 0.05) else "不一致")))
df_cmp = pd.DataFrame(rows)
log(df_cmp.to_string(index=False))
df_cmp.to_csv(os.path.join(OUT, "B3_与正文对照.csv"), index=False, encoding="utf-8-sig")

ok = (df_cmp["判定"] == "一致").sum()
log(f"\n一致 {ok} / {len(df_cmp)} 条；不一致 {int((df_cmp['判定'] == '不一致').sum())} 条；"
    f"无法计算 {int((df_cmp['判定'] == '无法计算').sum())} 条")

# 不等式型主张必须单列：它们是「其余细胞类型全部 < X」的断言，逐值表查不出来
log("\n-- 7b 不等式型主张核验（摘要与讨论的「其余 <X」）--")
OTHER4 = ["Macrophage", "T cell", "Endothelial", "Mast"]
piv_p = by_ct.pivot(index="cell_type", columns="基因", values="表达比例")
piv_cp = by_ct.pivot(index="cell_type", columns="基因", values="平均CP10K")
IMMUNE = ["Macrophage", "T cell"]
ineq = []
for gene, val, grp, gname, scale, src in [
        ("LOXL1", 1.6,   OTHER4, "其余（非成纤维）", "比例", "摘要与 §3.6：LOXL1「其余均 <1.6%」"),
        ("LOXL1", 1.0,   IMMUNE, "免疫（巨噬+T）",   "比例", "§3.6：LOXL1「免疫 <1%」"),
        ("VGLL4", 13.0,  IMMUNE, "免疫（巨噬+T）",   "比例", "§3.6：VGLL4「免疫 ≤13%」"),
        ("LOXL1", 0.022, OTHER4, "其余（非成纤维）", "平均CP10K", "讨论：LOXL1「其余均 <0.022」")]:
    tab = piv_cp if scale == "平均CP10K" else piv_p
    col = tab[gene].reindex([g for g in grp if g in tab.index])
    worst = col.idxmax()
    ineq.append({"基因": gene, "正文断言": f"{gname} 均 < {val}（{scale}）", "出处": src,
                 "实测最大者": CN.get(worst, worst), "实测最大值": round(float(col.max()), 4),
                 "是否成立": "成立" if float(col.max()) < val else "**不成立**",
                 "逐类型值": "; ".join(f"{CN.get(k,k)} {v:.3f}" for k, v in col.items())})
df_ineq = pd.DataFrame(ineq)
log(df_ineq[["基因", "正文断言", "实测最大者", "实测最大值", "是否成立", "出处"]].to_string(index=False))
for r in ineq:
    log(f"    {r['正文断言']} → {r['逐类型值']}")
df_ineq.to_csv(os.path.join(OUT, "B3_不等式型主张核验.csv"), index=False, encoding="utf-8-sig")

# 小细胞群里，「表达比例」会被少量细胞和环境污染主导。逐供者一看就露馅的数字，
# 不能直接当成「正文写错了」——必须先问它在最干净的细胞上是否还站得住。
log("\n-- 7c 小样本敏感性：只留标志基因明确的细胞 --")
CONF = {"T cell": ("CD3D", "PTPRC"), "Macrophage": ("LYZ", "CD14"), "Fibroblast": ("DCN", "LUM")}
sens = []
for ct, (m1, m2) in CONF.items():
    idx = full.obs_names[ct_of == ct]
    conf_idx = idx[(full[idx, m1].layers["lognorm"].toarray().ravel() > 0) &
                   (full[idx, m2].layers["lognorm"].toarray().ravel() > 0)]
    for gene in ["LOXL1", "SERPINE1"]:
        allv = XCP.loc[idx, gene].values
        cv = XCP.loc[conf_idx, gene].values
        pd_ = {}
        for d in ["SY_113", "SY_116", "SY_118"]:
            sub = idx[full.obs.loc[idx, "donor"] == d]
            pd_[d] = (round(float((XCP.loc[sub, gene].values > 0).mean() * 100), 2) if len(sub) else None)
        sens.append({"细胞类型": CN.get(ct, ct), "基因": gene,
                     "全部细胞_n": len(idx), "全部细胞_表达比例": round(float((allv > 0).mean() * 100), 2),
                     "全部细胞_平均CP10K": round(float(allv.mean()), 4),
                     "明确细胞_n": len(conf_idx),
                     "明确细胞_表达比例": round(float((cv > 0).mean() * 100), 2) if len(cv) else None,
                     "明确细胞_平均CP10K": round(float(cv.mean()), 4) if len(cv) else None,
                     "逐供者_表达比例": "; ".join(f"{k} {v}%" for k, v in pd_.items())})
df_sens = pd.DataFrame(sens)
log(df_sens.to_string(index=False))
df_sens.to_csv(os.path.join(OUT, "B3_小样本敏感性.csv"), index=False, encoding="utf-8-sig")

# ─────────────────────────────────────────────────────────────────────────────
log("\n=== ⑧ 出图：已移出本脚本 ===")
log("图 5 由 B3b_图5_共享色阶重绘_2026-09-28.py 从本脚本落盘的冻结源表画出，")
log("本脚本不再出图。原因与旧图的三处缺陷：")
log("  ① 旧版四个面板各自独立取色阶，色条却只由最后一个面板的散点生成 —— 读者会以为同尺；")
log("  ② 旧版首面板「All donors」是三供者比例等权平均（SERPINE1 内皮 74.60%），")
log("     而正文与图注引的是合并细胞值（68.31%）—— 两个口径在图上被混为一谈；")
log("  ③ 旧版 suptitle 写 \"dot size = % expressing\" 而图上没有面积图例。")
log("保留两个出图路径会让两份代码各自漂移，故此处只留指路，不复制绘图逻辑。")
log("源表：B3_面板表达_逐细胞类型.csv / B3_面板表达_逐供者逐细胞类型.csv（本脚本 ⑥ 段产出）")


# 逐细胞元数据落盘（供复核，不只留汇总）
full.obs[["donor", "leiden", "cell_type", "n_genes_by_counts", "pct_counts_mt", "total_counts"]] \
    .to_csv(os.path.join(OUT, "B3_逐细胞元数据.csv.gz"), compression="gzip")

summary = dict(
    细胞数_QC后=int(full.n_obs),
    逐供者细胞数={k: int(v) for k, v in full.obs["donor"].value_counts().items()},
    Leiden簇数=int(n_clu),
    细胞类型_逐供者={k: int(v) for k, v in full.obs["cell_type"].value_counts().items()},
    面板基因=pg, 面板缺失基因=missing,
    与h5ad对账一致=bool((df_chk["差"] == 0).all()),
    与正文一致条数=int(ok), 与正文对照总条数=int(len(df_cmp)),
)
with open(os.path.join(OUT, "B3_summary.json"), "w") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=2)
log("\n" + json.dumps(summary, ensure_ascii=False, indent=2))
log("\n产出目录：" + OUT)
log("B3 完成。")
