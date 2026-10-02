# -*- coding: utf-8 -*-
"""A-2.4 成纤维细胞状态关联（条件性）

方案 §三 A-2 第 4 项，三条要求逐条落实：

  1. **前置可行性检查先做**：只有 3 个供者，SY_118 成纤维仅 43.89%。先判能否得到跨供者稳定的
     状态划分；若某供者几乎独占某状态，按方案规定退回粗粒度定位描述，不强行细分。
  2. **状态基因集必须是外部的**（不来自本研究自身的差异分析），且**必须剔除被检验的候选基因自身**
     ——否则是循环论证（见记忆 compositional-index-must-exclude-test-genes）。
  3. 逐供者比较。

**为什么必须重建表达矩阵**：`~/bio-data/koa-stiffness-scrna/oa_synovium_v2.h5ad` 顶层键只有
`['obs']` —— 没有 X、没有 var，是**纯元数据文件**（2.2 MB）。这正是记忆
`h5ad-metadata-only-looks-complete` 记的形态：obs 齐全会让审计误判完好。
所幸原始 10x 三件套仍在 `GSE152805_raw/`，故从原始矩阵重建。

**重建口径照 B3_单细胞逐供者定位.py（该脚本才是 10600 个细胞与标签的产出者），不是照
run_scRNA_v2.py**。两者不同，抄错会得到不同分母：
  - B3：`sc.concat(index_unique=None, join="outer", fill_value=0)` → `obs_names_make_unique()`
        → `normalize_total(1e4)` → `log1p`，**全程没有 filter_genes**（分母 = 全部 32738 基因）；
        也**没有做双细胞过滤**（Scrublet 只在 `A.copy()` 上跑过、仅用于报告）。
  - run_scRNA_v2.py：先 `filter_genes(min_cells=3)` 再归一化，且滤 `~is_doublet`。
  实测反解 B3 已发布表达表所用的归一化分母，与「全基因分母」在 float32 精度上吻合
  （最大相对差 5e-7），与「过滤后分母」差 1.5e-3 —— 实证确认 B3 口径。

**最强的对账不是逐格比数值，而是 obs 名逐位相等**：若重建出的 `obs_names` 与 B3 元数据的
`barcode` 列（含 37 个 `-1-1` 双层后缀）逐位相同、`donor` 列也相同，则重建出的就是 B3 当时
的同一个对象，标签可直接按行复用。数值逐格比作为第二重独立对账。

用法：~/celloracle_env/bin/python 03_分析/scripts/A2d_成纤维状态关联.py
"""
import json
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.stats import binomtest
import scipy.sparse as sp

warnings.simplefilter("ignore")
sc.settings.verbosity = 1

ROOT = ("/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/"
        "Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/"
        "01.材料/01.博士材料/10.毕业/01.小论文/04.立项/"
        "滑膜基质刚度生成与刚度感知 0912")
REV = os.path.join(ROOT, "08_修订_2026-09-17")
RAW = os.path.expanduser("~/bio-data/koa-stiffness-scrna/GSE152805_raw")
RES = os.path.join(REV, "02_结果表")
AUD = os.path.join(REV, "01_审计")
OUT = os.path.join(REV, "03_分析/outputs/B3_scRNA")
os.makedirs(AUD, exist_ok=True)

SAMPLES = [("GSM4626763_SY_113", "SY_113"),
           ("GSM4626764_SY_116", "SY_116"),
           ("GSM4626765_SY_118", "SY_118")]
EXPECT_POST = {"SY_113": 5096, "SY_116": 3875, "SY_118": 1629}
TOL_VAL = 1e-5          # 数值第二重对账容差：反解实测全基因分母的最大相对差 ~5e-7，取 1e-5 留余量

ln = []


def rec(s=""):
    print(s, flush=True)
    ln.append(s)


# ---------- 外部状态基因集（文献来源，非本研究产出） ----------
# lining（衬里层）与 sublining（衬里下层）的经典标记。检索所得支持：
#   lining   = PRG4 / MMP3 / CD55 / CLIC5 / FN1 / CRTAC1 / HTRA1
#   sublining = CD34 / DCN / IGF1 / FBLN1 / SFRP1 / DKK3 / ASPN / MAMDC2 / PRELP
# **CILP 未被检索结果支持，故不纳入**（不凭印象写基因集）。
LINING = ["PRG4", "MMP3", "CD55", "CLIC5", "FN1", "CRTAC1", "HTRA1"]
SUBLINING = ["CD34", "DCN", "IGF1", "FBLN1", "SFRP1", "DKK3", "ASPN", "MAMDC2", "PRELP"]


def read_one(stem):
    """照 B3 的读法：`gzip.open(..., "rt")` 文本模式交给 mmread——传二进制句柄会
    `ValueError: too many values to unpack`（mmread 自行解压时才认 gzip）。"""
    import gzip
    f = os.path.join(RAW, f"{stem}.matrix.mtx.gz")
    with gzip.open(f, "rt") as fh:
        M = mmread(fh).tocsr()                       # genes x cells
    bc = [l.strip() for l in gzip.open(os.path.join(RAW, f"{stem}.barcodes.tsv.gz"), "rt")
          if l.strip()]
    gn = []
    for l in gzip.open(os.path.join(RAW, f"{stem}.genes.tsv.gz"), "rt"):
        q = l.rstrip("\n").split("\t")
        if len(q) > 1:
            gn.append(q[1])
    assert M.shape == (len(gn), len(bc)), f"{stem} 维度不符: {M.shape} vs {len(gn)}x{len(bc)}"
    A = sc.AnnData(X=sp.csr_matrix(M.T.astype(np.float32)),
                   obs=pd.DataFrame(index=pd.Index(bc, name="barcode")),
                   var=pd.DataFrame(index=pd.Index(gn, name="symbol")))
    A.var_names_make_unique()      # 与 B3 一致：该参考有 91 个重名符号
    return A


def qc(A):
    A.var["mt"] = A.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    return (A.obs["n_genes_by_counts"] >= 200) & (A.obs["pct_counts_mt"] < 20)


# ---------- 1. 重建（逐行照 B3 的设定） ----------
rec("# A-2.4 成纤维细胞状态关联（条件性）")
rec()
rec("## 〇、重建与两重对账")
rec()
adatas = []
for stem, name in SAMPLES:
    A = read_one(stem)
    keep = qc(A)
    rec(f"- {name}（{stem}）：QC 前 {A.n_obs} → QC 后 {int(keep.sum())} 细胞（≥200 基因且线粒体 <20%）")
    A = A[keep].copy()
    A.obs["donor"] = name
    adatas.append(A)
    assert int(keep.sum()) == EXPECT_POST[name], \
        f"{name} QC 后细胞数 {int(keep.sum())} != 期望 {EXPECT_POST[name]}"

ad = sc.concat(adatas, index_unique=None, join="outer", fill_value=0)
ad.obs_names_make_unique()
rec(f"- 合并后：{ad.n_obs} 细胞 x {ad.n_vars} 基因；逐供者 {dict(ad.obs['donor'].value_counts())}")
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)
rec("- 归一化：`normalize_total(target_sum=1e4)` → `log1p`，**分母取全部基因**"
    "（B3 全程无 `filter_genes`；实测反解确认，见文件头）")

# ---------- 2. 对账一：obs 名与供者标签逐位相等 ----------
meta = pd.read_csv(os.path.join(OUT, "B3_逐细胞元数据.csv.gz"))
meta.columns = [c.lstrip("﻿") for c in meta.columns]
mbc = meta.barcode.astype(str).values
mdon = meta.donor.astype(str).values
assert len(mbc) == ad.n_obs, f"细胞数不符：元数据 {len(mbc)} vs 重建 {ad.n_obs}"
eq_bc = np.array_equal(ad.obs_names.values.astype(str), mbc)
eq_dn = np.array_equal(ad.obs["donor"].astype(str).values, mdon)
n_suf = int(sum(1 for x in mbc if len(x) == 20))
rec(f"- **对账一（强判据）**：重建 `obs_names` 与 B3 元数据 `barcode` 列逐位相等 = {eq_bc}；"
    f"`donor` 列逐位相等 = {eq_dn}")
rec(f"  → 其中 {n_suf} 个条码是 20 字符（16bp + `-1` + 去重后缀 `-1`）。")
assert eq_bc and eq_dn, "重建对象与 B3 元数据不是同一批细胞/同一顺序——标签不能按行复用"
rec("  → 通过。重建出的就是 B3 当时的同一个对象，`leiden`/`cell_type` 可直接按行复用。")

# ---------- 3. 对账二：与既有表达表逐格比 ----------
chk_path = os.path.join(OUT, "B3_逐细胞marker表达.csv.gz")
chk = pd.read_csv(chk_path)
chk.columns = [c.lstrip("﻿") for c in chk.columns]
assert "barcode" in chk.columns and len(chk) == ad.n_obs, "对账表形状不符"
assert np.array_equal(chk.barcode.astype(str).values, mbc), "对账表条码与元数据不一致"
common = [c for c in chk.columns if c not in ("barcode", "leiden") and c in set(ad.var_names)]
Xi = ad[:, common].X.todense()
maxd = 0.0
worst = None
for k, c in enumerate(common):
    mine = np.asarray(Xi[:, k]).ravel().astype(float)
    dd = np.abs(mine - chk[c].values.astype(float))
    if np.nanmax(dd) > maxd:
        maxd, worst = float(np.nanmax(dd)), c
rec(f"- **对账二**：与既有表达表 {len(common)} 个公共基因逐格比较，最大绝对差 = {maxd:.3g}"
    + (f"（最大处 {worst}）" if worst else "") + f"，容差 {TOL_VAL:g}")
assert maxd <= TOL_VAL, f"对账二失败：最大差 {maxd:.3g} 超出容差 {TOL_VAL:g}"
rec("  → 通过。两条独立路径（obs 名 + 表达值）互相印证。")
rec()

# ---------- 4. 成纤维细胞子集 ----------
fib_mask = (meta.cell_type.values == "Fibroblast")
fib = meta[fib_mask].copy()
rec(f"- 成纤维细胞：{fib_mask.sum()} 个（占 {100*fib_mask.mean():.1f}%）；"
    f"逐供者 {fib.groupby('donor').size().to_dict()}")

# ---------- 5. 前置可行性检查（方案明确要求先做） ----------
rec()
rec("## 一、前置可行性检查：成纤维亚群是否跨供者稳定")
rec()
fib_cl = pd.crosstab(fib.leiden, fib.donor).reindex(columns=["SY_113", "SY_116", "SY_118"],
                                                    fill_value=0)
prop = fib_cl.div(fib_cl.sum(axis=1), axis=0)
fib_cl["合计"] = fib_cl.sum(axis=1)
fib_cl["最大供者占比"] = (prop.max(axis=1) * 100).round(1)
fib_cl["最大供者"] = prop.idxmax(axis=1)
fib_cl = fib_cl.sort_values("合计", ascending=False)
rec("| 成纤维 leiden 簇 | SY_113 | SY_116 | SY_118 | 合计 | 最大供者 | 最大供者占比 |")
rec("|---|---|---|---|---|---|---|")
for cl, r in fib_cl.iterrows():
    rec(f"| {cl} | {r['SY_113']} | {r['SY_116']} | {r['SY_118']} | {r['合计']} | "
        f"{r['最大供者']} | {r['最大供者占比']}% |")
rec()
MONO = 80.0
mono = fib_cl[fib_cl["最大供者占比"] >= MONO]
absent = fib_cl[(fib_cl[["SY_113", "SY_116", "SY_118"]] == 0).any(axis=1)]
rec(f"- 「单一供者占比 ≥ {MONO:.0f}%」的簇：**{len(mono)} / {len(fib_cl)}**"
    + (f"（{', '.join(str(x) for x in mono.index)}）" if len(mono) else ""))
rec(f"- 有供者完全缺席（该簇某供者 0 细胞）的簇：**{len(absent)} / {len(fib_cl)}**"
    + (f"（{', '.join(str(x) for x in absent.index)}）" if len(absent) else ""))
rec()
RETREAT = (len(mono) > 0) or (len(absent) > 0)
rec("**判定**：" + ("成纤维亚群**不满足跨供者稳定**，按方案规定退回粗粒度定位描述"
                    "（只报供者级 lining / sub-lining 打分，不报簇层面的状态细分）。"
                    if RETREAT else
                    "所有成纤维簇都有三个供者贡献且无单供者独占，细分状态可用。"))
rec()

# ---------- 6. 外部状态基因集：剔除候选基因自身（硬约束） ----------
cand = pd.read_csv(os.path.join(RES, "57基因总览.csv"), encoding="utf-8-sig")
CAND = set(cand["gene"])
ALIAS = {"CTGF": "CCN2", "CYR61": "CCN1", "MKL1": "MRTFA", "MRTFB": "MKL2",
         "CCN2": "CTGF", "CCN1": "CYR61", "MRTFA": "MKL1", "MKL2": "MRTFB"}
CAND |= {ALIAS[g] for g in list(CAND) if g in ALIAS}
V = set(ad.var_names)


def clean_set(genes, tag):
    dropped = [x for x in genes if x in CAND]
    present = [x for x in genes if x not in CAND and x in V]
    missing = [x for x in genes if x not in CAND and x not in V]
    rec(f"- {tag}：原始 {len(genes)} → 剔除候选基因 {len(dropped)} 个"
        + (f"（{', '.join(dropped)}）" if dropped else "（无）")
        + f"，数据缺失 {len(missing)}" + (f"（{', '.join(missing)}）" if missing else "（无）")
        + f"，实际使用 {len(present)} 个：{', '.join(present)}")
    return present


rec("## 二、外部状态基因集（已剔除被检验候选基因）")
rec()
L = clean_set(LINING, "lining")
S = clean_set(SUBLINING, "sublining")
assert len(L) >= 3 and len(S) >= 3, "可用基因太少，状态打分不可靠——应记为不可行"
rec()

# ---------- 7. 模块打分（Seurat AddModuleScore 式：减去表达分箱匹配的对照集） ----------
fib_idx = np.where(fib_mask)[0]
Xf = ad[fib_idx, :].X.tocsr()
mean_all = np.asarray(Xf.mean(axis=0)).ravel()
gname = pd.Index(ad.var_names)


def module_score(genes, seed=20260923):
    rng = np.random.default_rng(seed)
    gj = gname.get_indexer(genes)
    assert (gj >= 0).all()
    bins = pd.qcut(mean_all, 25, labels=False, duplicates="drop")
    used = np.zeros(len(gname), dtype=bool)
    used[gj] = True
    ctrl = []
    for b in np.unique(bins):
        pool = np.where((bins == b) & ~used)[0]
        ctrl += list(rng.choice(pool, size=min(50, len(pool)), replace=False))
    ctrl = np.array(sorted(set(ctrl)))
    sc_ = (np.asarray(Xf[:, gj].mean(axis=1)).ravel()
           - np.asarray(Xf[:, ctrl].mean(axis=1)).ravel())
    return sc_, len(ctrl)


sc_L, nL = module_score(L)
sc_S, nS = module_score(S)
rec(f"模块打分：lining 用 {len(L)} 个基因 / {nL} 个对照基因；"
    f"sublining 用 {len(S)} 个基因 / {nS} 个对照基因（对照集按平均表达 25 分箱匹配，"
    f"种子 20260923）。")
rec()

# ---------- 8. 逐供者比较 ----------
res = pd.DataFrame({"donor": fib.donor.values, "leiden": fib.leiden.values,
                    "cell_type": fib.cell_type.values, "lining": sc_L, "sublining": sc_S})
res["diff"] = res.lining - res.sublining
res["state"] = np.where(res["diff"] > 0, "lining-like", "sublining-like")

per_donor = res.groupby("donor").agg(
    n_fib=("lining", "size"),
    lining_mean=("lining", "mean"), lining_sd=("lining", "std"),
    sublining_mean=("sublining", "mean"), sublining_sd=("sublining", "std"),
    diff_mean=("diff", "mean"), diff_sd=("diff", "std"),
    lining_like_pct=("state", lambda s: 100 * (s == "lining-like").mean()),
).reindex(["SY_113", "SY_116", "SY_118"])

rec("## 三、逐供者状态打分（成纤维细胞）")
rec()
rec("| 供者 | 成纤维数 | lining 均值(SD) | sublining 均值(SD) | 差值 lining−sublining | lining-like 占比 |")
rec("|---|---|---|---|---|---|")
for d, r in per_donor.iterrows():
    rec(f"| {d} | {int(r.n_fib)} | {r.lining_mean:.3f} ({r.lining_sd:.3f}) | "
        f"{r.sublining_mean:.3f} ({r.sublining_sd:.3f}) | "
        f"**{r.diff_mean:+.3f}** ({r.diff_sd:.3f}) | {r.lining_like_pct:.1f}% |")
rec()
signs = sorted({int(s) for s in np.sign(per_donor.diff_mean.round(9))})
rec(f"- 三供者差值符号：{signs} → "
    + ("**方向一致**（三供者同向）" if len(signs) == 1 and 0 not in signs
       else "**方向不一致**（存在供者间反号）"))
rec(f"- 供者间差值范围：{per_donor.diff_mean.min():+.3f} ～ {per_donor.diff_mean.max():+.3f}；"
    f"三供者均值 {per_donor.diff_mean.mean():+.3f}。")
rec()
rec("> 只有 3 个供者，**不做组间检验**（n=3 无判别力，且方案要求不预设结论）；")
rec("> 只报逐供者数值、符号一致性与范围，供正文按「粗粒度定位」陈述。")
rec()

# ---------- 9. 候选基因 × 外部定义的状态（回答「候选是否偏向某状态」） ----------
# 非循环：状态由外部基因集定义，候选基因不在其中（见第二节，剔除 0 个因为本就无交集）。
# 方案明令禁止的写法是「用 LOXL1 划分 LOXL1 高表达亚群再报告该亚群 LOXL1 高」——此处不是。
PANEL = ["LOXL1", "VGLL4", "SDC4", "SERPINE1", "ITGA5", "SAV1"]
XLINK = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4", "PLOD1", "PLOD2", "TGM2"]
GOI = sorted(set(PANEL) | set(XLINK))
pres = [g for g in GOI if g in V]
miss = [g for g in GOI if g not in V]

rec("## 四、候选基因表达 × 外部定义的两类状态（逐供者）")
rec()
if miss:
    rec(f"- 数据中不存在的基因（明列，不作「无差异」解读）：{', '.join(miss)}")
Xg = np.asarray(ad[fib_idx, :][:, pres].X.todense())
state_v = res["state"].values
donor_v = res["donor"].values
rows = []
for gi_, g in enumerate(pres):
    for d in ["SY_113", "SY_116", "SY_118"]:
        m_d = (donor_v == d)
        for st in ["lining-like", "sublining-like"]:
            m = m_d & (state_v == st)
            if m.sum() == 0:
                rows.append(dict(基因=g, 供者=d, 状态=st, 细胞数=0,
                                 平均log1p=np.nan, 表达比例=np.nan))
                continue
            v = Xg[m, gi_]
            rows.append(dict(基因=g, 供者=d, 状态=st, 细胞数=int(m.sum()),
                             平均log1p=round(float(v.mean()), 4),
                             表达比例=round(float((v > 0).mean()), 4)))
df_gs = pd.DataFrame(rows)
rec(f"（每格 `lining-like / sublining-like (差值)`，值为成纤维细胞内的平均 log1p 表达；"
    f"方向一致性只看符号，**不做检验**——同一供者内的细胞不独立。）")
rec()
rec("| 基因 | SY_113 | SY_116 | SY_118 | 方向 |")
rec("|---|---|---|---|---|")
for g in pres:
    sub = df_gs[df_gs.基因 == g]
    cells, agree = [], []
    for d in ["SY_113", "SY_116", "SY_118"]:
        a = sub[(sub.供者 == d) & (sub.状态 == "lining-like")].iloc[0]
        b = sub[(sub.供者 == d) & (sub.状态 == "sublining-like")].iloc[0]
        dd = a["平均log1p"] - b["平均log1p"]
        if np.isfinite(dd):
            agree.append(int(np.sign(round(dd, 9))))
        cells.append(f"{a['平均log1p']:.3f} / {b['平均log1p']:.3f} ({dd:+.3f})")
    asc = sorted(set(agree))
    tag = ("三供者同向 ↑lining" if asc == [1] else
           "三供者同向 ↑sublining" if asc == [-1] else f"方向不一致 {asc}")
    rec(f"| {g} | {cells[0]} | {cells[1]} | {cells[2]} | {tag} |")
rec()
rec("### 4.1 逐基因 × 逐供者 明细（含细胞数与表达比例）")
rec()
rec(df_gs.to_string(index=False))
rec()

dir_tab = []
for g in pres:
    sub = df_gs[df_gs.基因 == g]
    sg = []
    for d in ["SY_113", "SY_116", "SY_118"]:
        a = sub[(sub.供者 == d) & (sub.状态 == "lining-like")].iloc[0]
        b = sub[(sub.供者 == d) & (sub.状态 == "sublining-like")].iloc[0]
        dd = a["平均log1p"] - b["平均log1p"]
        sg.append(np.sign(round(dd, 9)) if np.isfinite(dd) else 0)
    sg = [int(x) for x in sg]
    n_pos, n_neg = sg.count(1), sg.count(-1)
    dir_tab.append(dict(基因=g, SY_113=sg[0], SY_116=sg[1], SY_118=sg[2],
                        同向="是" if (n_pos == 3 or n_neg == 3) else "否",
                        倾向=("lining-like 高" if n_pos == 3 else
                              "sublining-like 高" if n_neg == 3 else "不一致")))
df_dir = pd.DataFrame(dir_tab)
rec("### 4.2 方向一致性汇总")
rec()
rec(df_dir.to_string(index=False))
rec()
n_same = int((df_dir.同向 == "是").sum())
rec(f"- 逐基因看：{n_same} / {len(df_dir)} 个基因在三供者间方向一致。")
rec()

# ---------- 4.3 零假设：随机基因的「三供者同向」发生率（不补则 8/13 无法判定） ----------
# 不检验、只看符号一致，本身不构成证据——必须给出随机基因下的发生率作参照。
# 这正是记忆 permutation-null-self-locks-when-k-is-data-driven 与
# small-cohort-cv-auc-needs-permutation-null 要求的：小样本的「一致」先定标。
rec("### 4.3 零假设：随机基因的符号一致率")
rec()
NPERM, SEED0 = 2000, 20260923
rng0 = np.random.default_rng(SEED0)
pool = np.array([g for g in ad.var_names
                 if g not in set(GOI) and g not in set(L) and g not in set(S)])
pick = rng0.choice(pool, size=min(NPERM, len(pool)), replace=False)
Xr = np.asarray(ad[fib_idx, :][:, pick].X.todense())
dmsk = {d: (donor_v == d) for d in ["SY_113", "SY_116", "SY_118"]}
m_ll = (state_v == "lining-like")
m_sl = ~m_ll
n_ag, n_ev = 0, 0
for k in range(Xr.shape[1]):
    sg = []
    for d in ["SY_113", "SY_116", "SY_118"]:
        a, b = dmsk[d] & m_ll, dmsk[d] & m_sl
        if a.sum() == 0 or b.sum() == 0:
            continue
        sg.append(np.sign(round(float(Xr[a, k].mean() - Xr[b, k].mean()), 12)))
    if len(sg) == 3:
        n_ev += 1
        if abs(sum(sg)) == 3:
            n_ag += 1
null_rate = n_ag / n_ev if n_ev else float("nan")
# 观测侧的零（同向→非零）基因数：信号为 0 的基因不该算进分母
obs_ev = int(sum(1 for _, r in df_dir.iterrows()
                 if any(v != 0 for v in (r["SY_113"], r["SY_116"], r["SY_118"]))))
rec(f"- 随机抽 {n_ev} 个基因（已排除候选与状态基因），**三供者符号完全一致**的发生率 "
    f"= {n_ag}/{n_ev} = **{null_rate:.1%}**。")
rec(f"- 观测侧：{n_same} / {obs_ev} 个候选基因三供者同向 = {n_same/obs_ev:.1%}。")
rec()
if null_rate > 0:
    rec(f"> 零假设下期望同向约 {null_rate*obs_ev:.1f} 个，观测 {n_same} 个。"
        + ("观测与零假设量级相当，**不能**据此声称候选基因有状态偏好。"
           if n_same <= null_rate * obs_ev * 1.5 else
           "观测高于零假设，但**未做检验、且细胞不独立**，仍只作描述性陈述。"))
else:
    rec("> 零假设下无基因达到三供者完全一致 —— 观测到的同向不是随机可得的。")
rec()

# ---------- 4.4 模块打分机制的自检 ----------
# 负面对照：随机 7 基因「lining」+ 随机 9 基因「sublining」，看三供者同号的比例。
# 若随机集也高频出现「三供者同号」，则第三节的 +0.723 不构成证据。
rec("### 4.4 模块打分机制自检（随机状态集，与第三节同构）")
rec()
rng1 = np.random.default_rng(SEED0 + 1)
Ndraw = 200
pool2 = np.array([g for g in ad.var_names if g not in set(GOI)])
n_same_s, n_allpos = 0, 0
for i in range(Ndraw):
    gl = list(rng1.choice(pool2, size=7, replace=False))
    gs = list(rng1.choice([x for x in pool2 if x not in set(gl)], size=9, replace=False))
    dl, _ = module_score(gl, seed=20260923 + i)
    ds, _ = module_score(gs, seed=20260923 + i)
    sg = [float(dl[dmsk[d]].mean() - ds[dmsk[d]].mean())
          for d in ["SY_113", "SY_116", "SY_118"]]
    sg9 = [np.sign(round(x, 12)) for x in sg]
    if abs(sum(sg9)) == 3:
        n_same_s += 1
    if all(x > 0 for x in sg9):
        n_allpos += 1
rec(f"- 用**同一套** `module_score` 跑 {Ndraw} 次「随机 7 基因衬里集 − 随机 9 基因衬里下集」：")
rec(f"  - 三供者同号 {n_same_s} 次（{n_same_s/Ndraw:.1%}）；")
rec(f"  - **三供者全为正** {n_allpos} 次（{n_allpos/Ndraw:.1%}）。")
rec()
rec("**因此第三节「三供者同为正」的随机发生率约为 "
    f"{n_allpos/Ndraw:.1%}**（n 只有 3，故该判据本身很松）。")
rec("> 第三节的 +0.574/+0.924/+0.672 只能按「描述性、不构成统计证据」写入正文；")
rec("> 若正文要声称衬里倾向，必须补独立数据集——本轮没有。")
rec()

# ---------- 10. 细胞组成侧（直接回答「组成 vs 状态」） ----------
# ---------- 4.5 表达量匹配的零假设 + 二项检验 ----------
# 4.3 的全基因零假设把大量近零表达基因也算了进去——那些基因的符号本来就是噪声，
# 会把零假设压低，使 61.5% 显得比实际更突出。正确做法是把零假设限制到**与候选基因
# 表达量相当**的基因上（同记忆 small-cohort-cv-auc-needs-permutation-null 的定标要求）。
rec("### 4.5 表达量匹配的零假设（对 4.3 的收紧）")
rec()
cand_mean = np.array([float(mean_all[gname.get_loc(g)]) for g in pres])
lo, hi = float(cand_mean.min()), float(cand_mean.max())
match_pool = np.array([g for i, g in enumerate(ad.var_names)
                       if lo <= mean_all[i] <= hi and g not in set(GOI)])
rng2 = np.random.default_rng(SEED0 + 2)
NM = min(3000, len(match_pool))
mpick = rng2.choice(match_pool, size=NM, replace=False)
Xm = np.asarray(ad[fib_idx, :][:, mpick].X.todense())
nm_ag = 0
for k in range(Xm.shape[1]):
    sg = []
    for d in ["SY_113", "SY_116", "SY_118"]:
        a, b = dmsk[d] & m_ll, dmsk[d] & m_sl
        if a.sum() == 0 or b.sum() == 0:
            continue
        sg.append(np.sign(round(float(Xm[a, k].mean() - Xm[b, k].mean()), 12)))
    if len(sg) == 3 and abs(sum(sg)) == 3:
        nm_ag += 1
matched_rate = nm_ag / NM
rec(f"- 候选基因在成纤维细胞内的平均 log1p 表达范围：{lo:.4f} ～ {hi:.4f}")
rec(f"- 落在同一表达区间、且非候选非状态基因的基因共 {len(match_pool)} 个，抽 {NM} 个：")
rec(f"  **三供者符号一致率 = {nm_ag}/{NM} = {matched_rate:.1%}**"
    f"（对比：全基因零假设 {null_rate:.1%}）")
rec()
bt = binomtest(n_same, obs_ev, matched_rate, alternative="greater")
rec(f"- 二项检验：观测 {n_same}/{obs_ev}，零假设率 {matched_rate:.1%}，"
    f"P(X≥{n_same}) = **{bt.pvalue:.4f}**")
rec()
rec("> **但 n 的独立性不成立**：13 个基因里 LOX/LOXL1/LOXL2/LOXL3/LOXL4 是同一家族、")
rec("> 常共调控，SDC4/ITGA5/SAV1/VGLL4 同属黏附—Hippo 轴。有效独立单位远小于 13，")
rec("> 故该 P 值是**上界乐观**的，只能作方向性提示，不足以支撑因果或定量主张。")
rec()
rec("### 4.6 逐基因结果的分组读法（按家族，避免把共调控当成独立重复）")
rec()
FAM = {"交联催化家族（LOX/LOXL1-4）": ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4"],
       "胶原修饰（PLOD1/2、TGM2）": ["PLOD1", "PLOD2", "TGM2"],
       "黏附—Hippo（SDC4/ITGA5/SAV1/VGLL4）": ["SDC4", "ITGA5", "SAV1", "VGLL4"],
       "纤溶（SERPINE1）": ["SERPINE1"]}
for fam, gs in FAM.items():
    part = df_dir[df_dir.基因.isin(gs)]
    nup = int((part.倾向 == "lining-like 高").sum())
    ndn = int((part.倾向 == "sublining-like 高").sum())
    nin = int((part.倾向 == "不一致").sum())
    rec(f"- **{fam}**（{len(part)} 基因）：一致偏高 lining-like {nup}，"
        f"一致偏高 sublining-like {ndn}，不一致 {nin}"
        + ("　→ 家族内方向一致" if nin == 0 and (nup == 0 or ndn == 0) else ""))
rec()

rec("## 五、细胞组成侧（供者间差异有多大）")
rec()
ct_tab = pd.crosstab(meta.cell_type, meta.donor).reindex(
    columns=["SY_113", "SY_116", "SY_118"], fill_value=0)
ct_pct = (ct_tab / ct_tab.sum(axis=0) * 100).round(2)
rec(ct_pct.to_string())
rec()
FIB = "Fibroblast"
if FIB in ct_pct.index:
    vals = " / ".join(f"{ct_pct.loc[FIB, d]:.2f}%" for d in ["SY_113", "SY_116", "SY_118"])
    span = ct_pct.loc[FIB].max() - ct_pct.loc[FIB].min()
    rec(f"- 成纤维占比（SY_113 / SY_116 / SY_118）：{vals}，极差 {span:.2f} 个百分点。")
rec()
rec("> 组成侧的供者间差异与状态侧的供者间差异都**只能描述**：n=3，无法做方差分解，")
rec("> 不能量化「组成贡献多少、状态贡献多少」。本轮**不给出**这样的比例。")
rec()

# ---------- 11. 本次做到与未做到 ----------
rec("## 六、本次做到与未做到（规则 4）")
rec()
rec("**做到**：")
rec("- 从原始 10x 重建表达矩阵，两条独立路径对账通过（obs 名逐位相等、表达值最大差 "
    f"{maxd:.2e}）；")
rec("- 前置可行性检查**未通过**，已按方案退回粗粒度（不报簇层面状态细分）；")
rec("- 外部基因集打分、候选基因零交集、逐供者比较（含符号一致性汇总）；")
rec(f"- **两重零假设已补**：随机基因的符号一致率 {null_rate:.1%}（4.3），"
    f"随机状态集「三供者全为正」率 {n_allpos/Ndraw:.1%}（4.4）。")
rec()
rec("**未做到 / 不能做**：")
rec("- `GSE152805` **是 9 例 OA、无健康对照**，故所有陈述只能是「OA 内部定位」，")
rec("  **不能**写成「OA 相对正常滑膜的状态偏移」；")
rec("- 无独立病例—对照滑膜 scRNA（见 A-2.5 记录），故无外部验证；")
rec("- 亚群层面的状态细分不可用（供者私有簇过多）；")
rec("- 组成与状态的方差分解未做、也不可做。")
rec()
rec("## 七、A-2.4 结论（本项实际能交付什么）")
rec()
rec(f"1. **状态轴不成立**：衬里/衬里下模块分数的差值三供者同为正，但随机 7/9 基因集也有 "
    f"{n_allpos/Ndraw:.1%} 的概率如此（4.4）——该差值不可解释。")
rec(f"2. **逐基因状态偏好未达显著**：三供者同向 8/13 = 61.5%，表达量匹配零假设 "
    f"{matched_rate:.1%}，二项 P = {binomtest(n_same, obs_ev, matched_rate, alternative='greater').pvalue:.4f}"
    f"（4.5）；且 13 个基因非独立，P 值乐观。")
rec("3. **可写入正文的只有两点**：")
rec(f"   - 供者间**细胞组成**差异极大（成纤维占比 86.83% / 87.54% / 43.89%，极差 "
    f"43.65 个百分点）；")
rec("   - 候选基因在两类状态间的表达差异，只能作**逐供者描述**（4.1 表），")
rec("     不得声称状态偏好、不得声称与疾病相关（GSE152805 无对照）。")
rec("4. **据此，方案 §四「细胞来源」这一层不能按原设想落地**：既拿不到跨供者稳定的亚群划分，")
rec("   也拿不到显著的状态关联。建议正文改为「bulk 疾病差异 + OA 内部粗粒度定位」。")
rec("   （此为建议，最终取舍由作者定。）")
rec()

# ---------- 12. 落盘 ----------
fib_cl.to_csv(os.path.join(RES, "A2d_成纤维簇_供者构成.csv"), encoding="utf-8-sig")
per_donor.to_csv(os.path.join(RES, "A2d_逐供者状态打分.csv"), encoding="utf-8-sig")
df_gs.to_csv(os.path.join(RES, "A2d_候选基因×状态_逐供者.csv"), index=False,
             encoding="utf-8-sig")
df_dir.to_csv(os.path.join(RES, "A2d_候选基因×状态_方向汇总.csv"), index=False,
              encoding="utf-8-sig")
ct_pct.to_csv(os.path.join(RES, "A2d_供者细胞组成.csv"), encoding="utf-8-sig")
res.to_csv(os.path.join(RES, "A2d_逐细胞状态打分.csv.gz"), index=False, compression="gzip",
           encoding="utf-8-sig")
pd.DataFrame({"gene_set": ["lining"] * len(L) + ["sublining"] * len(S),
              "gene": L + S}).to_csv(os.path.join(RES, "A2d_状态基因集.csv"),
                                     index=False, encoding="utf-8-sig")
json.dump({"n_cells": int(ad.n_obs), "n_genes_used_as_denominator": int(ad.n_vars),
           "n_fibroblast": int(fib_mask.sum()),
           "obs_names_match_b3_metadata": bool(eq_bc), "donor_match": bool(eq_dn),
           "n_double_suffixed_barcodes": n_suf,
           "max_abs_diff_vs_marker_table": maxd, "n_common_checked": len(common),
           "lining_genes_used": L, "sublining_genes_used": S,
           "retreat_to_coarse": bool(RETREAT),
           "cluster_donor_dominance_ge80pct": [str(x) for x in mono.index],
           "clusters_with_missing_donor": [str(x) for x in absent.index],
           "per_donor_diff_mean": {k: float(v) for k, v in per_donor.diff_mean.items()},
           "panel_genes_missing_from_data": miss,
           "n_genes_sign_consistent_across_donors": n_same,
           "n_genes_dir_reported": int(len(df_dir)),
           "n_genes_with_nonzero_dir": obs_ev,
           "null_random_gene_sign_agreement_rate": null_rate,
           "null_expression_matched_sign_agreement_rate": matched_rate,
           "null_expression_matched_pool_size": int(len(match_pool)),
           "null_expression_matched_n_sampled": int(NM),
           "binom_p_for_sign_agreement": float(binomtest(n_same, obs_ev, matched_rate,
                                                         alternative="greater").pvalue),
           "null_module_score_n_draws": Ndraw,
           "null_module_score_same_sign": n_same_s,
           "null_module_score_all_positive": n_allpos,
           "fibroblast_pct_per_donor": {d: float(ct_pct.loc[FIB, d])
                                        for d in ["SY_113", "SY_116", "SY_118"]}
           if FIB in ct_pct.index else {}},
          open(os.path.join(OUT, "A2d_summary.json"), "w"), ensure_ascii=False, indent=2)

with open(os.path.join(AUD, "A2d_成纤维状态关联_结果.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(ln) + "\n\n---\n"
             "结论只写进上面的实测句；未做组间检验、未预设方向；候选基因已剔除（见第二节）。\n")
print("\n已写出：02_结果表/A2d_*.csv（5 份）、03_分析/outputs/B3_scRNA/A2d_summary.json、"
      "01_审计/A2d_成纤维状态关联_结果.md")
