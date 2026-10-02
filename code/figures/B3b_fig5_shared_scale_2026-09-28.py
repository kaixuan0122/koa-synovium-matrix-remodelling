#!/usr/bin/env python3
# B3b 图 5 重绘：四面板共享色阶 + 首面板改池化细胞（2026-09-28）
#
# 为什么另起一个脚本，而不改 B3 里的绘图段：
#   B3 是从三份原始 10x 矩阵重建单细胞结果的**分析**脚本（约 10 分钟、依赖 ~/bio-data 下的
#   raw 数据）。图 5 的每一个数都已经冻结在交付包的源数据表里，重画图不该再依赖原始矩阵。
#   本脚本只读冻结源表，秒级可跑，交付包的读者可以据此复现图 5 —— 这也是 绘图代码_EN/ 里
#   其余四个图脚本的形态。
#
# 回应 Reviewer 3 的 R3-M1（Figure 5 re-export）三条：
#   ① 原图四个面板各自独立取色阶，却共用最后那个面板生成的色条 —— 色条上的刻度
#      描述的是 SY_118，读者会以为四个面板同尺。现改为 vmin/vmax 全局面板共享。
#   ② 原图 "All donors" 面板是**三供者比例等权平均**（SERPINE1 内皮 = 74.60%），
#      而正文与图注引的是**合并细胞**值（68.31%）。现首面板改读池化细胞表。
#   ③ 原图 suptitle 写 "dot size = % expressing" 但图上没有面积图例。现补图例。
#
# 逐值断言：图注（主稿 Figure Legends 区 Figure 5 段）引到的每一个数都在此复核，
# 断言失败即非零退出 —— 图注是待画的规格书，不是可以目测放过的散文。
#
# 输入：源数据_EN/scRNA_panel_expr_by_celltype.csv        （池化细胞，10,600 个）
#       源数据_EN/scRNA_panel_expr_by_donor_celltype.csv （逐供者 × 细胞类型）
# 输出：Figure5_scRNA_dotplot.pdf / .png
#       源数据_EN/Fig5_values.csv（逐点数值 + 口径说明，供审稿人逐格回核）
# 日志：调用方重定向到 03_分析/logs/B3b_图5重绘.log

import os, sys, csv, hashlib
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def find_root():
    """在若干候选根目录中找含 源数据_EN/ 的那一个。

    交付包被整体拷走后，脚本与源数据表的相对位置不变；本机工作树里则是 03_图表与源数据/。
    两处都试，找不到就明确报错，不猜、不用空表兜底。
    """
    cands = [
        os.path.join(HERE, "..", "..", ".."),                    # 工作树 scripts/ → 树根
        os.path.join(HERE, ".."),                                # 绘图代码_EN/ → 03_图表与源数据/
        os.path.join(HERE, "..", ".."),                          # scripts/ → 01_冻结结果与代码/
    ]
    for c in cands:
        for sub in ("源数据_EN", "03_图表与源数据/源数据_EN"):
            p = os.path.normpath(os.path.join(c, sub))
            if os.path.isfile(os.path.join(p, "scRNA_panel_expr_by_celltype.csv")):
                return os.path.normpath(c)
    sys.exit(f"✗ 找不到 源数据_EN/scRNA_panel_expr_by_celltype.csv；已试：\n  " +
             "\n  ".join(os.path.normpath(c) for c in cands))

ROOT = find_root()
# 两种布局：工作树（ROOT 下有 01_冻结结果与代码/，E3 从 04_图表/ 取图）与交付包
# （只有 绘图代码_EN/ + 源数据_EN/）。交付包里**一律写到同级的 重跑产出/**：交付的图件与
# 源数据是待核对的成品，读者跑一遍脚本不该把它们就地改掉 —— 旧版 delivered_fig() 在交付包
# 布局下正好落回 ROOT 本身，等于把交付的 Figure5_scRNA_dotplot.pdf 覆盖重写。
PACKAGE = not os.path.isdir(os.path.join(ROOT, "01_冻结结果与代码"))
OUTDIR = os.path.join(HERE, "重跑产出") if PACKAGE else ROOT
if PACKAGE:
    os.makedirs(os.path.join(OUTDIR, "源数据_EN"), exist_ok=True)


def src(name):
    for sub in ("源数据_EN", "03_图表与源数据/源数据_EN"):
        p = os.path.join(ROOT, sub, name)
        if os.path.isfile(p):
            return p
    sys.exit(f"✗ 源数据缺失：{name}（根目录 {ROOT}）")

# 落点两处，同一份字节：
#   ① 冻结件名  01_冻结结果与代码/04_图表/Fig_scRNA_dotplot.{pdf,png}
#      —— E3_图件交付组装 就是从这里取的；只写交付名会让 E3 重跑时把新图换回旧图。
#   ② 交付件名  03_图表与源数据/Figure5_scRNA_dotplot.{pdf,png}
# 以上是**工作树**里的落法。交付包布局（没有 01_冻结结果与代码/）下两处都改成
# 绘图代码_EN/重跑产出/ —— 交付件是待核对的成品，不能被读者的一次运行就地覆盖。
# 先写 ① 再 copy 到 ②（不 savefig 两次：matplotlib 每次落盘都带新的 /CreationDate 与 /ID，
# 两次 savefig 的字节必然不同，那样「两处一致」就成了一句空话）。
def frozen_fig(name):
    """冻结件落点。交付包里没有 01_冻结结果与代码/（E3 已把图组装成交付名），
    此时返回 None，只写交付件名 —— 不在别人的包里凭空造一个 01_… 目录出来。"""
    d = os.path.join(ROOT, "01_冻结结果与代码", "04_图表")
    if not os.path.isdir(os.path.join(ROOT, "01_冻结结果与代码")):
        return None
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)

def delivered_fig(name):
    if PACKAGE:
        return os.path.join(OUTDIR, name)          # 交付包：只写 重跑产出/，不碰交付件
    d = os.path.join(ROOT, "03_图表与源数据")
    if not os.path.isdir(d):
        d = ROOT
    return os.path.join(d, name)

fails = []
def must(ok, msg, detail=None):
    tag = "OK  " if ok else "FAIL"
    line = f"  {tag} {msg}"
    if detail is not None:
        line += f"   {detail}"
    print(line, flush=True)
    if not ok:
        fails.append(msg)

print("=" * 78)
print("B3b 图 5 重绘：四面板共享色阶 + 首面板改池化细胞")
print(f"根目录：{ROOT}")
print("=" * 78)

CT = pd.read_csv(src("scRNA_panel_expr_by_celltype.csv"))
DC = pd.read_csv(src("scRNA_panel_expr_by_donor_celltype.csv"))

GENES = ["LOXL1", "SERPINE1", "SDC4", "VGLL4", "ITGA5", "SAV1"]      # 行，主稿面板六gene
CTS   = ["Fibroblast", "Macrophage", "T cell", "Endothelial", "Mast"]  # 列，沿用原图列序
DONORS = ["SY_113", "SY_116", "SY_118"]

# ── 0. 输入完整性：先把「表本身是不是完整的」钉死，再去画 ──────────────────────
must(set(CT["cell_type"]) == set(CTS), "池化表恰含五个细胞群体", sorted(set(CT["cell_type"])))
must(set(CT["gene"]) == set(GENES), "池化表恰含六个面板gene", sorted(set(CT["gene"])))
must(len(CT) == 30, "池化表 5 × 6 = 30 行", len(CT))
must(set(DC["donor"]) == set(DONORS), "逐供者表恰含三个供者", sorted(set(DC["donor"])))
must(len(DC) == 90, "逐供者表 3 × 5 × 6 = 90 行", len(DC))

n_pool = CT.groupby("cell_type")["n_cells"].first()
must(int(n_pool.sum()) == 10600, "池化细胞总数 10,600", int(n_pool.sum()))
_mast = int(n_pool["Mast"])
must(_mast == 130, "图注所引「肥大细胞 130 个」", _mast)

# 逐供者n_cells：三表合起来必须等于池化总数（两个文件是同一批细胞的两种切法）
n_dc = DC[DC["gene"] == GENES[0]].set_index(["donor", "cell_type"])["n_cells"].unstack()
must(n_dc.sum().sum() == 10600, "逐供者表三供者相加 = 10,600", int(n_dc.sum().sum()))
for ct in CTS:
    must(int(n_dc[ct].sum()) == int(n_pool[ct]),
         f"{ct} 逐供者n_cells之和 = 池化表值", f"{int(n_dc[ct].sum())} vs {int(n_pool[ct])}")

# ── 1. 图注所引数值逐条复核 ───────────────────────────────────────────────────
print("\n-- 1. 图注引用的数值 --")
def cell(df, ct, gene, col):
    return float(df[(df.cell_type == ct) & (df.gene == gene)][col].iloc[0])

# 成纤维比例（分母 = 该供者全部细胞，含 Mast）
fib_share = {d: n_dc.loc[d, "Fibroblast"] / n_dc.loc[d].sum() * 100 for d in DONORS}
must([round(v, 2) for v in fib_share.values()] == [86.83, 87.54, 43.89],
     "三供者成纤维比例 86.83 / 87.54 / 43.89", {k: round(v, 2) for k, v in fib_share.items()})
_spread_r = round(fib_share["SY_116"] - fib_share["SY_118"], 2)
_spread_u = round(round(fib_share["SY_116"], 2) - round(fib_share["SY_118"], 2), 2)
print(f"     极差：先舍入再相减 = {_spread_u} 个百分点；用全精度相减 = {_spread_r} 个百分点")
print(f"     （图注写 43.65 = 报告值 87.54 − 43.89。差 0.01 纯属舍入次序，两读皆可；")
print(f"      本脚本登记该差别，不改写正文——见 07_核查与修改记录。）")
must(_spread_u == 43.65, "图注的 43.65 是报告值之差", _spread_u)

must(abs(cell(CT, "Fibroblast", "LOXL1", "expr_fraction") - 12.65) < 1e-9,
     "LOXL1 成纤维 12.65%", cell(CT, "Fibroblast", "LOXL1", "expr_fraction"))
_o = {ct: cell(CT, ct, "LOXL1", "expr_fraction") for ct in CTS if ct != "Fibroblast"}
must(abs(min(_o.values()) - 1.64) < 1e-9 and abs(max(_o.values()) - 6.08) < 1e-9,
     "LOXL1「其余四类 1.6%–6.1%」", {k: round(v, 2) for k, v in _o.items()})
must(abs(cell(CT, "Fibroblast", "SDC4", "expr_fraction") - 37.37) < 1e-9,
     "SDC4 成纤维 37.37%", cell(CT, "Fibroblast", "SDC4", "expr_fraction"))
must(abs(cell(CT, "Endothelial", "VGLL4", "expr_fraction") - 64.21) < 1e-9,
     "VGLL4 内皮 64.21%", cell(CT, "Endothelial", "VGLL4", "expr_fraction"))
must(abs(cell(CT, "Endothelial", "SERPINE1", "expr_fraction") - 68.31) < 1e-9,
     "SERPINE1 内皮 68.31%（池化细胞口径）", cell(CT, "Endothelial", "SERPINE1", "expr_fraction"))
_s = {ct: cell(CT, ct, "SERPINE1", "expr_fraction") for ct in CTS if ct != "Endothelial"}
must(abs(min(_s.values()) - 14.62) < 1e-9 and abs(max(_s.values()) - 39.70) < 1e-9,
     "SERPINE1「其余四类 14.62%–39.70%」", {k: round(v, 2) for k, v in _s.items()})

# R3-M1 的两个口径必须真的不同，否则「改口径」这一步是空转
_eqw = float(DC[(DC.gene == "SERPINE1") & (DC.cell_type == "Endothelial")]["expr_fraction"].mean())
must(abs(_eqw - 74.60) < 0.01, "等权供者均值确为 74.60%（原首面板画的就是它）", round(_eqw, 2))
must(abs(_eqw - 68.31) > 5, "池化与等权两个口径确实不同（改口径不是空转）",
     f"{round(_eqw,2)} vs 68.31")

# ── 2. 摆成四个面板 ──────────────────────────────────────────────────────────
# 面板 0 = 池化细胞（首面板）；其余三个按供者号。面板内 grid 一律 6 行 × 5 列。
panels = [("Pooled cells (all donors)", CT)]
for d in DONORS:
    panels.append((d, DC[DC["donor"] == d]))

GRID = {}
for title, df in panels:
    p = df.pivot_table(index="gene", columns="cell_type", values="expr_fraction").reindex(index=GENES, columns=CTS)
    m = df.pivot_table(index="gene", columns="cell_type", values="mean_log1p").reindex(index=GENES, columns=CTS)
    n = df.pivot_table(index="gene", columns="cell_type", values="n_cells").reindex(index=GENES, columns=CTS)
    GRID[title] = (p, m, n)
    must(p.notna().all().all(), f"{title} 面板 6×5 网格无缺格", int(p.isna().sum().sum()))
for title, (p, m, n) in GRID.items():
    row0 = n.iloc[0].values
    must(all((n.iloc[i].values == row0).all() for i in range(6)),
         f"{title} 面板内n_cells逐gene一致（n 是列属性不是单元格属性）")

# 共享色阶：vmin/vmax 取自四个面板的全部 120 个格
allv = np.concatenate([m.values.ravel() for _, (p, m, n) in GRID.items()])
VMIN, VMAX = float(np.nanmin(allv)), float(np.nanmax(allv))
must(VMIN == 0.0, "共享色阶下界为 0", VMIN)
must(abs(VMAX - 1.8664) < 1e-4, "共享色阶上界 = SY_113 内皮 SERPINE1 的 mean log1p（图注值 1.8664）", VMAX)
_113 = float(GRID["SY_113"][1].loc["SERPINE1", "Endothelial"])
must(abs(_113 - 1.8664) < 1e-4, "SY_113 内皮 SERPINE1 mean_log1p = 1.8664", _113)
for lbl, exp in (("SY_116", 1.5776), ("SY_118", 0.7625)):
    got = float(GRID[lbl][1].loc["SERPINE1", "Endothelial"])
    must(abs(got - exp) < 1e-4, f"{lbl} 内皮 SERPINE1 mean_log1p = {exp}", got)

# ── 3. 出图 ─────────────────────────────────────────────────────────────────
print("\n-- 3. 出图 --")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

S_MIN, S_MAX = 8.0, 300.0
def area(frac):
    """面积按expr_fraction线性映射。S_MIN 保证 0% 的格子仍可见（0% 不是「无数据」）。"""
    return S_MIN + (S_MAX - S_MIN) * (frac / 100.0)

cmap = plt.get_cmap("viridis")
norm = Normalize(vmin=VMIN, vmax=VMAX)

# 版式：[面板 1..4][面积图例条带]，色条再在最右。
# 为什么单给图例一条带：四个面板的右下角都被 Mast 列的格占着（Mast 是最右列），
# 图例放进任何一个面板都会盖住真实数据点。第一版就是那样，SY_118 面板的
# Mast × VGLL4/ITGA5/SAV1 三格被遮——盖住数据的图例本身就是一处缺陷。
fig, axes = plt.subplots(1, len(panels) + 1, figsize=(3.4 * len(panels) + 1.3, 5.0),
                         sharey=True, gridspec_kw={"width_ratios": [3.4] * len(panels) + [1.3]})
ax_leg = axes[-1]
axes = axes[:-1]
for ax, (title, _) in zip(axes, panels):
    p, m, n = GRID[title]
    xs, ys, ss, cs = [], [], [], []
    for i, g in enumerate(GENES):
        for j, c in enumerate(CTS):
            xs.append(j); ys.append(i)
            ss.append(area(float(p.loc[g, c])))
            cs.append(float(m.loc[g, c]))
    ax.scatter(xs, ys, s=ss, c=cs, cmap=cmap, norm=norm, edgecolors="k", linewidths=0.4)
    nrow = n.iloc[0].values.astype(int)
    ax.set_xticks(range(len(CTS)))
    ax.set_xticklabels([f"{c}\nn = {v:,}" for c, v in zip(CTS, nrow)], rotation=45, ha="right", fontsize=7.5)
    ax.set_yticks(range(len(GENES))); ax.set_yticklabels(GENES, fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.set_xlim(-0.6, len(CTS) - 0.4); ax.set_ylim(-0.6, len(GENES) - 0.4)
    ax.grid(alpha=0.25, linewidth=0.4)
    ax.invert_yaxis()

# 面积图例：Line2D 的 markersize 是直径(pt)，scatter 的 s 是面积(pt²)，故 ms = sqrt(s)
_h = [Line2D([], [], linestyle="none", marker="o", markerfacecolor="0.75",
             markeredgecolor="k", markeredgewidth=0.4, markersize=np.sqrt(area(v)), label=f"{v}%")
      for v in (5, 25, 50, 90)]
ax_leg.axis("off")
ax_leg.legend(handles=_h, title="% of cells expressing", loc="center",
              fontsize=7.5, title_fontsize=8, frameon=False, labelspacing=2.2,
              handletextpad=1.4, borderpad=0.2)

# 色条由全局面板共享的归一化造出，不再借某一个面板的 sct
cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=list(axes) + [ax_leg],
                  fraction=0.02, pad=0.015)
cb.set_label("Mean log1p(CP10K)  (shared across all four panels)", fontsize=8)

fig.suptitle("Panel gene expression across cell types", fontsize=11)

_bad = [t.get_text() for t in fig.findobj()
        if hasattr(t, "get_text") and any(ord(ch) > 127 for ch in t.get_text())]
must(not _bad, "图内文本全为 ASCII（本机无中文字形，非 ASCII 会静默变豆腐块）", _bad)

# 图例不得盖住任何数据点：在显示坐标系里逐点算外接方框，与图例外框求交。
# 「看图没问题」不算结案——第一版就是靠肉眼漏掉的（图例压住 SY_118 的 Mast 列三格）。
fig.canvas.draw()
_LB = ax_leg.get_legend().get_window_extent()
_occ = []
for ax, (title, _) in zip(axes, panels):
    p_, m_, n_ = GRID[title]
    for i, g in enumerate(GENES):
        for j, c in enumerate(CTS):
            xd, yd = ax.transData.transform((j, i))
            r = np.sqrt(area(float(p_.loc[g, c]))) / 2.0 * fig.dpi / 72.0   # pt → px
            if (_LB.x0 < xd + r and xd - r < _LB.x1 and
                    _LB.y0 < yd + r and yd - r < _LB.y1):
                _occ.append(f"{title}/{g}×{c}")
must(not _occ, "面积图例未覆盖任何数据点", _occ)

import shutil
DELIV_PDF, DELIV_PNG = (delivered_fig("Figure5_scRNA_dotplot.pdf"),
                        delivered_fig("Figure5_scRNA_dotplot.png"))
FRZ_PDF, FRZ_PNG = frozen_fig("Fig_scRNA_dotplot.pdf"), frozen_fig("Fig_scRNA_dotplot.png")
# 没冻结树（交付包里跑）时，冻结件名退回交付件名，下面的 copy 自动变成空操作
FRZ_PDF, FRZ_PNG = FRZ_PDF or DELIV_PDF, FRZ_PNG or DELIV_PNG
# /CreationDate 与 /ModDate 每次落盘都变，会让「连跑两遍 sha256 相同」这条固定点判据
# 在 PDF 上恒假（PNG 本来就没有时间戳）。PDF 是交付件，去掉时间戳才能做字节级固定点。
_META = {"CreationDate": None, "ModDate": None}
fig.savefig(FRZ_PDF, bbox_inches="tight", metadata=_META)
fig.savefig(FRZ_PNG, dpi=200, bbox_inches="tight")
plt.close(fig)
FIG_PDF, FIG_PNG = DELIV_PDF, DELIV_PNG
for a, b in ((FRZ_PDF, DELIV_PDF), (FRZ_PNG, DELIV_PNG)):
    if os.path.abspath(a) != os.path.abspath(b):
        shutil.copy2(a, b)
        must(sha256(a) == sha256(b),
             f"{os.path.basename(b)} 与冻结件 {os.path.basename(a)} 逐字节一致")
for p_ in (FIG_PDF, FIG_PNG):
    sz = os.path.getsize(p_)
    must(sz > 5000, f"{os.path.basename(p_)} 已写出且非空", f"{sz} 字节")
with open(FIG_PDF, "rb") as fh:
    must(fh.read(5) == b"%PDF-", "PDF 头正确")

# ── 4. PDF 文本层回读：图上真的写了这些字，而不是脚本自称写了 ─────────────────
# 图注必须与图对得上；只查「文件存在」是 E2 未运行却照常通过的老机制。
#
# 用仓内既有的 pdftext.pdf_text（D4b/E3 同款）。**不要**自己解流后按字面量搜：
# 流里的括号是转义的（\(CP10K\)），且 matplotlib 会按字距把一个标签切成多个 Tj 操作数，
# 自行拼字符串必然假失败——本脚本第一版就踩了这个坑（5 个标签假 FAIL）。
# 依赖定位：本脚本要能在两个布局里跑 —— 工作树的 03_分析/scripts/，和三份源表一起被拷走的
# 绘图代码_EN/。为此把 pdftext.py 逐字节replication到绘图代码_EN/ 随件走，并在此按目录表查找；
# 查不到就**判失败**而不是跳过：文本层核对正是 E2 未运行却照常通过的机制（计划 C5），
# 静默跳过等于把这道闸门关掉。失败后仍继续到第 5 节写逐点表，产物与失败同时留下。
def _load_pdftext():
    cands = [
        HERE,
        os.path.join(ROOT, "01_冻结结果与代码", "03_分析", "scripts"),
        os.path.join(HERE, "..", "..", "01_冻结结果与代码", "03_分析", "scripts"),
    ]
    for d in cands:
        if os.path.isfile(os.path.join(d, "pdftext.py")):
            sys.path.insert(0, os.path.normpath(d))
            import importlib
            return importlib.import_module("pdftext").pdf_text
    must(False, "找得到 pdftext.py（PDF 文本层核对依赖它；不得静默跳过）",
         "；已试 " + "、".join(os.path.normpath(c) for c in cands))
    return None

_pdf_text = _load_pdftext()
if _pdf_text is None:
    print("  —— 第 4 节 PDF 文本层核对**未执行**，见上方 FAIL（按纪律登记为未完成，不是通过）")
_recs, _raw = _pdf_text(FIG_PDF) if _pdf_text else ([], "")
TL_SP = " ".join(t for _, _, t in _recs)     # 字距切分过的标签能靠空格拼回
TL_NS = "".join(t for _, _, t in _recs)      # 有些切分不带空格，再去掉空格拼一遍
def in_pdf(s):
    return s in TL_SP or s.replace(" ", "") in TL_NS
if _pdf_text is None:
    pass       # 已在 _load_pdftext 里记了一次 FAIL，这里不再刷屏
else:
    must(len(_recs) > 40, "PDF 文本层可取到文本操作数", len(_recs))
    for need in ["Mean log1p(CP10K)", "shared across all four panels",
                 "Pooled cells (all donors)", "% of cells expressing",
                 "Panel gene expression across cell types",
                 "SY_113", "SY_116", "SY_118"] + GENES + CTS:
        must(in_pdf(need), f"PDF 文本层含「{need}」")
    # 逐面板n_cells：图上 n 与源表必须一致（图注引的「肥大细胞 130 个」就在这里）
    for lbl, df_ in panels:
        p, m, n = GRID[lbl]
        for j, c in enumerate(CTS):
            must(f"n = {int(n.iloc[0, j]):,}" in TL_SP or f"n = {int(n.iloc[0, j])}" in TL_SP,
                 f"{lbl} / {c} 的 n={int(n.iloc[0, j])} 出现在图上")
    # 图内不得出现旧口径的痕迹：四面板同尺之后，色条标题必须自陈「shared」
    must("shared across all four panels" in TL_SP,
         "色条自陈四面板共享（R3-M1 的核心要求）")

# ── 5. 逐点数值落盘（审稿人可逐格回核） ──────────────────────────────────────
rows = []
for title, _ in panels:
    p, m, n = GRID[title]
    for g in GENES:
        for c in CTS:
            rows.append(dict(panel=title, cell_type=c, gene=g,
                             n_cells=int(n.loc[g, c]),
                             expr_fraction=round(float(p.loc[g, c]), 2),
                             mean_log1p=round(float(m.loc[g, c]), 4),
                             dot_area_pt2=round(area(float(p.loc[g, c])), 1)))
df_out = pd.DataFrame(rows)
must(len(df_out) == 120, "逐点表 4 面板 × 30 = 120 行", len(df_out))
# 逐点表与 R 脚本的审计表同落一处：交付包里是 重跑产出/源数据_EN/（与交付的源数据同名同格，
# 供逐件对账），工作树里直接写进交付目录的源数据（那里它就是交付件本身）。
_dst = ([os.path.join(OUTDIR, "源数据_EN")] if PACKAGE else
        [os.path.join(ROOT, s) for s in ("源数据_EN", "03_图表与源数据/源数据_EN")])
for d in _dst:
    if os.path.isdir(d):
        df_out.to_csv(os.path.join(d, "Fig5_values.csv"), index=False)
        print(f"  已写出 {os.path.join(d, 'Fig5_values.csv')}")
        break

print("\n" + "=" * 78)
if fails:
    print(f"结果：{len(fails)} 项未通过")
    for f in fails:
        print("  ✗ " + f)
    sys.exit(1)
print("结果：全部通过")
print(f"图：{FIG_PDF}")
