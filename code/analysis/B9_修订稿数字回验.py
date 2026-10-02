#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""B9 修订稿数字回验

用途：把修订稿中本批修改所写入的每一个数字，逐格回验到分析产物 CSV。
目的不是复核分析，而是复核**转写**——2026-09-17 的 B3 结果表曾出现
「PTPRC 一列读错位置」的转写错误，故每批稿件修改同样不靠目测，全部机器核对。

输入：03_分析/outputs/B3_scRNA/*.csv（§1–§8）
      02_结果表/*.csv（§11）+ 01_审计/A2c_*.md（§11 的全基因组背景率）
      05_稿件/两稿（§9 §10 §11 的字面串与数字提取）
输出：stdout（通过/不通过逐条），退出码 0=全通过

────────────────────────────────────────────────────────────────────────
2026-09-23 重定位版对脚本本身的四处改动（须留痕，不得静默）
────────────────────────────────────────────────────────────────────────
(1) 路径改指重定位后的两稿（见 MS / MS_EN）。旧稿两份**保留在原处**，
    本脚本不再指向它们；§12 以「新稿标题须存在、旧稿标题须不存在」作反向守卫，
    使「指向被改回去」或「误拿旧稿当新稿」都会失败，而不是静默通过。
(2) §9/§10 的字面串清单按重定位稿重写。原清单断言的是 2026-09-17 稿的写法
    （如 EN 侧要求出现 "more discriminating, not because its evidence is stronger"、
    "ambient RNA contamination"），重定位稿的说法不同，若照旧会假失败。
    **改法不是让稿件去迁就旧字面串**：把「本轮明确不得再写」的表述放进 absent 清单
    作回归守卫，把「本轮必须写到」的实质表述放进 present 清单。改动逐条列于
    ABSENT_ADDED / PRESENT_ADDED，见 §12。
(3) 新增 §11：本轮新写入正文的数字（判据敏感性的名次与跨线、留一法、82107 诊断、
    单细胞组成、ML 三套 AUC）此前**没有任何环节核对过转写**——§1–§8 只核对产物与
    硬编码常量，不读稿件正文。§11 从稿件正文按上下文锚点提取数值，再回验到产物。
(4) §12 登记本轮的字面串增补与路径变更，并加反向守卫。
────────────────────────────────────────────────────────────────────────
"""
import os
import re
import sys
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PROJ = os.path.dirname(BASE)                                          # …/滑膜基质刚度生成与刚度感知 0912
OUT = os.path.join(BASE, "03_分析", "outputs", "B3_scRNA")
TAB = os.path.join(BASE, "02_结果表")
AUD = os.path.join(BASE, "01_审计")
SCR = os.path.join(BASE, "03_分析", "scripts")
MS = os.path.join(BASE, "05_稿件", "双语初稿_候选优选与组织来源_重定位2026-09-23.md")
MS_EN = os.path.join(BASE, "05_稿件", "English_manuscript_MDPI_重定位2026-09-23.md")
OLD_EN = os.path.join(BASE, "05_稿件", "English_manuscript_MDPI_修订2026-09-17.md")

CN2EN = {"成纤维细胞": "Fibroblast", "巨噬细胞": "Macrophage",
         "T 细胞": "T cell", "内皮细胞": "Endothelial", "肥大细胞": "Mast"}

# ---- §12 登记：本轮对字面串清单的增补（每项都写明理由）----
ABSENT_ADDED = {
    # 方案 §二 读数 3 的原始表述有误（LOXL4 的 |logFC| 小于 LOXL1，只有 P 更小），
    # 退回该表述即说明误读了核心表
    "LOXL4 的效应量大于 LOXL1": "方案读数 3 的误表述；实测 0.596 < 0.638",
    "LOXL4's effect is larger than LOXL1": "同上（英文侧）",
    # 82107 的诊断结论不得退回「低功效故可忽略」
    "低功效故可忽略": "A2c 已否证：该队列对这批基因并非无信息（P=0.035）",
    "systematically unreliable": "A2c 已否证：分歧多为测不出而非相矛盾",
    # 未找到独立病例—对照滑膜 scRNA（A-2.5），不得声称外部单细胞验证
    "外部单细胞验证": "A-2.5 检索后未找到，记为未完成",
    "独立单细胞验证": "同上",
    "independent single-cell validation": "同上（英文侧）",
    # 定位数据不得升格为功能或确定性断言
    "系统性地不可靠": "与 A2c 结论相反",
    "确定性地定位": "定位是比例，不是确定性",
}
PRESENT_ADDED_CN = [
    "判据依赖", "条件性排序", "不可合并", "不提供信息", "未做预注册",
    "43.65 个百分点", "0.0601 与 0.0586", "按 |log2FC| 排第 23，按 P 排第 8",
    "PLOD1 是唯一在两层随机效应下均显著的交联基因",
    "只有 LOXL1 与 SERPINE1 在每一次单队列剔除下保持方向且 P<0.01",
]
PRESENT_ADDED_EN = [
    "criterion-dependent", "conditional ordering", "not combinable",
    "not pre-registered", "43.65 percentage points",
    "0.0601 and 0.0586", "23rd by |log2FC| but 8th by P",
    "PLOD1 is the only cross-linker",
    "only LOXL1 and SERPINE1 preserve both direction and P<0.01",
    "not attributable to ambient RNA",
]
OVERRIDE = 0.05 / 7          # 面板 Bonferroni 阈值分母＝实际检验的 7 个基因

res = []


def chk(name, got, want, tol=0.005):
    """got 为实测值，want 为稿件中写的值；tol 为允许的舍入差。"""
    ok = abs(float(got) - float(want)) <= tol
    res.append((ok, name, got, want))


# ────────────────────────── §11 用的提取工具 ──────────────────────────
SUP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
       "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9"}
MINUS = {"−": "-", "–": "-", "—": "-"}


def tonum(s):
    """把稿件里的数值转成 float：支持 `0.0601`、`−0.596`、`3.9×10⁻⁴`、`21,355`。
    注意上标负号是 U+207B，不在 SUP 表内，须先取出定符号。"""
    s = s.strip().replace(",", "")
    for a, b in MINUS.items():
        s = s.replace(a, b)
    s = s.replace("+", "")
    m = re.match(r"^([-+]?[\d.]+)\s*[×x]\s*10([⁻⁺⁰¹²³⁴⁵⁶⁷⁸⁹]+)$", s)
    if m:
        e = m.group(2)
        sign = -1 if "⁻" in e else 1
        exp = sign * int("".join(SUP[c] for c in e.replace("⁻", "").replace("⁺", "")))
        return float(m.group(1)) * (10 ** exp)
    return float(s)


def grab(text, pattern, name, grp=1):
    """按上下文锚点从稿件提取一个数；锚点匹配不到即判失败（防句子被删后空转）。"""
    m = re.search(pattern, text)
    if not m:
        res.append((False, f"[提取失败] {name}", "锚点未匹配", repr(pattern[:60])))
        return None
    try:
        return tonum(m.group(grp))
    except ValueError:
        res.append((False, f"[提取失败] {name}", m.group(grp), "非数值"))
        return None


def transcribe(text, pattern, name, want, grp=1, tol=0.005, prefix="EN"):
    v = grab(text, pattern, f"{prefix} 转写 {name}", grp)
    if v is not None:
        chk(f"{prefix} 转写 {name}", v, want, tol)


def main():
    # ═══════════════════════ §1–§8：B3 产物核对（未改动） ═══════════════════════
    pct = pd.read_csv(os.path.join(OUT, "B3_面板表达_逐细胞类型.csv"))
    mean = pd.read_csv(os.path.join(OUT, "B3_平均表达口径排查.csv"))
    per_d = pd.read_csv(os.path.join(OUT, "B3_面板表达_逐供者逐细胞类型.csv"))
    small = pd.read_csv(os.path.join(OUT, "B3_小样本敏感性.csv"))
    comp = pd.read_csv(os.path.join(OUT, "B3_逐供者细胞类型_比例.csv"))
    scr = pd.read_csv(os.path.join(OUT, "B3_Scrublet逐供者.csv"))
    cmp13 = pd.read_csv(os.path.join(OUT, "B3_与正文对照.csv"))

    def P(gene, ct_en):
        return float(pct[(pct["基因"] == gene) &
                         (pct["cell_type"] == ct_en)]["表达比例"].iloc[0])

    def M(gene, ct_cn):
        return float(mean[(mean["基因"] == gene) &
                          (mean["细胞类型"] == ct_cn)]["该口径值"].iloc[0])

    def PD(gene, ct_en, donor):
        return float(per_d[(per_d["基因"] == gene) & (per_d["cell_type"] == ct_en) &
                           (per_d["donor"] == donor)]["表达比例"].iloc[0])

    def S(ct_cn, gene, col):
        return float(small[(small["细胞类型"] == ct_cn) &
                           (small["基因"] == gene)][col].iloc[0])

    # ---- 1. 正文写入的比例 ----
    # LOXL1 保留原文 12.80（实测 12.65，差 0.15pp，稿件中已说明复现至 1.3pp 以内）
    chk("LOXL1 成纤维 比例=12.80", P("LOXL1", "Fibroblast"), 12.80, tol=0.20)
    chk("LOXL1 内皮 1.64", P("LOXL1", "Endothelial"), 1.64)
    chk("LOXL1 巨噬 2.73", P("LOXL1", "Macrophage"), 2.73)
    chk("LOXL1 肥大 3.08", P("LOXL1", "Mast"), 3.08)
    chk("LOXL1 T 细胞 6.08", P("LOXL1", "T cell"), 6.08)
    chk("SDC4 成纤维 37.48", P("SDC4", "Fibroblast"), 37.48, tol=0.20)
    chk("SDC4 巨噬 20.27", P("SDC4", "Macrophage"), 20.27)
    chk("SDC4 T 细胞 19.34", P("SDC4", "T cell"), 19.34)
    chk("VGLL4 内皮 65.45", P("VGLL4", "Endothelial"), 65.45, tol=1.30)
    chk("VGLL4 成纤维 45.31", P("VGLL4", "Fibroblast"), 45.31, tol=0.20)
    chk("VGLL4 T 细胞 24.86", P("VGLL4", "T cell"), 24.86)
    chk("VGLL4 巨噬 15.53", P("VGLL4", "Macrophage"), 15.53)
    chk("SERPINE1 内皮 69.10", P("SERPINE1", "Endothelial"), 69.10, tol=0.90)
    chk("SERPINE1 成纤维 39.60", P("SERPINE1", "Fibroblast"), 39.60, tol=0.20)
    chk("SERPINE1 巨噬 23.22", P("SERPINE1", "Macrophage"), 23.22)
    chk("SERPINE1 T 细胞 17.13", P("SERPINE1", "T cell"), 17.13)

    # ---- 2. 平均表达（全细胞 mean CP10K）----
    chk("LOXL1 成纤维 平均 0.15", M("LOXL1", "成纤维细胞"), 0.15, tol=0.01)
    chk("SERPINE1 内皮 平均 5.32", M("SERPINE1", "内皮细胞"), 5.32, tol=0.03)
    chk("SERPINE1 成纤维 平均 0.90", M("SERPINE1", "成纤维细胞"), 0.90, tol=0.02)
    chk("SERPINE1 巨噬 平均 0.49", M("SERPINE1", "巨噬细胞"), 0.49, tol=0.05)

    # 稿件写「其余四类 1.6%–6.1%」「平均 0.009–0.066 CP10K」——须为真实上下界
    others = ["Endothelial", "Macrophage", "Mast", "T cell"]
    lo, hi = min(P("LOXL1", c) for c in others), max(P("LOXL1", c) for c in others)
    chk("LOXL1 其余四类 下界 1.6", lo, 1.6, tol=0.05)
    chk("LOXL1 其余四类 上界 6.1", hi, 6.1, tol=0.02)
    mo = [float(pct[(pct["基因"] == "LOXL1") & (pct["cell_type"] == c)]["平均CP10K"].iloc[0])
          for c in others]
    chk("LOXL1 其余四类 平均下界 0.009", min(mo), 0.009, tol=0.0005)
    chk("LOXL1 其余四类 平均上界 0.066", max(mo), 0.066, tol=0.001)

    # 稿件写「成纤维较其余类型高 2–8 倍」
    fib = P("LOXL1", "Fibroblast")
    chk("LOXL1 富集倍数 下界 2", fib / hi, 2, tol=0.15)
    chk("LOXL1 富集倍数 上界 8", fib / lo, 8, tol=0.35)

    # 稿件写 SERPINE1「其余三类 17%–40%」（成纤维/巨噬/T）
    o3 = [P("SERPINE1", c) for c in ["Fibroblast", "Macrophage", "T cell"]]
    chk("SERPINE1 其余三类 下界 17", min(o3), 17.1, tol=0.10)
    chk("SERPINE1 其余三类 上界 40", max(o3), 39.6, tol=0.15)

    # ---- 3. 逐供者稳定性 ----
    for d, v in [("SY_113", 10.78), ("SY_116", 14.12), ("SY_118", 17.20)]:
        chk(f"LOXL1 成纤维 {d}", PD("LOXL1", "Fibroblast", d), v, tol=0.02)
    for d, v in [("SY_113", 28.09), ("SY_116", 11.11), ("SY_118", 4.62)]:
        chk(f"SERPINE1 T {d}", PD("SERPINE1", "T cell", d), v, tol=0.02)

    # ---- 4. 明确细胞子集 ----
    chk("LOXL1 成纤维 明确细胞 12.74", S("成纤维细胞", "LOXL1", "明确细胞_表达比例"),
        12.74, tol=0.02)
    chk("SERPINE1 T 明确细胞 18.97", S("T 细胞", "SERPINE1", "明确细胞_表达比例"),
        18.97, tol=0.02)

    # ---- 5. 供者细胞构成 ----
    for d, v in [("SY_113", 86.83), ("SY_116", 87.54), ("SY_118", 43.89)]:
        got = float(comp[comp["donor"] == d]["Fibroblast"].iloc[0])
        chk(f"{d} 成纤维占比 {v}", got, v, tol=0.02)
    chk("SY_118 巨噬占比 40.64",
        float(comp[comp["donor"] == "SY_118"]["Macrophage"].iloc[0]), 40.64, tol=0.02)
    chk("SY_118 内皮占比 10.80",
        float(comp[comp["donor"] == "SY_118"]["Endothelial"].iloc[0]), 10.80, tol=0.02)

    # ---- 6. Scrublet ----
    chk("SY_113 双细胞 0", float(scr[scr["供者"] == "SY_113"]["检出双细胞"].iloc[0]), 0)
    chk("SY_116 双细胞 0", float(scr[scr["供者"] == "SY_116"]["检出双细胞"].iloc[0]), 0)
    chk("SY_118 双细胞 1", float(scr[scr["供者"] == "SY_118"]["检出双细胞"].iloc[0]), 1)
    chk("SY_118 双细胞比例 0.06",
        float(scr[scr["供者"] == "SY_118"]["双细胞比例"].iloc[0]), 0.06)

    # ---- 7. 「比例复现至 1.3pp」这句本身 ----
    ok13 = cmp13[cmp13["指标"] == "表达比例"]
    tgt = ok13[ok13["细胞类型"].isin(["成纤维细胞", "内皮细胞"])]
    bad = tgt[tgt["判定"] != "一致"]
    res.append((len(bad) == 0,
                "涉及成纤维/内皮的表达比例全部判定为一致",
                f"{len(bad)} 条不一致", "0 条"))
    mx = tgt["差"].abs().max()
    res.append((mx <= 1.3, "上述比例最大偏差 ≤1.3pp", round(mx, 2), "≤1.3"))
    okm = cmp13[cmp13["指标"] == "平均CP10K"]
    mxm = okm["差"].abs().max()
    res.append((mxm <= 0.05, "平均表达最大偏差 ≤0.05 CP10K", round(mxm, 3), "≤0.05"))
    res.append((len(okm) == 4, "平均表达对照共 4 条", len(okm), 4))
    nfail = len(cmp13[cmp13["判定"] != "一致"])
    res.append((nfail == 3, "逐值对照不一致恰为 3 条", nfail, 3))

    # ---- 8. 细胞总数 ----
    tot = int(pct.groupby("cell_type")["细胞数"].first().sum())
    res.append((tot == 10600, "四类+肥大 细胞合计 = 10,600", tot, 10600))

    # ═══════════════════════ §9 §10：字面串（本轮重写） ═══════════════════════
    txt = open(MS, encoding="utf-8").read()
    txt_en = open(MS_EN, encoding="utf-8").read()

    for pat in ["其余均 <1.6%", "免疫细胞 <1%", "免疫 <1%", "免疫 ≤13%",
                "均 <1.6%、<0.022", "的表达限于成纤维细胞",
                "通过 18 个组合校正", "通过跨组合校正"] + list(ABSENT_ADDED):
        res.append((pat not in txt, f"[CN] 旧/禁用断言已消失：{pat!r}",
                    "仍存在" if pat in txt else "无", "无"))
    for pat in ["其余四类细胞 1.6%–6.1%", "全细胞平均 CP10K", "43.89%", "18.97%"] \
            + PRESENT_ADDED_CN:
        res.append((pat in txt, f"[CN] 新表述已写入：{pat!r}",
                    "有" if pat in txt else "缺", "有"))

    EN_ABSENT = ["12.8% vs <1.6% elsewhere", "12.80% vs <1% immune cells",
                 "immune ≤13%", "37.48% vs 17.71% macrophages",
                 "<1.6% and <0.022 in every other cell type",
                 "its cellular localization is determinate",
                 "the **determinacy of its cellular localization",
                 "four annotated cell types;",
                 "more discriminating, not because its evidence is stronger",
                 "simply underpowered", "LOXL4's effect is larger than LOXL1"]
    for pat in EN_ABSENT:
        res.append((pat not in txt_en, f"[EN] 旧/禁用断言已消失：{pat!r}",
                    "仍存在" if pat in txt_en else "无", "无"))
    EN_PRESENT = ["1.6%–6.1% in the other four cell types",
                  "whole-cell mean CP10K", "43.89%", "18.97%"] + PRESENT_ADDED_EN
    for pat in EN_PRESENT:
        res.append((pat in txt_en, f"[EN] 新表述已写入：{pat!r}",
                    "有" if pat in txt_en else "缺", "有"))

    # ═══════ §11：本轮新写入正文的数字 —— 从正文提取，回验到 A2/A2c/A2e/B6 产物 ═══════
    cross = pd.read_csv(os.path.join(TAB, "A2_交联家族_判据敏感性.csv")).set_index("gene")
    over = pd.read_csv(os.path.join(TAB, "57基因总览.csv")).set_index("gene")
    perm = pd.read_csv(os.path.join(TAB, "A2c_置换零分布.csv"))
    diag = pd.read_csv(os.path.join(TAB, "A2c_队列级诊断.csv")).set_index("cohort")
    cls = pd.read_csv(os.path.join(TAB, "A2c_82107_方向相左分类.csv"))
    per = pd.read_csv(os.path.join(TAB, "A2c_57基因_逐队列.csv"))
    # A2d 是 cell_type 为行的宽表（供者在列），不是 donor 为行
    comp_don = pd.read_csv(os.path.join(TAB, "A2d_供者细胞组成.csv")).set_index("cell_type")
    ml12 = pd.read_csv(os.path.join(TAB, "F11下游_ml12_定标新旧对照.csv")).set_index("队列")
    ml11 = pd.read_csv(os.path.join(TAB, "F11下游_ml11_GSE89408行.csv")).iloc[0]
    ml10 = pd.read_csv(os.path.join(TAB, "F11下游_ml10_新旧对照.csv"), header=None)
    deltap = pd.read_csv(os.path.join(TAB, "A2_显著性跨线计数.csv"))
    loo = pd.read_csv(os.path.join(TAB, "A2e_留一法_逐基因名次变动.csv"))
    aud = open(os.path.join(AUD, "A2c_全基因组一致性诊断_结果.md"), encoding="utf-8").read()

    # 11.1 判据敏感性的核心数字（发现层）
    transcribe(txt_en, r"under random effects, SDC4 falls from [^\s]+ to P=([\d.]+)",
               "发现层 SDC4 的 RE-DL P", over.loc["SDC4", "disc_re_P"])
    transcribe(txt_en, r"ITGA5 from 0\.022 to P=([\d.]+)",
               "发现层 ITGA5 的 RE-DL P", over.loc["ITGA5", "disc_re_P"])
    transcribe(txt_en, r"heterogeneity is the highest in the panel \(I²=([\d.]+) and ([\d.]+)\)",
               "发现层 SDC4 的 I²", over.loc["SDC4", "disc_I2"], tol=0.05)
    transcribe(txt_en, r"heterogeneity is the highest in the panel \(I²=[\d.]+ and ([\d.]+)\)",
               "发现层 ITGA5 的 I²", over.loc["ITGA5", "disc_I2"], tol=0.05, grp=1)
    transcribe(txt_en, r"LOXL1 \+([\d.]+), P=([^\s]+) in discovery",
               "发现层 LOXL1 的 log2FC", over.loc["LOXL1", "disc_logFC"])
    transcribe(txt_en, r"LOXL1 \+[\d.]+, P=([^\s]+) in discovery",
               "发现层 LOXL1 的 IVW-FE P", over.loc["LOXL1", "disc_P"])
    transcribe(txt_en,
               r"in discovery and \+([\d.]+), P=([\d.]+) in replication",
               "复制层 LOXL1 的 pooled log2FC", over.loc["LOXL1", "rep_logFC"])
    transcribe(txt_en,
               r"in discovery and \+[\d.]+, P=([\d.]+) in replication",
               "复制层 LOXL1 的 pooled P", over.loc["LOXL1", "rep_P"])
    transcribe(txt_en, r"LOX \+([\d.]+), P=([\d.]+) and",
               "发现层 LOX 的 log2FC", over.loc["LOX", "disc_logFC"])
    transcribe(txt_en, r"LOX \+[\d.]+, P=([\d.]+) and",
               "发现层 LOX 的 IVW-FE P", over.loc["LOX", "disc_P"])

    # 11.2 复制层 pooled 全集 + RE-DL 反转（「条件性首位」的核心证据）
    transcribe(txt_en, r"under it LOXL1 \(P=([\d.]+)\)", "复制层 LOXL1 的 pooled P",
               over.loc["LOXL1", "rep_P"])
    for g, pat in [("SDC4", r"SDC4 \(([\d.]+)\), SAV1"),
                   ("SAV1", r"SAV1 \(([\d.]+)\), ITGA5"),
                   ("ITGA5", r"ITGA5 \(([\d.]+)\), SERPINE1"),
                   ("SERPINE1", r"SERPINE1 \(([\d.]+)\) and VGLL4"),
                   ("VGLL4", r"VGLL4 \(([\d.]+)\) all reach")]:
        transcribe(txt_en, pat, f"复制层 {g} 的 pooled P", over.loc[g, "rep_P"])
    transcribe(txt_en, r"LOXL1 rises to P=([\d.]+)", "复制层 LOXL1 的 RE-DL P",
               over.loc["LOXL1", "rep_re_P"])
    transcribe(txt_en, r"both lose replication significance \(P=0\.0601 and ([\d.]+)\)",
               "复制层 LOX 的 RE-DL P", over.loc["LOX", "rep_re_P"])
    for g, pat in [("SDC4", r"SDC4 to ([\d.]+),"), ("VGLL4", r"VGLL4 to ([\d.]+),"),
                   ("ITGA5", r"ITGA5 to ([\d.]+) and"), ("SAV1", r"SAV1 to ([\d.]+),")]:
        transcribe(txt_en, pat, f"复制层 {g} 的 RE-DL P", over.loc[g, "rep_re_P"])
    # 注意：同一句里 pooled 值写作 `SERPINE1 (0.0254)`，RE-DL 值写作 `SERPINE1 (P=0.0351)`，
    # 故 RE-DL 的锚点必须带 `P=`，否则会静默取到 pooled 那一格。
    transcribe(txt_en, r"SERPINE1 \(P=([\d.]+)\) and, among the cross-linkers",
               "复制层 SERPINE1 的 RE-DL P", over.loc["SERPINE1", "rep_re_P"])
    transcribe(txt_en, r"PLOD1 \(P=([\d.]+)\) retain it", "复制层 PLOD1 的 RE-DL P",
               over.loc["PLOD1", "rep_re_P"])
    transcribe(txt_en, r"LOX loses it in discovery as well \(P=([\d.]+)\)",
               "发现层 LOX 的 RE-DL P", over.loc["LOX", "disc_re_P"])

    # 11.3 PLOD1 = 两层 RE 均显著；名次随判据变
    transcribe(txt_en, r"\(discovery P=([\d.]+) and replication P=([\d.]+), I²=0 in both\)",
               "发现层 PLOD1 的 RE-DL P", over.loc["PLOD1", "disc_re_P"])
    transcribe(txt_en, r"\(discovery P=[\d.]+ and replication P=([\d.]+), I²=0 in both\)",
               "复制层 PLOD1 的 RE-DL P", over.loc["PLOD1", "rep_re_P"])
    transcribe(txt_en, r"it is (\d+)rd by \|log2FC\| but (\d+)th by P",
               "复制层 PLOD1 的 |log2FC| 名次", over.loc["PLOD1", "rank_abs_rep"])
    transcribe(txt_en, r"it is \d+rd by \|log2FC\| but (\d+)th by P",
               "复制层 PLOD1 的 P 名次", over.loc["PLOD1", "rank_P_rep"])

    # 11.4 LOXL4：家族中最小 P 但方向相反且不可合并（读数 3 的正确写法）
    transcribe(txt_en, r"LOXL4 reaches P=([^\s]+) in the discovery layer",
               "发现层 LOXL4 的 IVW-FE P", over.loc["LOXL4", "disc_P"])
    transcribe(txt_en, r"smaller than LOXL1's ([^\s]+)",
               "发现层 LOXL1 的 IVW-FE P（对照）", over.loc["LOXL1", "disc_P"])
    transcribe(txt_en, r"opposite direction \(([^)]+)\)",
               "发现层 LOXL4 的 log2FC", over.loc["LOXL4", "disc_logFC"])
    transcribe(txt_en, r"larger in both layers \(([^)]+) discovery, ([^)]+) replication\)",
               "发现层 SERPINE1 的 log2FC", over.loc["SERPINE1", "disc_logFC"])
    transcribe(txt_en, r"larger in both layers \([^)]+ discovery, ([^)]+) replication\)",
               "复制层 SERPINE1 的 pooled log2FC", over.loc["SERPINE1", "rep_logFC"], grp=1)
    # 跨线计数：本轮实测发现层 2 个、复制层 3 个，正文未逐一列举，仅核对计数表自洽
    res.append((int(deltap.iloc[0]["n_crossed"]) == 2,
                "发现层 FE↔RE 跨线计数 = 2（LOX、PLOD2）",
                int(deltap.iloc[0]["n_crossed"]), 2))
    res.append((int(deltap.iloc[1]["n_crossed"]) == 3,
                "复制层 pooled↔RE 跨线计数 = 3（LOX、LOXL1、TGM2）",
                int(deltap.iloc[1]["n_crossed"]), 3))

    # 11.5 留一队列剔除（数值一律从 A2e 产物读取，不写死）
    DISC_EST, REP_EST = "IVW-FE", "pooled(~0+group+cohort)"

    def looP(gene, dropped, layer="发现", est=DISC_EST):
        sub = loo[(loo["layer"] == layer) & (loo["estimator"] == est) &
                  (loo["gene"] == gene) & (loo["dropped"] == dropped)]
        return float(sub["P"].iloc[0])

    transcribe(txt_en, r"LOX loses significance when GSE89408 is removed "
                       r"\(P=[^\s]+ → ([\d.]+)\)",
               "剔除 GSE89408 后 LOX 的 P", looP("LOX", "GSE89408"), tol=0.001)
    transcribe(txt_en, r"PLOD1 when GSE206848 is removed \(→([\d.]+)\)",
               "剔除 GSE206848 后 PLOD1 的 P", looP("PLOD1", "GSE206848"), tol=0.001)
    transcribe(txt_en, r"PLOD2 when either GSE89408 or GSE1919 is removed "
                       r"\(→([\d.]+) and ([\d.]+)\)",
               "剔除 GSE89408 后 PLOD2 的 P", looP("PLOD2", "GSE89408"), tol=0.001)
    transcribe(txt_en, r"PLOD2 when either GSE89408 or GSE1919 is removed "
                       r"\(→[\d.]+ and ([\d.]+)\)",
               "剔除 GSE1919 后 PLOD2 的 P", looP("PLOD2", "GSE1919"), tol=0.001)

    def flipped(dropped, layer="发现", est=DISC_EST, full_col="disc_logFC"):
        """剔除后的估计与全集主估计符号相反的面板/交联基因。
        该表没有全集列，故全集效应量取自 57基因总览 的对应层主估计列。"""
        sub = loo[(loo["layer"] == layer) & (loo["estimator"] == est) &
                  (loo["dropped"] == dropped) & (loo["in_panel"] | loo["is_crosslink"])]
        if sub.empty:
            return []
        sub = sub.assign(full=over.loc[sub["gene"], full_col].values)
        return sorted(sub[sub["logFC"] * sub["full"] < 0]["gene"])

    fl19 = flipped("GSE1919")
    res.append((fl19 == ["ITGA5", "LOXL2", "TGM2"],
                "剔除 GSE1919 时翻号的面板/交联基因 = ITGA5/LOXL2/TGM2",
                "、".join(fl19), "ITGA5、LOXL2、TGM2"))
    for d in ["GSE89408", "GSE206848"]:
        fl = flipped(d)
        res.append((fl == [], f"发现层剔除 {d} 时无面板/交联基因翻号",
                    "、".join(fl) or "无", "无"))
    rep_flip = set()
    for d in loo[loo["layer"] == "复制"]["dropped"].unique():
        rep_flip |= set(flipped(d, layer="复制", est=REP_EST, full_col="rep_logFC"))
    res.append((not rep_flip, "复制层任何单队列剔除均无面板/交联基因翻号",
                "、".join(sorted(rep_flip)) or "无", "无"))

    # 11.6 GSE82107 诊断（不得写成「低功效故可忽略」）
    disc = cls[cls["class"] != "方向一致"]
    res.append((len(disc) == 17, "82107 方向相左候选数 = 17", len(disc), 17))
    res.append(((cls["class"] == "a_有功效仍反向").sum() == 1,
                "82107「有功效仍反向」= 1（CYR61）",
                int((cls["class"] == "a_有功效仍反向").sum()), 1))
    res.append(((cls["class"] == "b_标称显著却反向").sum() == 3,
                "82107「标称显著却反向」= 3",
                int((cls["class"] == "b_标称显著却反向").sum()), 3))
    res.append(((cls["class"] == "c_测不出_不提供信息").sum() == 13,
                "82107「测不出」= 13",
                int((cls["class"] == "c_测不出_不提供信息").sum()), 13))
    cyr = cls[cls["gene"] == "CYR61"].iloc[0]
    cyr_pat = (r"CCN1/CYR61 \(\+([\d.]+) versus ([^\s,]+), exceeding that cohort's "
               r"MDE of ([\d.]+)\)")
    transcribe(txt_en, cyr_pat, "82107 CYR61 的反向 logFC", cyr["logFC"])
    transcribe(txt_en, cyr_pat, "82107 CYR61 的参照层 logFC", cyr["ref_logFC"], grp=2)
    transcribe(txt_en, cyr_pat, "82107 CYR61 的 MDE", cyr["MDE"], grp=3)
    # 「标称显著却反向」的基因名须与分类表逐名一致
    m_b = re.search(r"nominally significant yet reversed \(([^)]+)\) and thirteen",
                    txt_en)
    if m_b is None:
        res.append((False, "[提取失败] EN 中「标称显著却反向」基因名", "锚点未匹配", "有"))
    else:
        got_b = sorted(g.strip() for g in m_b.group(1).replace(" and ", ", ").split(",")
                       if g.strip())
        want_b = sorted(cls[cls["class"] == "b_标称显著却反向"]["gene"])
        res.append((got_b == want_b, "EN「标称显著却反向」基因名与分类表一致",
                    "、".join(got_b), "、".join(want_b)))
    n_same = int((cls["class"] == "方向一致").sum())
    n_eval = len(cls)
    transcribe(txt_en, r"for these 52 candidates, ([\d.]+)%, against a label-permutation "
                       r"null whose median is ([\d.]+) of 52 \(empirical P=([\d.]+)\)",
               "82107 的 52 基因同号率（表 2 脚注）", n_same / n_eval * 100, tol=0.05)
    # 同一组数字在 §3.2 正文又以另一种句式出现一次，必须一并核对：
    # 只核一处时，改掉另一处仍会全绿（负面对照已实测出这一漏洞）。
    transcribe(txt_en, r"for the 52 candidates it is \*\*([\d.]+)%\*\* \((\d+) of 52, "
                       r"against a null median of ([\d.]+); empirical P=([\d.]+)\)",
               "82107 的 52 基因同号率（§3.2 正文）", n_same / n_eval * 100, tol=0.05)
    transcribe(txt_en, r"for the 52 candidates it is \*\*[\d.]+%\*\* \((\d+) of 52, "
                       r"against a null median of ([\d.]+); empirical P=([\d.]+)\)",
               "82107 的 52 基因同号数（§3.2 正文）", n_same, tol=0, grp=1)
    n_ge = int((perm["g57"] >= n_same).sum())
    p57 = (n_ge + 1) / (len(perm) + 1)
    transcribe(txt_en, r"null whose median is ([\d.]+) of 52 \(empirical P=([\d.]+)\)",
               "82107 57 基因置换零分布中位", float(perm["g57"].median()), tol=0.05)
    transcribe(txt_en, r"against a null median of ([\d.]+); empirical P=([\d.]+)\)",
               "82107 57 基因置换零分布中位（§3.2 正文）", float(perm["g57"].median()),
               tol=0.05)
    transcribe(txt_en, r"against a null median of [\d.]+; empirical P=([\d.]+)\)",
               "82107 57 基因置换经验 P（§3.2 正文）", p57, tol=0.001, grp=1)
    transcribe(txt_en, r"null whose median is [\d.]+ of 52 \(empirical P=([\d.]+)\)",
               "82107 57 基因置换经验 P", p57, tol=0.001)
    # 57 基因子集的两项同样须与 A2c 审计记录对账（该记录此处是散文而非表行）
    row57 = re.search(r"57 基因子集（可评 (\d+) 基因）[^；]*；与零分布（([\d.]+) 同号）"
                      r"相比经验 ?P ?= ?([\d.]+)", aud)
    if row57 is None:
        res.append((False, "[提取失败] A2c 审计记录的 57 基因行", "锚点未匹配", "有"))
    else:
        chk("A2c 审计记录 57 基因可评数 = 52", float(row57.group(1)), n_eval, tol=0)
        # 该句里的第二个数是零分布中位（26.0 同号），不是实测同号率（实测另见分类表）
        chk("A2c 审计记录 57 基因零分布中位", float(row57.group(2)),
            float(perm["g57"].median()), tol=0.05)
        chk("A2c 审计记录 57 基因经验 P 可由置换 CSV 重现", p57, float(row57.group(3)),
            tol=0.0006)
    # 全基因组背景率：置换 CSV 无实测列，故实测值取自 A2c 审计记录，
    # 但必须先与置换 CSV 对账（P 由 CSV 现算须重现审计记录的 P），不能只抄记录。
    row_all = re.search(r"\|\s*全基因组交集\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*([\d.]+)%"
                        r"\s*\|[^|]*\|\s*([\d.]+)\s*\|", aud)
    if row_all is None:
        res.append((False, "[提取失败] A2c 审计记录的全基因组行", "锚点未匹配", "有"))
        obs_all, n_all, rate_all, p_all = 7272, 13039, 55.8, 0.025
    else:
        n_all, obs_all, rate_all, p_all = (int(row_all.group(1)), int(row_all.group(2)),
                                           float(row_all.group(3)), float(row_all.group(4)))
    chk("A2c 审计记录 全基因组同号率 = 实测/交集", obs_all / n_all * 100, rate_all, tol=0.05)
    n_ge_all = int((perm["all"] >= obs_all).sum())
    p_all_csv = (n_ge_all + 1) / (len(perm) + 1)
    chk("A2c 审计记录的 全基因组经验 P 可由置换 CSV 重现", p_all_csv, p_all, tol=0.0006)
    # 正文该句括号内还有 "against a label-permutation null"，故锚点须止于 "against"
    transcribe(txt_en, r"genome-wide sign concordance with the other two cohorts is "
                       r"([\d.]+)% \(empirical P=([\d.]+) against",
               "82107 全基因组同号率", rate_all, tol=0.05)
    transcribe(txt_en, r"genome-wide sign concordance with the other two cohorts is "
                       r"[\d.]+% \(empirical P=([\d.]+) against",
               "82107 全基因组置换经验 P", p_all, tol=0.0006, grp=1)
    plat = r"\(([\d,]+) features vs ([\d,]+)\)"
    transcribe(txt_en, plat, "82107 平台特征数", diag.loc["GSE82107", "n_gene_platform"])
    transcribe(txt_en, plat, "GPL96 平台特征数", diag.loc["GSE55235", "n_gene_platform"],
               grp=2)
    n82107 = int(per[(per["cohort"] == "GSE82107") & (per["P"] < 0.05)].shape[0])
    n55457 = int(per[(per["cohort"] == "GSE55457") & (per["P"] < 0.05)].shape[0])
    transcribe(txt_en, r"only (\d+) of the 57 candidates reach nominal P<0.05 there, "
                       r"against (\d+) in GSE55457",
               "82107 名义 P<0.05 的候选数", n82107)
    transcribe(txt_en, r"only \d+ of the 57 candidates reach nominal P<0.05 there, "
                       r"against (\d+) in GSE55457",
               "GSE55457 名义 P<0.05 的候选数", n55457, grp=1)
    # 两个 GPL96 队列的 MDE 必须**各自**取锚点。此前的写法是
    #   want=max(GSE55235, GSE55457), tol=0.002
    # 而两队列实测为 0.641999 与 0.642740（差 0.0007），故 |实测−max| 恒小于 tol：
    # 稿件把其中一个队列的值写成「两个队列共用一个值」时它照样全绿——该断言结构上不可能失败。
    # 现改为锚点缺任一数即判失败，且逐队列分别核对。
    transcribe(txt_en, r"genome-wide median minimum detectable effect is ([\d.]+) log2FC ",
               "82107 中位 MDE（全基因组）", diag.loc["GSE82107", "MDE_logFC_med"],
               tol=0.0005)
    transcribe(txt_en,
               r"log2FC against ([\d.]+) and ([\d.]+) in the two GPL96 cohorts",
               "GSE55235 中位 MDE", diag.loc["GSE55235", "MDE_logFC_med"],
               tol=0.0005, grp=1)
    transcribe(txt_en,
               r"log2FC against [\d.]+ and ([\d.]+) in the two GPL96 cohorts",
               "GSE55457 中位 MDE", diag.loc["GSE55457", "MDE_logFC_med"],
               tol=0.0005, grp=1)
    # 同一句给出的 d 尺度值也须核对：审计记录（A2c 结果 §一 表注）自己规定
    # 「跨队列比功效必须用 d 尺度」，故正文若只给 log2FC 尺度值即为漏报该尺度。
    transcribe(txt_en, r"\(([\d.]+) against [\d.]+ on the common d scale\)",
               "82107 中位 MDE（d 尺度）", diag.loc["GSE82107", "MDE_d_med"], tol=0.005)
    # 两 GPL96 队列的 d 尺度值在产物中完全相同（均 1.252907），故正文共用一个值是正确的，
    # 这里也只需核对一次——与上面 log2FC 尺度必须逐队列分开的情形不同。
    transcribe(txt_en, r"\([\d.]+ against ([\d.]+) on the common d scale\)",
               "GPL96 中位 MDE（d 尺度）", diag.loc["GSE55235", "MDE_d_med"],
               tol=0.005, grp=1)

    # 11.7 单细胞：组成极差与逐供者比例（文本侧）
    fib_row = comp_don.loc["Fibroblast"]          # 行=细胞类型，列=供者
    spread = float(fib_row.max() - fib_row.min())
    transcribe(txt_en, r"a spread of \*\*([\d.]+) percentage points\*\*",
               "成纤维细胞占比极差", spread, tol=0.02)
    transcribe(txt_en, r"in fibroblasts \(([\d.]+)%, against 1\.6%–6\.1%",
               "LOXL1 成纤维比例（文本侧）", P("LOXL1", "Fibroblast"), tol=0.20)
    transcribe(txt_en, r"fibroblasts \(([\d.]+)%, against 20\.27%",
               "SDC4 成纤维比例（文本侧）", P("SDC4", "Fibroblast"), tol=0.20)
    transcribe(txt_en, r"endothelial cells \(([\d.]+)%\) but is present in the other "
                       r"three at ([\d.]+)%–([\d.]+)%",
               "SERPINE1 内皮比例（文本侧）", P("SERPINE1", "Endothelial"), tol=0.90)
    transcribe(txt_en, r"but is present in the other three at ([\d.]+)%–[\d.]+%",
               "SERPINE1 其余三类下界", min(o3), tol=0.10)
    transcribe(txt_en, r"but is present in the other three at [\d.]+%–([\d.]+)%",
               "SERPINE1 其余三类上界", max(o3), tol=0.15, grp=1)
    transcribe(txt_en, r"restricting T cells to CD3D/PTPRC double-positive cells raised "
                       r"the fraction to ([\d.]+)%",
               "SERPINE1 T 细胞明确细胞比例（文本侧）",
               S("T 细胞", "SERPINE1", "明确细胞_表达比例"), tol=0.02)

    # 11.8 ML 三套 AUC（本轮新增的第三套 0.727 必须与另两套并列写出）
    ml = ml12.loc["GSE89408"]
    transcribe(txt_en, r"joint model AUC \*\*([\d.]+)\*\* by repeated internal 5-fold",
               "GSE89408 内部重复 CV AUC（权威值）", ml["面板CV_AUC.new"])
    transcribe(txt_en, r"the same cohort, median ([\d.]+), 95th percentile ([\d.]+); "
                       r"empirical P=([\d.]+)\)",
               "零分布中位", ml["零分布中位.new"])
    transcribe(txt_en, r"median [\d.]+, 95th percentile ([\d.]+); empirical P=([\d.]+)\)",
               "零分布 95 分位", ml["零分布95分位.new"], grp=1)
    transcribe(txt_en, r"95th percentile [\d.]+; empirical P=([\d.]+)\)",
               "定标经验 P", ml["经验P.new"], grp=1)
    # 单次 5 折的两值与「过滤前基因数」**不在** ml10（那是候选基因计数表），
    # 出自 B6 重跑记录；先与 ml10 对账候选基因数，再从 B6 取另三项。
    b6 = open(os.path.join(TAB, "B6_F11下游重跑_结果.md"), encoding="utf-8").read()
    m_b6 = re.search(r"候选基因 14,510（([\d,]+) 个中保留）\s*\|\s*([\d,]+)\s*\|\s*"
                     r"\*\*([\d,]+)\*\*", b6)
    if m_b6 is None:
        res.append((False, "[提取失败] B6 记录的候选基因行", "锚点未匹配", "有"))
        n_pre, n_old, n_new = 24894.0, 14510.0, 14130.0
    else:
        n_pre, n_old, n_new = (tonum(m_b6.group(1)), tonum(m_b6.group(2)),
                               tonum(m_b6.group(3)))
    chk("B6 与 ml10 的候选基因数一致（旧）", float(ml10.iloc[1, 1]), n_old, tol=0)
    chk("B6 与 ml10 的候选基因数一致（新）", float(ml10.iloc[1, 2]), n_new, tol=0)
    m_cv = re.search(r"\|\s*11 号[^|]*\|[^|]*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|", b6)
    if m_cv is None:
        res.append((False, "[提取失败] B6 记录的单次 5 折行", "锚点未匹配", "有"))
        cv_old, cv_new = 0.575, 0.727
    else:
        cv_old, cv_new = float(m_cv.group(1)), float(m_cv.group(2))
    transcribe(txt_en, r"gives \*\*([\d.]+)\*\* on the pre-correction expression matrix "
                       r"and \*\*([\d.]+)\*\* on the corrected one",
               "单次 5 折 AUC（旧）", cv_old, tol=0.001)
    transcribe(txt_en, r"gives \*\*[\d.]+\*\* on the pre-correction expression matrix "
                       r"and \*\*([\d.]+)\*\* on the corrected one",
               "单次 5 折 AUC（新）", cv_new, tol=0.001, grp=1)
    transcribe(txt_en, r"single-gene AUC ([\d.]+)–([\d.]+)",
               "单基因 AUC 下界", float(ml11["单基因AUC下限"]), tol=0.001)
    transcribe(txt_en, r"single-gene AUC [\d.]+–([\d.]+)",
               "单基因 AUC 上界", float(ml11["单基因AUC上限"]), tol=0.001, grp=1)
    transcribe(txt_en, r"\*\*([\d,]+)\*\* candidate genes, from ([\d,]+) retained",
               "CPM 规则候选基因数", n_new)
    transcribe(txt_en, r"\*\*[\d,]+\*\* candidate genes, from ([\d,]+) retained",
               "过滤前基因数", n_pre, grp=1)
    m = re.search(r"the cross-platform adaptive rule used to make all five cohorts "
                  r"comparable gives \*\*([\d,]+)\*\*", txt_en)
    res.append((m is not None and tonum(m.group(1)) == float(ml11["候选基因"]),
                "[EN] 转写 自适应规则候选基因数",
                m.group(1) if m else "锚点未匹配", ml11["候选基因"]))

    # 11.9 CN 侧平行提取（同一批数字，中文稿必须一致）
    transcribe(txt, r"成纤维细胞占比相差 ([\d.]+) 个百分点", "成纤维细胞占比极差",
               spread, tol=0.02, prefix="CN")
    transcribe(txt, r"随机效应下 P=([\d.]+) 与 ([\d.]+)", "复制层 LOXL1 的 RE-DL P",
               over.loc["LOXL1", "rep_re_P"], prefix="CN")
    transcribe(txt, r"随机效应下 P=[\d.]+ 与 ([\d.]+)", "复制层 LOX 的 RE-DL P",
               over.loc["LOX", "rep_re_P"], grp=1, prefix="CN")
    transcribe(txt, r"按 \|log2FC\| 排第 (\d+)，按 P 排第 (\d+)",
               "复制层 PLOD1 的 |log2FC| 名次", over.loc["PLOD1", "rank_abs_rep"],
               prefix="CN")
    transcribe(txt, r"按 \|log2FC\| 排第 \d+，按 P 排第 (\d+)",
               "复制层 PLOD1 的 P 名次", over.loc["PLOD1", "rank_P_rep"], grp=1,
               prefix="CN")
    # 与 EN 侧同理：原先 want=max(两队列), tol=0.002，因两队列仅差 0.0007 而恒真。
    transcribe(txt, r"全基因组中位可测下限为 ([\d.]+) log2FC",
               "82107 中位 MDE（全基因组）", diag.loc["GSE82107", "MDE_logFC_med"],
               tol=0.0005, prefix="CN")
    transcribe(txt, r"log2FC，两个 GPL96 队列为 ([\d.]+) 与 ([\d.]+)",
               "GSE55235 中位 MDE", diag.loc["GSE55235", "MDE_logFC_med"],
               tol=0.0005, prefix="CN")
    transcribe(txt, r"log2FC，两个 GPL96 队列为 [\d.]+ 与 ([\d.]+)",
               "GSE55457 中位 MDE", diag.loc["GSE55457", "MDE_logFC_med"],
               tol=0.0005, grp=1, prefix="CN")
    transcribe(txt, r"共同 d 尺度上为 ([\d.]+) 对 [\d.]+",
               "82107 中位 MDE（d 尺度）", diag.loc["GSE82107", "MDE_d_med"],
               tol=0.005, prefix="CN")
    transcribe(txt, r"共同 d 尺度上为 [\d.]+ 对 ([\d.]+)",
               "GPL96 中位 MDE（d 尺度）", diag.loc["GSE55235", "MDE_d_med"],
               tol=0.005, prefix="CN")
    transcribe(txt, r"该队列 57 个候选基因中仅 (\d+) 个达名义 P<0\.05，GSE55457 则为 (\d+) 个",
               "82107 名义 P<0.05 的候选数", n82107, prefix="CN")
    transcribe(txt, r"GSE55457 则为 (\d+) 个", "GSE55457 名义 P<0.05 的候选数",
               n55457, grp=1, prefix="CN")
    transcribe(txt, r"表达比例在成纤维细胞中最高（([\d.]+)%，其余四类细胞 1\.6%–6\.1%",
               "LOXL1 成纤维比例（文本侧）", P("LOXL1", "Fibroblast"), tol=0.20,
               prefix="CN")
    transcribe(txt, r"以内皮细胞最高（([\d.]+)%），但在其余三类中仍有 ([\d.]+)%–([\d.]+)%",
               "SERPINE1 内皮比例（文本侧）", P("SERPINE1", "Endothelial"), tol=0.90,
               prefix="CN")
    transcribe(txt, r"但在其余三类中仍有 ([\d.]+)%–[\d.]+%",
               "SERPINE1 其余三类下界", min(o3), tol=0.10, prefix="CN")
    transcribe(txt, r"但在其余三类中仍有 [\d.]+%–([\d.]+)%",
               "SERPINE1 其余三类上界", max(o3), tol=0.15, grp=1, prefix="CN")
    transcribe(txt, r"该比例升至 ([\d.]+)%", "SERPINE1 T 细胞明确细胞比例（文本侧）",
               S("T 细胞", "SERPINE1", "明确细胞_表达比例"), tol=0.02, prefix="CN")
    transcribe(txt, r"下为 \*\*([\d.]+)\*\*，对照\*\*随机基因集零分布\*\*",
               "GSE89408 内部重复 CV AUC（权威值）", ml["面板CV_AUC.new"], prefix="CN")
    transcribe(txt, r"中位 ([\d.]+)，95 百分位 ([\d.]+)；经验 P=([\d.]+)）",
               "零分布中位", ml["零分布中位.new"], prefix="CN")
    transcribe(txt, r"中位 [\d.]+，95 百分位 ([\d.]+)；经验 P=([\d.]+)）",
               "零分布 95 分位", ml["零分布95分位.new"], grp=1, prefix="CN")
    transcribe(txt, r"95 百分位 [\d.]+；经验 P=([\d.]+)）", "定标经验 P",
               ml["经验P.new"], grp=1, prefix="CN")
    transcribe(txt, r"校正前表达矩阵上为 \*\*([\d.]+)\*\*、在校正后为 \*\*([\d.]+)\*\*",
               "单次 5 折 AUC（旧）", cv_old, tol=0.001, prefix="CN")
    transcribe(txt, r"校正前表达矩阵上为 \*\*[\d.]+\*\*、在校正后为 \*\*([\d.]+)\*\*",
               "单次 5 折 AUC（新）", cv_new, tol=0.001, grp=1, prefix="CN")

    # ═══════════════════════ §12 路径与字面串变更的反向守卫 ═══════════════════════
    res.append((os.path.exists(MS) and os.path.exists(MS_EN),
                "§12 重定位后的两稿均存在", "存在", "存在"))
    res.append((os.path.exists(OLD_EN), "§12 旧稿仍在原处（未被覆盖）", "存在", "存在"))
    res.append(("重定位2026-09-23" in MS_EN and "重定位2026-09-23" in MS,
                "§12 回验对象确为重定位稿", "是", "是"))
    res.append(("Cross-Cohort Expression Profile and Cellular Localization"
                not in txt_en,
                "§12 旧稿标题未混入新稿（防误拿旧稿当新稿）", "无", "无"))
    res.append((not MS.endswith("修订2026-09-17.md"),
                "§12 MS 路径未退回旧稿", MS.split(os.sep)[-1], "非 09-17 稿"))
    res.append((abs(OVERRIDE - 0.007142857) < 1e-8,
                "§12 面板 Bonferroni 分母确为 7（非 6）", round(OVERRIDE, 8),
                round(0.05 / 7, 8)))

    # ---- 输出 ----
    npass = sum(1 for r in res if r[0])
    for ok, name, got, want in res:
        if not ok:
            print(f"[不通过] {name}: 实测={got} 稿件={want}")
    print(f"\n共 {len(res)} 项，通过 {npass} 项，不通过 {len(res) - npass} 项")
    print(f"（§12 登记：字面串 absent 增补 {len(ABSENT_ADDED)} 条、"
          f"CN present 增补 {len(PRESENT_ADDED_CN)} 条、EN present 增补 {len(PRESENT_ADDED_EN)} 条；"
          f"路径由 2026-09-17 稿改指 2026-09-23 重定位稿）")
    return 0 if npass == len(res) else 1


if __name__ == "__main__":
    sys.exit(main())
