#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V3 — 图注的结构描述回图现查（2026-09-29）

对应计划 F2 的第二句：「图注的结构描述（行数、栏数、所画队列数）回图现查」。

为什么单开一件：E3 已经把图 1（内嵌文本 vs 终稿图注）、图 3（行标签三方对账）、
图 4（内嵌数值 vs Table 3）、图 5（画出的细胞类别数 vs 图注计数）对过账；
**图 2 一处未覆盖**——它的图注里有一串结构声明（六个面板基因、八个交联基因、
五个无复制层值、点色/点面积的所指），而 E3 对图 2 只查了「图内无图号残留」。
本件把五幅图的**结构计数**统一成一组判据，逐条对撞图面文本层。

判据的形态：图注自称的数 N ↔ 图面实际画出的东西。凡「图注写 N、图面是 M」即失败。
基因分组不写死在脚本里，从冻结表 `02_结果表/57基因总览.csv` 的 `in_panel` /
`is_crosslink` 两列现取——写死会在换基因集时静默失效（见记忆：核心基因集一变
所有派生图都要改）。

图内文本用 `03_分析/scripts/pdftext.py` 的取词器（E3 同款）：matplotlib 会把一个
标签切成多个 Tj 操作数，自行解流后按字面量搜会对长串假报「不存在」。

用法：python3 V3_图注结构描述对图现查_2026-09-29.py ; echo $?
"""

import csv
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parent
SCRIPTS = OUT / "01_冻结结果与代码" / "03_分析" / "scripts"
FIGS = OUT / "03_图表与源数据"
MS = OUT / "02_共用主稿" / "English_manuscript_共用清洁稿_2026-09-28i.md"
OVERVIEW = OUT / "01_冻结结果与代码" / "02_结果表" / "57基因总览.csv"

sys.path.insert(0, str(SCRIPTS))
from pdftext import pdf_text as _house_pdf_text   # noqa: E402

fails = []
oks = [0]


def must(ok, msg, detail=None):
    if ok:
        oks[0] += 1
        print(f"  OK   {msg}")
    else:
        fails.append(msg)
        print(f"  FAIL {msg}" + (f"  | {detail}" if detail is not None else ""))
    return ok


# 交付命名的图号 → 词干。**必须按交付件读**：`01_冻结结果与代码/04_图表/Fig_*.pdf`
# 是绘图脚本的落点，交付目录里的是 `Figure{N}_{stem}.pdf`；读前者会漏掉 E3 组装时
# 才发生的重命名与替换（此前 E2 未运行却照常通过，正是这类「读了另一份」的产物）。
DELIVERED = {
    1: "study_design",
    2: "candidate_overview",
    3: "forest_logFC",
    4: "evidence_matrix",
    5: "scRNA_dotplot",
}


def fig_text(n):
    p = FIGS / f"Figure{n}_{DELIVERED[n]}.pdf"
    if not p.is_file():
        raise FileNotFoundError(f"交付图缺：{p}")
    recs, _ = _house_pdf_text(p)
    txt = " ".join(t for _, _, t in recs).replace("−", "-")
    # 压掉连续空白：图面是排版文字，`1  Candidate scope` 在源里是**两个**空格
    # （与稿件图注的单空格不同）。不压会让逐字比对假报「图面缺该标签」。
    return " ".join(txt.split())


def captions(md):
    """按**段落**取图注：段落首为 `**Figure N:` 的那一段。

    不用 `re.search(r'Figure 5[:.](.+?)(?=\\n\\n)')`：正文里的「(Figure 5).」
    也满足该模式，会先命中正文而取到一段错误的「图注」。段落切分没有这个歧义。
    """
    out = {}
    for para in md.split("\n\n"):
        m = re.match(r"^\*\*Figure (\d)[:.]", para.strip())
        if m:
            out[int(m.group(1))] = " ".join(para.split())
    return out


WORDS = {"three": 3, "four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8}


def claimed_count(caption, pat):
    """图注自称的个数 → 数字。取不到返回 (None, None)。

    **这是判据本体**：正查与负面对照都必须走它。否则对照只是在改过的测试数据上
    重抄一遍比较式，证不了判据本身会失败（本户已吃过一次空转对照的亏）。
    全词表也是必须的：缺 "five" 时 `.get("five", 0) == 5` 恒假，会把真结论误报成
    失败——本次实得 E3 失败即此，故词表与判据一并放在这里。
    """
    m = re.search(pat, caption, re.I)
    if not m:
        return None, None
    return WORDS.get(m.group(1).lower()), m.group(1)


def count_matches(caption, pat, actual):
    """图注自称的个数是否等于实际个数。"""
    n, _ = claimed_count(caption, pat)
    return n == actual


def names_match(caption, named_pat, expected):
    """图注点名的集合是否与冻结表实际的集合**同一组**（不是包含关系）。"""
    # re.I 不可省：图注写 "Five candidates (...)"，模式是全小写；漏了就对不上，
    # 且调用处若用 `if m:` 接住会**静默跳过**整条判据（判据空转的典型形态）。
    m = re.search(named_pat, caption, re.I)
    if not m:
        return False
    return sorted(x.strip() for x in m.group(1).split(",")) == sorted(expected)


# 图注侧的模式常量：正查与负面对照**共用同一份**。两处各写一遍会漂移，
# 而漂移后对照仍会「通过」。
PAT_LAYERS = r"\b(three|four|five)\b evidence layers are reported in order"
PAT_PANEL_GENES = r"the (six|seven|eight) panel genes"
PAT_XLINK_GENES = r"the (six|seven|eight) cross-linking genes"
PAT_NOREP = r"(five|four|six) candidates \([^)]+\) have no replication-layer value"
PAT_NOREP_NAMES = r"(?:five|four|six) candidates \(([^)]+)\) have no replication-layer value"
PAT_COHORTS = r"all (six|five|seven) cohorts"
PAT_POPULATIONS = r"the (four|five|six) annotated OA synovial populations"


def main():
    for p in (MS, OVERVIEW):
        if not p.is_file():
            print(f"FATAL 缺 {p}")
            return 2
    md = MS.read_text(encoding="utf-8")
    caps = captions(md)
    must(sorted(caps) == [1, 2, 3, 4, 5], "五条图注按段落取全", f"实得 {sorted(caps)}")

    rows = list(csv.DictReader(open(OVERVIEW, encoding="utf-8-sig")))
    PANEL = [r["gene"] for r in rows if r["in_panel"].strip().upper() == "TRUE"]
    XLINK = [r["gene"] for r in rows if r["is_crosslink"].strip().upper() == "TRUE"]
    NOREP = [r["gene"] for r in rows if r["rep_k"].strip() in ("NA", "")]
    print(f"  冻结表：面板 {len(PANEL)} 个、交联 {len(XLINK)} 个、复制层不可汇总 {len(NOREP)} 个")
    must(len(PANEL) == 6 and len(XLINK) == 8, "冻结表给出 6 面板 / 8 交联",
         f"{len(PANEL)}/{len(XLINK)}")

    print("\n[A] 图 1：图注自称的证据层数 vs 图面编号层标签")
    t1 = fig_text(1)
    LAYERS = ["1 Candidate scope", "2 Repeatability", "3 Cross-linker family", "4 Cellular source"]
    miss = [s for s in LAYERS if s not in t1]
    must(not miss, f"图 1 图面画出四个编号层标签", f"缺 {miss}")
    n_cap, w_cap = claimed_count(caps[1], PAT_LAYERS)
    must(w_cap is not None, "A2 图注含「N evidence layers are reported in order」句")
    must(count_matches(caps[1], PAT_LAYERS, len(LAYERS)),
         f"A3 图注所报层数（{n_cap}）等于图面画出的层数（{len(LAYERS)}）",
         f"图注 {w_cap}，图面 {len(LAYERS)}")
    n_enum = len(re.findall(r"\((\d)\) ", caps[1]))
    must(n_enum == len(LAYERS), f"A4 图注逐层列举了 {len(LAYERS)} 项", f"实得 {n_enum}")

    print("\n[B] 图 2：图注自称的基因数 vs 图面标签")
    t2 = fig_text(2)
    miss_p = [g for g in PANEL if g not in t2]
    miss_x = [g for g in XLINK if g not in t2]
    must(not miss_p, f"B1 图面画出全部 {len(PANEL)} 个面板基因标签", f"缺 {miss_p}")
    must(not miss_x, f"B2 图面画出全部 {len(XLINK)} 个交联基因标签", f"缺 {miss_x}")
    union = set(PANEL) | set(XLINK)
    must(len(union) == 13, "B3 两个集合的并集为 13（LOXL1 同属两侧）", f"实得 {len(union)}")
    np_, w_p = claimed_count(caps[2], PAT_PANEL_GENES)
    nx_, w_x = claimed_count(caps[2], PAT_XLINK_GENES)
    must(w_p is not None and w_x is not None, "B4 图注含「N panel genes and M cross-linking genes」句")
    must(count_matches(caps[2], PAT_PANEL_GENES, len(PANEL)),
         f"B5a 图注所报面板基因数（{w_p}）等于冻结表的 {len(PANEL)}")
    must(count_matches(caps[2], PAT_XLINK_GENES, len(XLINK)),
         f"B5b 图注所报交联基因数（{w_x}）等于冻结表的 {len(XLINK)}")
    # 图注点名的「无复制层值」几个，必须与冻结表的 rep_k 缺失组**同一组**
    ng_, w_n = claimed_count(caps[2], PAT_NOREP)
    must(ng_ is not None, "B6 图注点名了无复制层值的候选")
    got_names = re.search(PAT_NOREP_NAMES, caps[2], re.I)
    must(got_names is not None, "B6b 图注的名单可解析（取不到即失败，不得静默跳过）")
    named = [x.strip() for x in got_names.group(1).split(",")] if got_names else []
    must(names_match(caps[2], PAT_NOREP_NAMES, NOREP),
         f"B7 图注点名的 {len(named)} 个与冻结表 rep_k 缺失的 {len(NOREP)} 个是同一组",
         f"图注 {sorted(named)} vs 表 {sorted(NOREP)}")
    must(ng_ == len(named),
         f"B8 图注所报个数（{w_n}）与点名个数（{len(named)}）一致")
    must("5 have no replication-layer value" in t2 or "not combinable" in t2,
         "B9 图面同句标注了不可汇总", f"图面未出现该标注")

    print("\n[C] 图 3：图注自称的队列数 vs 图面队列号")
    t3 = fig_text(3)
    DISC = ["GSE89408", "GSE206848", "GSE1919"]
    REP = ["GSE55235", "GSE55457", "GSE82107"]
    miss3 = [g for g in DISC + REP if g not in t3]
    must(not miss3, "C1 图面画出全部六个队列号", f"缺 {miss3}")
    must(all(g in caps[3] for g in DISC + REP), "C2 图注也点名了同样六个队列")
    n3, w3 = claimed_count(caps[3], PAT_COHORTS)
    must(w3 is not None, "C3 图注含「all N cohorts」句")
    must(count_matches(caps[3], PAT_COHORTS, len(DISC + REP)),
         f"C4 图注所报队列数（{w3}）等于图面画出的（{len(DISC + REP)}）")

    print("\n[D] 图 4：图注的列清单 vs 图面列头（栏数与次序）")
    t4 = fig_text(4)
    # 图面六组列头，其子标签按图面印刷次序。组名里的换行是绘图时的条带折行，
    # 取词器把它并成一个空格（见 fig_text 的空白压缩），故此处用空格写法。
    GROUPS = {
        "Discovery primary (REML+HK)": ["log2FC", "P", "I2 (%)"],
        "Discovery sensitivity (IVW-FE)": ["log2FC", "P"],
        "Replication primary (REML+HK)": ["log2FC", "P", "I2 (%)"],
        "Replication sensitivity (pooled)": ["log2FC", "P"],
        "Rank (REML+HK, both layers)": ["|log2FC| disc", "P disc", "|log2FC| rep", "P rep"],
        "Deletion (discovery, IVW-FE)": ["GSE1919", "GSE206848", "GSE89408"],
    }
    for g, subs in GROUPS.items():
        must(g in t4, f"D1 图面有列组「{g}」")
        missd = [s for s in subs if s not in t4]
        must(not missd, f"D2 列组「{g}」的 {len(subs)} 个子标签齐备", f"缺 {missd}")
    n_cols = sum(len(v) for v in GROUPS.values())
    must(n_cols == 17, f"D3 图面共 {n_cols} 栏（3+2+3+2+4+3）", n_cols)
    # 图注侧：图注按「每层一主一敏」成句（主估计量 3 栏 ＋ 敏感性 2 栏），不逐栏罗列。
    # 故不数栏，改为逐项断言图注确实写出了每层的可辨别特征——**含估计量归属**：
    # 敏感性那两半句写明了谁配 fixed-effect、谁配 pooled，与图面条带一一对应；
    # 若图改了而图注留旧，D4 必失败（这正是本轮之前 15 栏版与图注脱节的那类事故）。
    NEED = ["discovery-layer log2FC, P and I² under the declared primary estimator",
            "the discovery-layer inverse-variance-weighted fixed-effect log2FC and P"
            " as its sensitivity counterpart",
            "the same five columns for the replication layer with the pooled log2FC and P"
            " as its sensitivity counterpart",
            "the rank of each gene under each layer's |log2FC| and P criteria",
            "the outcome of the three single-cohort deletions"]
    missn = [s for s in NEED if s not in caps[4]]
    must(not missn, f"D4 图注逐项列出了 {len(NEED)} 个列组成分", f"缺 {missn}")
    # D4b：图注与图面的**层—敏感性估计量配对**必须一致（图中两层主估计同名 REML+HK，
    # 敏感性各不同），故逐层对本：图面条带里的简称 ↔ 图注里的全称。
    SENS = [("Discovery", "inverse-variance-weighted fixed-effect", "IVW-FE"),
            ("Replication", "pooled", "pooled")]
    for lay, long_, short_ in SENS:
        strip = f"{lay} sensitivity ({short_})"
        must(long_ in caps[4] and strip in t4,
             f"D4b 图注的 {lay} 层敏感性全称「{long_}」与图面条带「{strip}」配对")
    must("three single-cohort deletions" in caps[4] and len(GROUPS["Deletion (discovery, IVW-FE)"]) == 3,
         "D5 图注所称「three single-cohort deletions」与图面三栏一致")
    must("Rows are the eight cross-linking genes" in caps[4] and len(XLINK) == 8,
         "D6 图注所称「eight cross-linking genes」行数与冻结表一致")
    order = ["LOXL1", "LOX", "PLOD1", "PLOD2", "LOXL2", "TGM2", "LOXL4", "LOXL3"]
    must(set(order) == set(XLINK), "D7 图面行序里的八个基因与冻结表 is_crosslink 同一集合",
         f"图面序 {sorted(order)} vs 表 {sorted(XLINK)}")

    print("\n[E] 图 5：图注自称的类别数 vs 图面细胞类型（复核 E3 [E]）")
    t5 = fig_text(5)
    CELLS = ["Macrophage", "Fibroblast", "Endothelial", "T cell", "Mast"]
    miss5 = [c for c in CELLS if c not in t5]
    must(not miss5, f"E1 图面画出全部 {len(CELLS)} 类细胞", f"缺 {miss5}")
    n5, w5 = claimed_count(caps[5], PAT_POPULATIONS)
    must(w5 is not None, "E2 图注含「the N annotated OA synovial populations」句")
    must(count_matches(caps[5], PAT_POPULATIONS, len(CELLS)),
         f"E3 图注所报类别数（{w5}）等于图面画出的（{len(CELLS)}）")
    # 图 5 的面板数：合并细胞 + 三供者
    panels = [p for p in ("Pooled cells (all donors)", "SY_113", "SY_116", "SY_118") if p in t5]
    must(len(panels) == 4, "E4 图面四个面板（合并细胞 + 三供者）齐备", f"实得 {panels}")
    must("per donor" in caps[5] and "pooled" in caps[5],
         "E5 图注说明了「逐供者」与「合并细胞」两种口径")
    must("mean log1p(CP10K)" in caps[5] and "Mean log1p(CP10K)" in t5,
         "E6 图注与图面的均值口径同为 mean log1p(CP10K)")

    print("\n[F] 负面对照：判据必须能失败")
    # F1 图注侧：把图 2 图注所报的面板基因数从 six 改成 eight。先自证改写生效，
    # 再**调用判据本体** count_matches —— 只重抄比较式证不了判据会失败。
    neg = caps[2].replace("the six panel genes", "the eight panel genes", 1)
    must("the eight panel genes" in neg and "the six panel genes" not in neg,
         "F1 前提：改写确已生效（six→eight，旧串归零）")
    must(not count_matches(neg, PAT_PANEL_GENES, len(PANEL)),
         f"F1 负面对照：图注数改成 8 后 count_matches 判据翻转为 False",
         f"仍返回 {count_matches(neg, PAT_PANEL_GENES, len(PANEL))}")
    # F2 图面侧：把图 4 文本层里的一个列组名改掉，D1 必须翻转
    t4neg = t4.replace("Deletion (discovery, IVW-FE)", "Deletion", 1)
    must("Deletion (discovery, IVW-FE)" in t4 and "Deletion (discovery, IVW-FE)" not in t4neg,
         "F2 负面对照：图面列组名被改后 D1 判据翻转为失败（改前为真、改后为假）")
    # F3 命名集合侧：把「无复制层值」的名单里去掉一个，判据本体 names_match 必须翻转
    got3 = re.search(PAT_NOREP_NAMES, caps[2], re.I)
    must(got3 is not None, "F3 前提：原图注的名单可解析")
    named = [x.strip() for x in got3.group(1).split(",")] if got3 else []
    neg3 = caps[2].replace(", ".join(named), ", ".join(named[:-1]), 1)
    must(neg3 != caps[2] and names_match(caps[2], PAT_NOREP_NAMES, NOREP),
         "F3 前提：负例确已改写（少一个名），且原图注在判据下为真")
    must(not names_match(neg3, PAT_NOREP_NAMES, NOREP),
         "F3 负面对照：少点一个名后 names_match 判据翻转为 False",
         f"仍返回 {names_match(neg3, PAT_NOREP_NAMES, NOREP)}")
    # F4 前提：上述三条对照确实作用在**真实产物**上，不是空串
    must(len(caps[2]) > 200 and len(t4) > 500, "F4 对照前提：图注与图面文本层均非空",
         f"图注 {len(caps[2])} 字符 / 图面 {len(t4)} 字符")

    print(f"\n结果：{'全部通过' if not fails else '存在失败'}")
    print(f"通过 {oks[0]} 项，失败 {len(fails)} 项")
    for f in fails:
        print(f"  - {f}")
    if not fails:
        print("\n判读：五幅图的结构计数（层数、基因标签数、队列数、栏数、细胞类别数）"
              "全部与图注自称的一致；其中图 2 的「六个面板 / 八个交联 / 五个无复制层值」"
              "三项为 E3 未覆盖者，本件补齐。")
    else:
        print(f"\n判读：有 {len(fails)} 条结构声明与图面不符，见上。判读句不结论。")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
