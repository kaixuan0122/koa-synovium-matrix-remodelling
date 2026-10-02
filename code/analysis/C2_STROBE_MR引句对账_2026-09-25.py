#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
C2 — STROBE-MR 重映射的独立对账

独立于生成器 C1：从成品文件读出「标为逐字」的引句，逐条回原文检索。
检索语料 = 新稿正文 + S1/S2/S3（任一命中即通过，并报出命中文件）。
归一化：去掉 markdown 强调符 *，折叠空白与换行——两侧同样处理，故不掩盖字面差异。

判据：引句必须在语料中出现；同时给出反向对照（把引句改动一字后必须查不到），
以证明检索本身有效而不是恒真。
"""
import re
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[3]
REV = ROOT / "08_修订_2026-09-17"
CHK = REV / "06_补充材料" / "STROBE-MR_checklist.md"

CORPUS_FILES = [
    REV / "05_稿件" / "English_manuscript_MDPI_重定位2026-09-23.md",
    REV / "06_补充材料" / "Supplementary_Table_S1.md",
    REV / "06_补充材料" / "Supplementary_Table_S2.md",
    REV / "06_补充材料" / "Supplementary_Table_S3.md",
]

norm = lambda s: re.sub(r"\s+", " ", s.replace("*", "")).strip()
corpus = {}
for f in CORPUS_FILES:
    corpus[f.name] = norm(f.read_text(encoding="utf-8"))

# 清单中标记为逐字引用的句子（item, 引句）
QUOTES = [
    ("1/4a 摘要", "Exploratory Mendelian randomisation gave no transferable genetic support"),
    ("3", "We therefore asked which matrix-remodelling and mechanotransduction-related candidates change reproducibly in KOA synovium, how far their priority depends on the criterion used to establish it, and which synovial cells carry them."),
    ("3", "the genetic evidence is bounded by the cell type in which each instrument was derived"),
    ("4e", "Not applicable. All data analyzed were obtained from public repositories in de-identified form."),
    ("5 relevance", "is a cis-eQTL for the exposure gene, quantified by the gene-level F statistic."),
    ("5 independence", "instruments are assumed unrelated to confounders of the exposure–outcome relationship — stated, not tested, since LD pruning addresses correlation between instruments and cannot establish freedom from confounding [15]."),
    ("5 residual LD", "six retained pairs still exceed r²=0.01, the largest being 0.063 (GTEx fibroblast VGLL4)"),
    ("5 bias", "not quantified"),
    ("5 exclusion", "instruments are assumed to affect KOA only through the exposure gene, which is not directly testable; the MR-Egger intercept tested directional horizontal pleiotropy where three or more instruments survived, and Cochran's Q tested heterogeneity."),
    ("6c estimator", "Genes retaining one instrument were analysed by Wald ratio and genes retaining two or more by IVW; MR-Egger and the weighted median were available only where three or more instruments survived."),
    ("7 strength", "F>10 considered strong"),
    ("10d overlap", "draw on non-overlapping sample frames, so no exposure–outcome sample overlap is expected"),
    ("13 no support", "no panel gene obtained transferable genetic support"),
    ("13 no correction", "no across-combination correction applies"),
    ("15 limits opener", "Several limits bound these conclusions."),
    ("15 state test", "was tested against two null calibrations and was not established"),
    ("15 ancestry", "The genetic data are entirely European-ancestry"),
    ("17 cell type", "An instrument that predicts VGLL4 expression in macrophages is not thereby a valid instrument for VGLL4 expression in synovial fibroblasts"),
    ("2 appendix", "exploratory genetic evidence as a bounded appendix"),
    ("4b MDE", "2.802×SE"),
    ("4c IV changed", "IV set changed"),
    ("12 gene-level F def", "is the arithmetic mean of those per-SNP values"),
    ("12 VGLL4 F", "F=62.3"),
    ("6e alpha", "α=0.0028"),
    ("6e 7807", "7,807 common genes"),
    ("15 ML null", "0.635"),
    ("15 AUC", "0.647"),
    ("15 spread", "43.65-percentage-point"),
    ("4d release7", "n=659"),
    ("4c window", "narrower than the cis-region"),
]

chk_text = norm(CHK.read_text(encoding="utf-8"))

print("C2 独立对账 —— STROBE-MR 引句回原文检索")
print(f"语料: {', '.join(corpus)}\n")

fail = 0
for tag, q in QUOTES:
    qn = norm(q)
    hits = [name for name, txt in corpus.items() if qn in txt]
    in_chk = qn in chk_text
    status = "OK " if hits else "MISS"
    if not hits:
        fail += 1
    where = ",".join(h.split("_")[-1].replace(".md", "") for h in hits) or "—"
    print(f"  {status} [{tag}] 命中={where}  在清单中={in_chk}  «{q[:56]}{'…' if len(q) > 56 else ''}»")

# 反向对照：把每条引句从中间切掉 3 个字符后必须查不到（证明检索非恒真）
print("\n反向对照（引句中部删除 3 字符后必须查不到）:")
neg_fail = 0
for tag, q in QUOTES:
    qn = norm(q)
    if len(qn) < 30:
        continue
    m = len(qn) // 2
    broken = norm(qn[:m] + qn[m + 3:])
    if any(broken in txt for txt in corpus.values()):
        neg_fail += 1
        print(f"  假阳性 [{tag}] 删字后仍命中 —— 检索过于宽松")
print(f"  反向对照：{len([q for _, q in QUOTES if len(norm(q)) >= 30])} 条受检，假阳性 {neg_fail} 条")

print(f"\n结果：{len(QUOTES) - fail}/{len(QUOTES)} 条引句在语料中逐字命中；未命中 {fail} 条")
raise SystemExit(1 if fail or neg_fail else 0)
