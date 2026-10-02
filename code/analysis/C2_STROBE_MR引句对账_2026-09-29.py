#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
C2 — STROBE-MR 引句对账（2026-09-29 移植版）

与 `C2_STROBE_MR引句对账_2026-09-25.py` 的差别：
  1. 旧版 `ROOT = parents[3]` 之后接的是旧基线目录 `08_修订_2026-09-17/`，其语料（09-23 稿、
     旧 S1–S3）在本树不存在，故不可直接复跑；本版语料改为当前树的主稿（09-28d）与 S1–S5。
  2. 增一条**生成式**判据：清单里所有「…」引句由清单自身现抽（不再只靠手抄的 QUOTES 表），
     按省略号切段、逐段回语料检索。清单改一次，本对账自动覆盖新句——手写清单永远会滞后。
  3. 保留旧版手抄清单中的短引句与**数值锚**（F=62.3、α=0.0028、0.647 等），这些短串
     生成式抽取捞不到（长度 < 阈值），但正是最容易被改写悄悄漂移的部分。

判据：
  · 每条引句必须在语料中逐字命中（归一化：去 markdown 强调符 *、折叠空白）；
  · 生成式抽取的结果必须非空且含长句（防抽取器自身退化后「零条待检」式空转）；
  · 反向对照：每条引句中部删 3 字符后必须查不到。全部对照「通过」是空转的签名。

首轮（2026-09-29）本脚本查出清单侧三处、已由 F6 的 D24/D24b/D25 修掉，三处都不是语义分歧，
而是**清单指向的对象已改名或退场**：item 7 引号内的句子在主稿里已不存在（`considered`→`taken`）、
item 12 的「gene-level F」这称法已随 §2.6 退场、item 4c 点名的列在 S3 里叫 `Set changed`。
本版把这三处固化为判据：新措辞须命中、旧措辞须查不到（下接 [C] 的反向项）。

语料 = 主稿 09-28d + Supplementary_Table_S1..S5。
"""
import re
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[3]
CHK = ROOT / "04_补充材料" / "STROBE-MR_checklist.md"
MS = ROOT / "02_共用主稿" / "English_manuscript_共用清洁稿_2026-09-28d.md"
SUP = ROOT / "04_补充材料"
CORPUS_FILES = [MS] + [SUP / f"Supplementary_Table_S{i}.md" for i in "12345"]

norm = lambda s: re.sub(r"\s+", " ", s.replace("*", "")).strip()


def seg_hits(seg, corpus):
    """省略号切出的单段：归一化后回语料检索，返回命中文件名列表。"""
    return [name for name, txt in corpus.items() if seg in txt]


fails = []


def must(ok, msg, detail=""):
    print(f"  {'OK  ' if ok else 'FAIL'} {msg}" + (f"   {detail}" if detail else ""))
    if not ok:
        fails.append(msg)


for f in CORPUS_FILES:
    if not f.is_file():
        raise SystemExit(f"语料缺失：{f}")
corpus = {f.name: norm(f.read_text(encoding="utf-8")) for f in CORPUS_FILES}
chk_text = norm(CHK.read_text(encoding="utf-8"))

print("C2 独立对账 —— STROBE-MR 引句回语料检索（2026-09-29 移植版）")
print(f"清单: {CHK.name}（{len(chk_text)} 字符归一化后）")
print(f"语料: {', '.join(corpus)}\n")

# ---------------------------------------------------------------- 生成式抽取
print("[A] 生成式抽取：清单自身的「…」引句（省略号切段，逐段检索）")
raw = re.findall(r'"([^"]{1,600})"', CHK.read_text(encoding="utf-8"))
AUTO, SKIPPED = [], []
for q in raw:
    segs = [norm(s) for s in re.split(r"…", q)]
    segs = [s for s in segs if len(s) >= 15]
    if segs:
        AUTO.append((q, segs))
    else:
        SKIPPED.append(q)

must(len(AUTO) >= 15, "抽取到足够多的可检引句（防抽取器退化后空转）", f"{len(AUTO)} 条")
must(any(len(s) >= 100 for _, segs in AUTO for s in segs),
     "抽取结果中含长句（≥100 字符），确认抽的不是标题碎片",
     f"最长 {max((len(s) for _, segs in AUTO for s in segs), default=0)}")
must(len(SKIPPED) <= 10, "被跳过的短引串数量有限（跳过的都是 <15 字符的标题/序数）",
     f"{len(SKIPPED)} 条：" + " | ".join(s[:24] for s in SKIPPED))

auto_bad = []
for q, segs in AUTO:
    for s in segs:
        if not seg_hits(s, corpus):
            auto_bad.append((q, s))
must(not auto_bad, "清单里的每条引句（含省略号切段）都在语料中逐字命中",
     f"{len(AUTO)} 条引句 / {sum(len(s) for _, s in AUTO)} 段，未命中 {len(auto_bad)} 段")
for q, s in auto_bad:
    print(f"       MISS «{s[:80]}…»  （出自清单引句 «{q[:60]}…»）")

# --------------------------------------------------------- 手抄锚（短串/数值）
print("\n[B] 手抄锚：短引句与数值（生成式抽取的长度阈值捞不到）")
QUOTES = [
    ("1/4a 摘要", "Exploratory Mendelian randomisation gave no transferable genetic support"),
    ("3 问题句", "We therefore asked which matrix-remodelling and mechanotransduction-related candidates change reproducibly in OA synovium"),
    ("3 边界句", "the genetic evidence is bounded by the cell type in which each instrument was derived"),
    ("2 附录定位", "exploratory genetic evidence as a bounded appendix"),
    ("4c window", "narrower than the cis-region"),
    ("4b MDE", "2.802"),
    ("4d release7", "n=659"),
    ("5 strength", "F>10"),
    ("5 相关性", "six retained pairs still exceed r²=0.01, the largest being 0.063 (GTEx fibroblast VGLL4)"),
    ("5 无关性", "stated, not tested"),
    ("5 偏倚", "not quantified"),
    ("6e alpha", "α=0.0028"),
    ("6e 7807", "7,807 common genes"),
    ("12 VGLL4 F", "F=62.3"),
    ("13 无支持", "no panel gene obtained transferable genetic support"),
    ("6e 跨组合校正（现措辞）", "no across-combination correction is applied"),
    ("4c 列名", "Set changed"),
    ("7 强度引句（现措辞）", "F>10 taken as strong"),
    ("15 开头", "Several limits bound these conclusions."),
    ("15 零标定", "was tested against two null calibrations and was not established"),
    ("15 祖源", "The genetic data are entirely European-ancestry"),
    ("15 AUC", "0.647"),
    ("15 ML null", "0.635"),
    ("15 spread", "43.65-percentage-point"),
    ("17 细胞类型", "An instrument that predicts VGLL4 expression in macrophages is not thereby a valid instrument for VGLL4 expression in synovial fibroblasts"),
]
manual_bad = []
for tag, q in QUOTES:
    qn = norm(q)
    hit = seg_hits(qn, corpus)
    ok = bool(hit)
    if not ok:
        manual_bad.append((tag, q))
    where = ",".join(h.split("Supplementary_Table_")[-1].replace(".md", "") for h in hit) or "—"
    print(f"  {'OK  ' if ok else 'MISS'} [{tag}] 命中={where}  «{q[:58]}{'…' if len(q) > 58 else ''}»")
must(not manual_bad, "手抄锚全部命中", f"{len(QUOTES)} 条，未命中 {len(manual_bad)} 条")

# ------------------------------------------------------------------ 反向对照
print("\n[C] 反向对照：每条引句中部删 3 字符后必须查不到")
NEG = [(f"auto{i}", s) for i, (_, segs) in enumerate(AUTO) for s in segs if len(s) >= 30]
NEG += [(f"manual:{t}", norm(q)) for t, q in QUOTES if len(norm(q)) >= 30]
neg_fp = 0
for tag, s in NEG:
    m = len(s) // 2
    broken = s[:m] + s[m + 3:]
    if seg_hits(broken, corpus):
        neg_fp += 1
        print(f"      假阳性 [{tag}] 删字后仍命中 —— 检索过于宽松  «{s[:60]}…»")
must(len(NEG) >= 20, "反向对照覆盖足够条目", f"{len(NEG)} 条受检")
must(neg_fp == 0, "反向对照：删 3 字符后无一条仍命中（证明检索非恒真）", f"假阳性 {neg_fp} 条")

# ------------------------------------------------- 退役措辞（须查不到，不是须命中）
print("\n[D] 退役措辞：本轮修掉的旧称法不得在清单或语料里复活")
RETIRED = [
    ("旧强度引句", "considered strong"),
    ("旧跨组合校正句", "no across-combination correction applies"),
    ("旧列名", "IV set changed"),
]
for tag, s in RETIRED:
    n_chk, n_cor = chk_text.count(s), sum(txt.count(s) for txt in corpus.values())
    must(n_chk == 0 and n_cor == 0, f"退役措辞零命中：{tag} {s!r}",
         f"清单 {n_chk} / 语料 {n_cor}")
# 「gene-level F」不能一刀切禁掉：清单 item 10d 是历史记录，记的正是上一版 Table 3 的列名
# （旧稿 09-17 版 Table 3 表头逐字为 Gene-level F）。故按行归属判：清单里每一处命中都必须
# 落在含「previous version」的句子里；语料侧则一处都不许有。
_gl_chk = [l for l in CHK.read_text(encoding="utf-8").split("\n") if "gene-level F" in l]
must(all("previous version" in l for l in _gl_chk),
     "「gene-level F」在清单里只出现在「上一版」的历史记录句中（逐行归属，非全篇禁词）",
     f"清单命中 {len(_gl_chk)} 行")
_n_cor_gl = sum(txt.count("gene-level F") for txt in corpus.values())
must(_n_cor_gl == 0, "语料里「gene-level F」零命中（现行产物已不用该称法）", f"实得 {_n_cor_gl}")
# 负面对照：把一条退役措辞注回语料，判据必须翻转为失败。
# 注意注入对象是**语料**（主稿），不是清单——退役句的现措辞住在主稿里；注到清单上会因
# 「清单本来就没有这句」而假通过（首版正是这么写的，前提断言当场失败）。
_GL_OLD = "no across-combination correction applies"
_GL_NEW = "no across-combination correction is applied"
_probe_ms = corpus[MS.name].replace(_GL_NEW, _GL_OLD, 1)
must(_probe_ms != corpus[MS.name] and _probe_ms.count(_GL_OLD) == 1,
     "负面对照前提：注入生效（主稿里确有可替换的现措辞）",
     f"主稿现措辞 {corpus[MS.name].count(_GL_NEW)} 处")
_n_old_cor = sum(txt.count(_GL_OLD) for txt in corpus.values())
must(_n_old_cor == 0 and _probe_ms.count(_GL_OLD) == 1,
     "负面对照：同一判据在真语料上通过、在注入旧措辞后失败（证明 [D] 能失败）")

print(f"\n结果：生成式 {len(AUTO)} 条 + 手抄锚 {len(QUOTES)} 条逐字命中；"
      f"未命中 {len(auto_bad) + len(manual_bad)} 条；反向对照假阳性 {neg_fp} 条；"
      f"退役措辞 {len(RETIRED)} 条零命中")
print("结果：全部通过" if not fails else f"结果：失败 {len(fails)} 项")
raise SystemExit(0 if not fails else 1)
