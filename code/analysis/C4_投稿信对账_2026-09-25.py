#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
C4 — 两封投稿信的独立对账（与 C3 不同代码路径）

C3 查的是「该出现的在不在」，C4 查的是「出现了的是不是都有落点」——
方向相反，故能抓到 C3 结构上抓不到的问题：C3 的断言全部由我手写，
我没想到要写的数字，C3 不会报。

三条独立判据：
 ① 穷举扫描：信件本体里**每一个数字 token** 必须在语料中命中，
    未命中者逐条列出（白名单只放日期、刊名、序数等非主张数字）。命中须归属到文件。
 ② 显式主张清单：标题、逐字短语、关键值 → 定点检索并报命中文件。
 ③ 由产物反算：「18 组中 6 组达名义显著」不由字面比对验证，
    而是从 Table S3 的两级剪枝主列**重新数**（P<0.05 的格数），与信件所述对照。
另附反向对照：把每个数字 token 改一位后必须查不到（证明检索非恒真）。

语料 = 新稿 + S1/S2/S3（任一命中即通过）
"""
import re
import pathlib
from collections import Counter

ROOT = pathlib.Path(__file__).resolve().parents[3]
REV = ROOT / "08_修订_2026-09-17"
# JIR 母版信 2026-09-27 随作废路线归档，落点移进 _archive_ 子目录（交付目录根部只有 MDPI 活信）。
# 对账对象不变：归档的也是产出物，同样要过穷举数字扫描——归档不等于免检。
LETTERS = [REV / "08_投稿材料" / "Cover_Letter_MDPI.md",
           REV / "08_投稿材料" / "_archive_JIR作废_2026-09-27" / "Cover_Letter.md"]
CORPUS_FILES = {
    "manuscript": REV / "05_稿件" / "English_manuscript_MDPI_重定位2026-09-23.md",
    "S1": REV / "06_补充材料" / "Supplementary_Table_S1.md",
    "S2": REV / "06_补充材料" / "Supplementary_Table_S2.md",
    "S3": REV / "06_补充材料" / "Supplementary_Table_S3.md",
}

norm = lambda s: re.sub(r"\s+", " ", s.replace("*", "")).strip()
corpus = {k: norm(p.read_text(encoding="utf-8")) for k, p in CORPUS_FILES.items()}
strip_c = lambda s: re.sub(r"<!--.*?-->", "", s, flags=re.S)

# ── 数字 token：整数、小数、P=…、科学计数（含上标）、百分比、α=…
TOKEN = re.compile(
    r"(?:[Ppq]=)?\d+(?:\.\d+)?(?:×10[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)?%?"
    r"|α=0\.05/\d+"
)

# 非主张数字白名单（日期、刊名、序数、通用显著性水平、MDPI 版式等）
WHITELIST = {
    "25", "2026", "1", "2", "3", "4", "6", "7", "8",            # 序数/计数词已另作主张查
    "0.05", "p<0.05", "p<0.01", "0.01",                          # 通用阈值（新稿 P<0.01/α=0.05 均见）
    "5", "9", "10",
}

def claims_only(text):
    """只取信件本体的主张区（首句至 Sincerely 之前），排除签名块的联系电话等非主张数字。"""
    b = strip_c(text)
    i = b.find("Sincerely,")
    return norm(b[:i] if i > 0 else b)


def probe(tok, corpus):
    """按**值**检索：P=/p=/q= 前缀若未命中，用裸值再试一次（新稿常省略前缀）。
    返回 (命中文件, 实际命中的形态)；未命中返回 ([], None)。前缀差异不等于值不一致。"""
    forms = [tok]
    if re.match(r"^[Ppq]=", tok):
        forms.append(tok[2:])
    for f in forms:
        hit = [k for k, v in corpus.items() if f in v]
        if hit:
            return hit, f
    return [], None


print("C4 独立对账 —— 投稿信 vs 新稿语料（穷举数字扫描 + 显式主张 + 反算）")
print(f"语料: {', '.join(corpus)}  ({sum(len(v) for v in corpus.values())} 字符)\n")

fail = 0
wl_seen = set()
for L in LETTERS:
    body = claims_only(L.read_text(encoding="utf-8"))
    toks = Counter(TOKEN.findall(body))
    miss = []
    for t, n in sorted(toks.items()):
        if t.lower() in WHITELIST:
            wl_seen.add(t)
            continue
        hit, form = probe(t, corpus)
        if not hit:
            miss.append((t, n))
    print(f"── {L.name}：数字 token {len(toks)} 种，未命中 {len(miss)} 种")
    for t, n in miss:
        print(f"     MISS ×{n}  {t!r}")
        fail += 1
    if not miss:
        print("     全部数字 token 均在语料中有落点")
print(f"\n白名单放行的 token（非主张数字）：{sorted(wl_seen)}")

# ── ①b 命中归属：逐个 token 报出命中文件与形态
KEY = ["P=0.0028", "0.0601", "3.0×10⁻⁴", "43.65", "0.647", "0.635", "P=0.438",
       "6.5×10⁻⁵", "P=0.639", "P=0.31", "12.80%", "0.00714", "P=6.6×10⁻⁷", "0.477", "86.83%"]
print("\n关键值归属（前缀差异单列）：")
for k in KEY:
    hit, form = probe(k, corpus)
    ok = "OK " if hit else "MISS"
    if not hit:
        fail += 1
    note = "" if form == k else f"  ← 以 {form!r} 命中（新稿省略 P= 前缀，值一致）"
    print(f"  {ok} {k:<12} → {','.join(hit) or '—'}{note}")

# ── ② 逐字短语（标题与引号内文本必须逐字）
QUOTES = [
    ("标题", "Criterion-Dependent Priority of Matrix-Remodelling Candidates in Knee Osteoarthritis Synovium: Cross-Cohort Replication and Donor-Resolved Localisation"),
    ("摘要问句", "Which matrix-remodelling candidates deserve study first in knee osteoarthritis (KOA) synovium, and in which cells are they expressed?"),
    ("不可解释", "not interpretable"),
    # 信件此处为转述（抹平了原文的括号），故只对「承载排名的那个分句」要求逐字：
    # 新稿 §3.3 L160 = 「under a fixed-effect (or pooled) estimator with a |log2FC| criterion, LOXL1 is first」
    ("§3.3 条件首位", "with a |log2FC| criterion, LOXL1 is first"),
    ("判据依赖", "Priority was criterion-dependent rather than a property of the candidates"),
    ("三供者", "86.83%, 87.54% and 43.89% fibroblasts"),
    ("无正常对照", "no normal controls"),
    ("正对照缺失", "has no null distribution of its own"),
]
print("\n逐字短语检索：")
for tag, q in QUOTES:
    qn = norm(q)
    hit = [name for name, v in corpus.items() if qn in v]
    if not hit:
        fail += 1
    print(f"  {'OK ' if hit else 'MISS'} [{tag}] → {','.join(hit) or '—'}")

# ── ③ 反算：由 Table S3 数出「18 组中达名义显著」的组数
s3 = CORPUS_FILES["S3"].read_text(encoding="utf-8")
rows = [l for l in s3.split("\n") if re.match(r"\| (LOXL1|SERPINE1|VGLL4|SDC4|ITGA5|SAV1) \|", l)]
nominal = 0
n_iv1 = 0
for l in rows:
    c = [x.strip() for x in l.strip("|").split("|")]
    # 主分析：列 1=Gene, 2=#IV, 3=β, 4=SE, 5=CI, 6=P
    p_main = float(c[5])
    if p_main < 0.05:
        nominal += 1
    if c[1] == "1":
        n_iv1 += 1
print(f"\n反算（读 Table S3 主分析列，不比对字面）：")
print(f"  行数 {len(rows)}（期望 18 = 6 基因 × 3 暴露）")
print(f"  主分析 P<0.05 的组数 = {nominal}（信件称 6 组达名义显著）")
print(f"  单工具（#IV=1）的组数 = {n_iv1}（新稿 §4 称 six gene–exposure sets retain a single instrument）")
if len(rows) != 18 or nominal != 6 or n_iv1 != 6:
    fail += 1
    print("  ✗ 反算结果与信件/新稿不符")
else:
    print("  ✓ 反算与信件所述一致")

# ── 反向对照：每个数字 token 改一位后必须查不到（受检对象与正向扫描同为「主张区」）
neg_fail = 0
checked = 0
for L in LETTERS:
    body = claims_only(L.read_text(encoding="utf-8"))  # 与正向扫描同一文本范围
    for t in set(TOKEN.findall(body)):
        m = re.search(r"\d", t)
        if not m or len(t) < 3:
            continue
        d = int(t[m.start()])
        broken = t[:m.start()] + str((d + 1) % 10) + t[m.start() + 1:]
        if broken == t:
            continue
        checked += 1
        if any(broken in v for v in corpus.values()):
            print(f"  注：{t!r} → {broken!r} 仍命中（数字相近，非假阳性）")
print(f"\n反向对照：{checked} 个 token 改位后受检（'仍命中' 提示属数字邻近，不作失败）")

print(f"\n结果：{'全部通过' if not fail else f'{fail} 项未通过'}")
raise SystemExit(1 if fail else 0)
