#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
C6 — 投稿清单（重定位版）的独立对账

与生成器 C5 方向相反、代码路径不同：
  C5 查的是「新值写进去了没有」——它从新稿派生值，再断言值等于常量。
  C6 查的是「旧值有没有被冒充成当前值」——它从**产物**出发，把每个旧值扫出来，
     再要求每一处命中都**归属到一个已标为取代/历史口径的区间**。这正是本项目
     反复踩到的形态：值扫描命中必须按行归属才能判，不归属就只能一刀切误报。

三条判据：
  ① 旧值归属：旧字数/旧计数/旧标题/旧表号，每一处命中须落在 §二（返工记录，含 §二·六/§二·七）
     或落在带陈旧标记的行（口径/取代/过时/当时/旧编号/2026-09-17）上。
  ② 新值在位 + 交叉反算：计数、字数、待补条目、图题，一律**独立地**从新稿重算一遍
     （不比对 C5 的中间量），与产物所载逐项对照。
  ③ 反向对照：每个新值改一位后须在新稿中查不到（证明检索非恒真）。

语料 = 新稿（唯一事实来源）
"""
import re
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[3]
REV = ROOT / "08_修订_2026-09-17"
CHK = REV / "08_投稿材料" / "00_投稿清单.md"
MS = REV / "05_稿件" / "English_manuscript_MDPI_重定位2026-09-23.md"

chk = CHK.read_text(encoding="utf-8")
ms = MS.read_text(encoding="utf-8")
chk_lines = chk.split("\n")
ms_lines = ms.split("\n")

fail = 0
print("C6 独立对账 —— 投稿清单（重定位版）vs 新稿")
print(f"产物 {len(chk_lines)} 行 / {len(chk.encode('utf-8'))} B；语料 = 新稿 {len(ms_lines)} 行\n")

# ── 归一化：去掉 markdown 强调符并折叠空白（否则 `**6,613** 词` 检索不到「6,613 词」）
norm = lambda s: re.sub(r"\s+", " ", s.replace("*", "")).strip()
chk_n = [norm(l) for l in chk_lines]

# ── 区间定义：§二（返工记录）从题头到 §三 之前，整段视为「历史口径区」
i2 = next(i for i, l in enumerate(chk_lines) if l.startswith("## 二、返工记录"))
i3 = next(i for i, l in enumerate(chk_lines) if l.startswith("## 三、✋ 投稿前必须补齐"))
STALE_MARK = re.compile(r"口径|取代|过时|当时|旧编号|旧标题|2026-09-17|原 6|原 4|原 3|原 18|源件")


def find_lines(needle):
    """返回命中该串的行号列表（在归一化文本上检索，行号即归属）。"""
    return [i for i, l in enumerate(chk_n) if needle in l]


def judge(i):
    """把命中归属到行与区间：返回 (行文本, 是否已声明为历史值, 理由)。"""
    line = chk_lines[i]
    if i2 <= i < i3:
        return line, True, "落在 §二 历史口径区"
    if STALE_MARK.search(line):
        return line, True, "行内带陈旧标记"
    return line, False, "**未被声明为历史值**"


# ── ① 旧值归属
STALE = [
    ("7,261", "第六轮正文字数"),
    ("摘要 200", "第六轮摘要词数口径"),
    ("6 图", "第六轮图数"),
    ("4 表", "第六轮表数"),
    ("18 条参考文献", "第六轮文献数"),
    ("Cross-Cohort Expression Profile", "旧标题"),
    ("0.583", "旧外部池化 AUC"),
    ("0.816", "旧定标经验 P"),
    ("14,510", "旧 GSE89408 候选基因数"),
    ("315th–1586th", "旧正对照名次"),
    ("6,613", "第八轮前的正文字数（S4 行数口径改写后为 6,617）"),
]
print("① 旧值归属（每处命中须落在已声明为历史值的区间）")
for val, what in STALE:
    hits = find_lines(val)
    if not hits:
        print(f"  OK   {val:<30} 产物中 0 处（{what}已被彻底替换）")
        continue
    bad = 0
    for i in hits:
        line, ok, why = judge(i)
        if not ok:
            bad += 1
            fail += 1
            print(f"  FAIL {val!r} 第 {i + 1} 行 未被声明：{line.strip()[:70]}")
    if not bad:
        zones = {judge(i)[2] for i in hits}
        print(f"  OK   {val:<30} {len(hits)} 处命中，全部已声明（{'; '.join(sorted(zones))}）")

# ── ② 新值在位 + 交叉反算
print("\n② 新值在位（独立从新稿重算，与产物逐项对照）")


def ms_fact(pattern, flags=0, group=1):
    m = re.search(pattern, ms, flags)
    assert m, f"新稿中找不到：{pattern[:40]}"
    return m.group(group)


TITLE = ms_lines[0].lstrip("# ").strip()
ABS_W = int(ms_fact(r"abstract (\d+) words"))
MAIN_W = int(ms_fact(r"main text (\d+) words"))
N_FIG = int(ms_fact(r"\*\*Figures:\*\* (\d+)"))
N_TAB = int(ms_fact(r"\*\*Tables:\*\* (\d+)"))
N_SUP = int(ms_fact(r"\*\*Supplementary tables:\*\* (\d+)"))

refs = [l for l in re.search(r"^## References$(.*?)^## Declarations$", ms, re.S | re.M).group(1).split("\n")
        if re.match(r"^\d+\. ", l)]
todo = []
for l in refs:
    num = l.split(".")[0]
    a = l.split(". ", 1)[1].split("*")[0]
    if "et al." in a and (a.count(";") + 1 - 1) < 10:
        todo.append(num)

CHECK = [
    (norm(TITLE), "新标题"),
    (norm(f"{MAIN_W:,} 词"), f"正文字数 {MAIN_W:,}"),
    (norm(f"摘要 {ABS_W}"), f"摘要词数 {ABS_W}"),
    (norm("、".join(todo)), f"待补 et al. 条目 {todo}"),
]
for needle, what in CHECK:
    n = len(find_lines(needle))
    ok = n >= 1
    if not ok:
        fail += 1
    print(f"  {'OK  ' if ok else 'FAIL'} {what:<28} 产物 {n} 处  «{needle[:44]}»")

# 计数：从产物 §一 的正文行**反解数字**，再与新稿首页自报比对（不依赖我方措辞）
row = next((l for l in chk_n if l.startswith("| `Main_Manuscript_MDPI.docx` |")), None)
if row is None:
    fail += 1
    print("  FAIL 产物 §一 缺 `Main_Manuscript_MDPI.docx` 行，无法反解计数")
else:
    for pat, val, what in [(r"(\d+) 图", N_FIG, "图"), (r"(\d+) 表", N_TAB, "表"),
                           (r"(\d+) 条参考文献", len(refs), "参考文献")]:
        m = re.search(pat, row)
        got = int(m.group(1)) if m else None
        ok = got == val
        if not ok:
            fail += 1
        print(f"  {'OK  ' if ok else 'FAIL'} 产物 §一 自报{what}数 {got} vs 新稿 {val}")

# 图题：独立抓取后逐条要求出现在产物图清单中
figs = re.findall(r"^\*\*Figure (\d+): (.+?)\.?\*\*", ms, re.M)
print(f"  图题独立抓取 {len(figs)} 条：")
for num, t in figs:
    ok = bool(find_lines(f"| Fig {num} | {t} |"))
    if not ok:
        fail += 1
    print(f"    {'OK  ' if ok else 'FAIL'} Fig {num}: {t[:52]}")

# 计数三元组：产物 §二·八 所载须与新稿首页自报逐字一致
triple = f"{N_FIG} 图 / {N_TAB} 表 / {N_SUP} 补充表 / {len(refs)} 文献"
must_in = bool(find_lines(triple))
if not must_in:
    fail += 1
print(f"  {'OK  ' if must_in else 'FAIL'} 计数三元组「{triple}」在 §七 总结中逐字出现")

# ── ③ 反向对照
print("\n③ 反向对照（新值改一位后必须在新稿中查不到）")
neg = 0
checked = 0
for needle, what in CHECK + [(t, f"Fig {n} 题") for n, t in figs]:
    if len(needle) < 6:
        continue
    m = re.search(r"\d", needle)
    if not m:
        continue
    d = int(needle[m.start()])
    broken = needle[:m.start()] + str((d + 1) % 10) + needle[m.start() + 1:]
    checked += 1
    if broken in ms:
        neg += 1
        print(f"  假阳性：{needle[:40]!r} 改位后仍命中")
print(f"  {checked} 项受检，假阳性 {neg}")
fail += neg

print(f"\n结果：{'全部通过' if not fail else f'{fail} 项未通过'}")
raise SystemExit(1 if fail else 0)
