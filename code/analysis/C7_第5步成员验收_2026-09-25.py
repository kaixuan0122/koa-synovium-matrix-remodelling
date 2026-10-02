#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
C7 — 第 5 步「五处同步」的成员验收审计

底稿 `第5步_五处同步_清单与增量_2026-09-23.md` §三·7 的方法要求：
「凡成员文件里的 § 引用须逐条回新稿核对后重写，**不得批量替换**——旧 §3.x 有 7 节、
新 §3.x 只有 4 节，编号不能按偏移推。」

因此本审计不做「编号是否 ≤ 最大节号」这类偏移判据，而是：
  ① 把新稿的节标题集合当作**白名单**，逐成员列出每一处 § 引用及其所在行，
     逐个判存在性；成员文件里出现白名单外的节号即失败。
  ② 陈旧字串扫描：硬失败集（无条件错误）+ 复核集（合法与否要看上下文，逐行打印供判读）。
  ③ 产物 md5 台账落盘，供回执引用。

两处刻意的排除（避免假阳性）：
  · `08_投稿材料/00_投稿清单.md` 的 § 引用是**自引用**（它自己的 §2.2、§1.2、§二·八），
    不是指新稿，故逐条判存在性对它无意义——单列并说明。
  · 投稿信经实测**不含任何 §**；若有命中即为失败（此前一轮 BSD grep 曾误报，见回执）。
"""
import hashlib
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[3]
REV = ROOT / "08_修订_2026-09-17"
MS = REV / "05_稿件" / "English_manuscript_MDPI_重定位2026-09-23.md"

# 成员清单（本轮在册者；JIR 稿与 JIR 投稿清单按底稿 §五·二 不同步，不在册）
MEMBERS = [
    ("1  MDPI 英文稿", REV / "05_稿件" / "English_manuscript_MDPI_重定位2026-09-23.md"),
    ("2  中文稿", REV / "05_稿件" / "双语初稿_候选优选与组织来源_重定位2026-09-23.md"),
    ("6  S1 英文", REV / "06_补充材料" / "Supplementary_Table_S1.md"),
    ("6  S1 中文", REV / "06_补充材料" / "Supplementary_Table_S1_中文.md"),
    ("7  S2 英文", REV / "06_补充材料" / "Supplementary_Table_S2.md"),
    ("7  S2 中文", REV / "06_补充材料" / "Supplementary_Table_S2_中文.md"),
    ("8  S3", REV / "06_补充材料" / "Supplementary_Table_S3.md"),
    ("9  STROBE-MR", REV / "06_补充材料" / "STROBE-MR_checklist.md"),
    # JIR 母版信 2026-09-27 归档处置：成员身份不变（REGISTER 的 tag 照旧），只在册路径移到
    # _archive_ 子目录。审计照样覆盖它——归档件同样要过 § 引用与陈旧字串扫描。
    ("4  投稿信 JIR 母版", REV / "08_投稿材料" / "_archive_JIR作废_2026-09-27" / "Cover_Letter.md"),
    ("5  投稿信 MDPI 母版", REV / "08_投稿材料" / "Cover_Letter_MDPI.md"),
]
SELF_REF_ONLY = [("10 投稿清单", REV / "08_投稿材料" / "00_投稿清单.md")]

ms = MS.read_text(encoding="utf-8")
# 注意：顶层写 `## 1. Introduction`（编号后有句点），子节写 `### 2.1 Data sources`（编号后无句点），
# 故句点必须可选，否则子节全部漏进白名单、把每一处子节引用都误判为不存在。
SECS = set(re.findall(r"^#{2,3} (\d+(?:\.\d+)?)\.?\s", ms, re.M))
assert len(SECS) >= 17, f"节标题抓取只得到 {sorted(SECS)}，正则或稿件结构有变，勿据此判成员"
# 裸 §2（指 Methods 整层）在 STROBE-MR item 9b 中是合法用法，已逐条核实（C1 记录）
SECS_OK = SECS | {"2"}

fail = 0
print("C7 第 5 步成员验收审计")
print(f"新稿节标题白名单：{sorted(SECS, key=lambda s: [int(x) for x in s.split('.')])}")
print(f"（另加裸 §2 = Methods 整层，STROBE-MR item 9b 的合法用法，逐条核实见 C1 日志）\n")

# ── ① § 锚点逐条判存在性
print("① § 锚点存在性（每一处引用逐条判，按行归属）")
for tag, p in MEMBERS:
    if not p.exists():
        print(f"  FAIL [{tag}] 文件不存在：{p.name}")
        fail += 1
        continue
    lines = p.read_text(encoding="utf-8").split("\n")
    hits = []
    for i, l in enumerate(lines):
        for tok in re.findall(r"§(\d+(?:\.\d+)?)", l):
            hits.append((i + 1, tok, l.strip()))
    if not hits:
        print(f"  OK   [{tag}] 无 § 引用（0 处）")
        continue
    bad = [(n, t, l) for n, t, l in hits if t not in SECS_OK]
    uniq = sorted({t for _, t, _ in hits}, key=lambda s: [int(x) for x in s.split(".")])
    if bad:
        fail += len(bad)
        for n, t, l in bad:
            print(f"  FAIL [{tag}] L{n} §{t} 不在新稿节标题中：{l[:70]}")
    else:
        print(f"  OK   [{tag}] § 引用 {len(hits)} 处 / {len(uniq)} 个编号 "
              f"{uniq} 全部存在于新稿")

# 投稿信的强判据：不含任何 §（此前 BSD grep 误报过，故用 python 复核）
for tag, p in MEMBERS:
    if "投稿信" in tag:
        n = p.read_text(encoding="utf-8").count("§")
        if n:
            fail += 1
        print(f"  {'OK  ' if n == 0 else 'FAIL'} [{tag}] `§` 命中 {n} 处（预期 0）")

for tag, p in SELF_REF_ONLY:
    lines = p.read_text(encoding="utf-8").split("\n")
    n = sum(len(re.findall(r"§\d", l)) for l in lines)
    refs = sorted({t for l in lines for t in re.findall(r"§(\d+(?:\.\d+)?)", l)})
    print(f"  ——   [{tag}] § 引用 {n} 处为**自引用**（本文件自己的 §2.2 / §1.2 / §二·八），"
          f"非指新稿，故不判存在性；涉及编号 {len(refs)} 个")

# ── ② 陈旧字串扫描
HARD = [
    ("Cross-Cohort Expression Profile and Cellular Localization", "2026-09-17 版英文标题"),
    ("Genetically Causal", "已撤下的因果主张（§2.2 返工）"),
    ("connective-tissue cis-eQTL", "已更正的暴露标识（§2.2 返工）"),
    ("PP.H4=0.783", "旧共定位值（§2.3 已改）"),
]
REVIEW = [
    ("PP.H4=0.742", "旧共定位值——新稿未报告共定位结果，信件中应已删"),
    ("a priori", "面板改称 knowledge-guided exploratory 后不应再出现"),
    ("permutation calibration", "已改称 gene-set calibration；但『不是标签置换检验』一句合法"),
    ("cell-type-restricted", "已被推翻的措辞（2026-09-17 第五轮）"),
    ("positive control", "新稿命中 0；信件中应已随整段替换删除"),
    ("315th", "旧正对照名次，新稿未报告"),
]
print("\n② 陈旧字串扫描（硬失败集 / 复核集）")
for s, why in HARD:
    hits = [(tag, i + 1) for tag, p in MEMBERS
            for i, l in enumerate(p.read_text(encoding="utf-8").split("\n")) if s in l]
    ok = not hits
    if not ok:
        fail += len(hits)
    print(f"  {'OK  ' if ok else 'FAIL'} {s[:44]:<46} {len(hits)} 处  ({why})")
print("  —— 复核集（登记制）")
# 复核集的 14 处命中已逐条人工判读，结论固化在下面的登记表里。
# 判读一次就打印一次是没有守卫的：文件日后被改动、命中数变了，没人会重看。
# 故改为登记制——未登记的命中、登记锚不符的命中、登记了却消失的命中，三者都失败。
# 列义：(成员, 陈旧串, 行号, 行内锚, 判读结论)
REGISTER = [
    ("4  投稿信 JIR 母版", "PP.H4=0.742", 21,
     "the revised manuscript reports no colocalization result",
     "合法：回应审稿人时援引旧值以说明「新稿不再报共定位」，非把旧值当现值"),
    ("5  投稿信 MDPI 母版", "PP.H4=0.742", 25,
     "the revised manuscript reports no colocalization result",
     "合法：同上"),
    # 行号随文件改动漂移，**不参与配对**（配对看锚）。下列行号仅供人读：
    # 2026-09-27 第八轮两次改写两稿文件头注释，三处由 86/158/69 → 89/161/74 → 94/166/79。
    # 锚点自始至终未变，判读结论不变。
    ("1  MDPI 英文稿", "a priori", 94, "A curated",
     "合法：显式否定用法（\"a priori\" is not used）"),
    ("1  MDPI 英文稿", "a priori", 166, "Third, the candidate with the smallest discovery P",
     "合法：**子串假阳性**——该行是 a priority 一词，恰好含 a priori 字面"),
    ("2  中文稿", "a priori", 79, "机械/刚度基因集",
     "合法：同英文稿那处的否定句，中文「不使用『a priori』一词」"),
    ("7  S2 英文", "positive control", 1, "Supplementary Table S2",
     "合法：S2 题名自带该小节名，非陈旧引用"),
    ("7  S2 英文", "positive control", 12, "S2a. Per-cohort unbiased feature selection",
     "合法：S2 自有小节标题"),
    ("7  S2 英文", "positive control", 34, "GSE55235 is the only cohort that can be interpreted",
     "合法：S2 内部对阳性对照队列的说明"),
    ("7  S2 中文", "positive control", 1, "Supplementary Table S2",
     "合法：同 S2 英文 L1（中文版题名沿用英文小节名）"),
    ("4  投稿信 JIR 母版", "positive control", 22, "machine-learning positive control in GSE55235",
     "合法：陈述旧稿把 GSE55235 当阳性对照这一缺陷，非主张现状"),
    ("5  投稿信 MDPI 母版", "positive control", 26, "machine-learning positive control in GSE55235",
     "合法：同上"),
    ("7  S2 英文", "315th", 35, "The most likely reason is",
     "合法：S2 自己报告的旧判据名次"),
    ("4  投稿信 JIR 母版", "315th", 22, "machine-learning positive control in GSE55235",
     "合法：同上（引旧名次说明该判据的缺陷）"),
    ("5  投稿信 MDPI 母版", "315th", 26, "machine-learning positive control in GSE55235",
     "合法：同上"),
]
seen = set()
for s, why in REVIEW:
    for tag, p in MEMBERS:
        for i, l in enumerate(p.read_text(encoding="utf-8").split("\n")):
            if s in l:
                # 配对键是**内容锚**，不是行号：判定是否过期取决于该处内容是否变了，行号只是位置。
                # 2026-09-27 前本键含行号，一轮之内被文件头注释的增删触发了两次假失败（同一处
                # 内容、结论完全没变），故去掉行号；行号仍打印，漂移时显式提示，便于人判读。
                reg = next((r for r in REGISTER
                            if r[0] == tag and r[1] == s and r[3] in l), None)
                if reg is None:
                    near = [r for r in REGISTER if r[0] == tag and r[1] == s]
                    fail += 1
                    hint = f"（该文件另有已登记 {s!r} {len(near)} 处，但行内锚均不符）" if near else ""
                    print(f"     ✗ 未登记 [{tag}] L{i + 1} {s!r}：{l.strip()[:90]}{hint}")
                    continue
                seen.add((tag, s, reg[3]))
                drift = "" if reg[2] == i + 1 else f"  ← 行号由 L{reg[2]} 漂移到 L{i + 1}（锚相符，判定仍成立）"
                print(f"     ✓ [{tag}] L{i + 1} {s!r} — {reg[4]}{drift}")
for tag, s, n, anchor, verdict in REGISTER:
    if (tag, s, anchor) not in seen:
        fail += 1
        print(f"     ✗ 登记项已消失 [{tag}] {s!r}（锚 {anchor[:40]!r}）——该处内容已变或被删，须重新判读")

# ── ③ 产物 md5 台账
LEDGER = REV / "01_审计" / "第5步产物md5.tsv"
rows = []
for tag, p in MEMBERS + SELF_REF_ONLY:
    if p.exists():
        rows.append((tag, p.relative_to(REV).as_posix(), hashlib.md5(p.read_bytes()).hexdigest(),
                     p.stat().st_size))
LEDGER.write_text("成员\t路径\tmd5\t字节\n" + "\n".join("\t".join(map(str, r)) for r in rows) + "\n",
                  encoding="utf-8")
print(f"\n③ 产物 md5 台账：{LEDGER.relative_to(REV)}（{len(rows)} 行）")
for tag, path, h, sz in rows:
    print(f"  {h}  {sz:>7} B  {path}")

print(f"\n结果：{'全部通过' if not fail else f'{fail} 项未通过'}")
raise SystemExit(1 if fail else 0)
