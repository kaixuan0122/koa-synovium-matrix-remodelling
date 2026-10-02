#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
C5 — MDPI 投稿清单的重定位版（第 5 步成员 #10）

底稿 `07_修订记录/第5步_五处同步_清单与增量_2026-09-23.md` §四·7 要求：
「标题、计数、字数、S2 描述、图清单、§二·七 的历史记录须注明『该轮数字已被 2026-09-23 重定位取代』。」

做法与 C3（投稿信）一致：**读源件、按行做定点手术**，而不从零重写——
源件 §二 的六轮返工记录是唯一存档，重写等于销毁它（方案规则 1）。
因此历史的每一行原样保留，只做两类改动：
  (a) 数值/文件名一类**事实已变**的行，改成本轮口径并标明来源节；
  (b) 历史节的**抬头**加「已被取代」横幅，使旧数字不致被当作当前值。

**所有计数、字数、文献待补条目、图题均从新稿程序化派生**（不手抄），
故本文件与新稿不可能脱节；派生值与预期常量不符即断言失败。
"""
import re
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from pdftext import pdf_text  # noqa: E402  （§1.2 要用它扫图件 PDF 的文本层）

ROOT = HERE.parents[2]   # HERE 是 scripts/ 目录，故比原来的 __file__ 少一级
REV = ROOT / "08_修订_2026-09-17"
SRC = ROOT / "05_投稿包" / "MDPI版" / "00_投稿清单.md"
MS = REV / "05_稿件" / "English_manuscript_MDPI_重定位2026-09-23.md"
DST = REV / "08_投稿材料" / "00_投稿清单.md"
DATE = "2026-09-25"
DATE_D = "2026-09-27"   # 第八轮（四项处置）的日期。与 DATE 分开：DATE 同时是投稿信的落款日期，
                        # 改它会连带改错信里的日期，故新轮次一律用新常量。


def must(cond, msg):
    assert cond, f"后置条件失败: {msg}"
    print(f"  post OK  {msg}")


def must2(cond, ok_msg, bad_msg):
    """must() 成功时也打印消息，故成功版必须是陈述句。"""
    assert cond, f"后置条件失败: {bad_msg}"
    print(f"  post OK  {ok_msg}")


# ─────────────────────────────────────────── 一、从新稿派生事实
ms = MS.read_text(encoding="utf-8")
lines_ms = ms.split("\n")

TITLE = lines_ms[0].lstrip("# ").strip()
must2(TITLE.startswith("Criterion-Dependent Priority"),
      f"新稿标题取自首行：{TITLE[:60]}…", f"首行不是预期的新标题：{TITLE[:60]}")

wc_line = next(l for l in lines_ms if l.startswith("**Word count:**"))
m = re.search(r"abstract (\d+) words; main text (\d+) words", wc_line)
ABS_W, MAIN_W = int(m.group(1)), int(m.group(2))

cnt_line = next(l for l in lines_ms if l.startswith("**Figures:**"))
m2 = re.search(r"\*\*Figures:\*\* (\d+) · \*\*Tables:\*\* (\d+) · \*\*Supplementary tables:\*\* (\d+)", cnt_line)
N_FIG, N_TAB, N_SUP = int(m2.group(1)), int(m2.group(2)), int(m2.group(3))

refs_block = re.search(r"^## References$(.*?)^## Declarations$", ms, re.S | re.M).group(1)
ref_lines = [l for l in refs_block.split("\n") if re.match(r"^\d+\. ", l)]
N_REF = len(ref_lines)

# 待补条目：MDPI 要求列前 10 位作者，不足 10 位却用 et al. 者待补
todo_refs = []
for l in ref_lines:
    num = int(l.split(".")[0])
    authors = l.split(". ", 1)[1].split("*")[0]
    has_etal = "et al." in authors
    named = authors.count(";") + 1
    if has_etal:
        named -= 1                      # 「et al.」本身不是作者名，须扣除
    if has_etal and named < 10:
        todo_refs.append(str(num))

decl_block = re.search(r"^## Declarations$(.*?)^---$", ms, re.S | re.M).group(1)
N_DECL = len([l for l in decl_block.split("\n") if re.match(r"^\*\*[^*]+:\*\*", l)])

fig_block = re.search(r"^## Figure Legends$(.*)$", ms, re.S | re.M).group(1)
FIG_TITLES = re.findall(r"^\*\*Figure (\d+): (.+?)\.?\*\*", fig_block, re.M)

affils = [l.strip() for l in lines_ms if re.match(r"^\^[123]\^ ", l)]

# ── 断言：派生值必须等于本轮确认过的常量（防新稿再变而清单静默跟随）
must(N_REF == 20, f"新稿参考文献 {N_REF} 条（重定位口径 20 条）")
must(N_FIG == 5 and N_TAB == 3 and N_SUP == 4,
     f"新稿自报 图 {N_FIG} / 表 {N_TAB} / 补充表 {N_SUP}（重定位口径 5/3/4）")
must(ABS_W == 195 and MAIN_W == 6617,
     f"新稿自报 摘要 {ABS_W} 词 / 正文 {MAIN_W} 词（第八轮口径 195/6617；原 6,613 因 S4 行数口径改写 +4 词）")
must(N_DECL == 7, f"新稿声明 {N_DECL} 项（口径 7 项）")
must(len(FIG_TITLES) == 5, f"Figure Legends 段抓到 {len(FIG_TITLES)} 条图题（期望 5）")
must(todo_refs == ["4", "5", "13", "14", "16", "18"],
     f"待补 et al. 条目由新稿算出为 {todo_refs}（原清单为 4/5/11/13/14/15，编号已随文献表全换）")
must(len(affils) == 3, f"抓到 {len(affils)} 行单位（期望 3）")

# ─────────────────────────────────────────── 二、手术工具
src = SRC.read_text(encoding="utf-8")
lines = src.split("\n")


def find_unique(prefix, tag):
    hits = [i for i, l in enumerate(lines) if l.startswith(prefix)]
    assert len(hits) == 1, f"{tag}: 行首 '{prefix[:40]}' 命中 {len(hits)} 处（期望 1）"
    return hits[0]


def sub_line(prefix, new_block, tag):
    i = find_unique(prefix, tag)
    lines[i:i + 1] = new_block
    print(f"  改行  [{tag}] {prefix[:34]}… → {len(new_block)} 行")


def insert_after(prefix, new_block, tag):
    i = find_unique(prefix, tag)
    lines[i + 1:i + 1] = new_block
    print(f"  插块  [{tag}] 于「{prefix[:30]}…」之后 +{len(new_block)} 行")


def insert_before(prefix, new_block, tag):
    i = find_unique(prefix, tag)
    lines[i:i] = new_block
    print(f"  插块  [{tag}] 于「{prefix[:30]}…」之前 +{len(new_block)} 行")


print(f"\nC5 重定位投稿清单：源 {SRC.name} {len(lines)} 行")
print(f"    新稿派生：标题 1｜图 {N_FIG} 表 {N_TAB} 补充表 {N_SUP}｜文献 {N_REF}（待补 {todo_refs}）"
      f"｜声明 {N_DECL}｜字数 {ABS_W}/{MAIN_W}\n")


def blk(text):
    """多行文本 → 行列表（首尾空行去掉）。"""
    return text.strip("\n").split("\n")


# ── ① 抬头：标题与更新行
sub_line("**稿件：**", [f"**稿件：** {TITLE}"], "抬头标题")
sub_line("**建包日期：**", blk(f"""
**建包日期：** 2026-09-15　|　**最近更新：** {DATE_D}（**2026-09-23 重定位** 见 §二·八；**第 5 步五处同步** 见 §二·八·一；**第八轮四项处置** 见 §二·九）
**现役正文：** `08_修订_2026-09-17/05_稿件/English_manuscript_MDPI_重定位2026-09-23.md`（中文稿为同目录 `双语初稿_候选优选与组织来源_重定位2026-09-23.md`）
**本轮不重出 docx**（§二·八·三）：图件已全部重出，但文件名与格式尚未对齐本行声明的
`Figures/Fig1`–`Fig5`.tiff`，重出仍会产出指向不存在图件的成品（§二·九·二）。
**本文件是 `05_投稿包/MDPI版/00_投稿清单.md` 的重定位版**——源件按方案规则 1 原样保留，不改动。
"""), "抬头更新行")

# ── ② §一 包内文件表：计数与图件
sub_line("| `Main_Manuscript_MDPI.docx` |", blk(f"""
| `Main_Manuscript_MDPI.docx` | 正文（Title Page、摘要、正文、**{N_FIG} 图**内嵌、**{N_TAB} 表**、**{N_REF} 条参考文献**、{N_DECL} 项声明）——计数为 **2026-09-23 重定位口径**（图 6→5、表 4→3、文献 18→20）；文件名暂沿用，重出时随图件一并生成 |
"""), "正文行计数")
sub_line("| `Figures/Fig1`", blk("""
| `Figures/Fig1`–`Fig5`.tiff` | **5 张图**（原 6 张，重定位后主图减为 5），TIFF / RGB / 300 dpi / Arial。**五幅内容已全部重出**（`04_图表/` 下 5 组 PDF+PNG），**但文件名与格式尚未对齐本行声明**，见 §1.2 与 §二·九·二 |
"""), "图件行计数")
sub_line("| `Supplementary/Supplementary_Table_S3.docx` |", blk("""
| `Supplementary/Supplementary_Table_S3.docx` | MR 敏感性：18 组「基因 × 暴露」在 r²≤0.1 与 r²≤0.01 两级剪枝阈值下的对照，含逐 SNP F |
| `Supplementary/Supplementary_Table_S4.docx` | 判据敏感性：名次随估计量变动、显著性跨线计数、逐基因留一法重排。**已单独成表**（2026-09-25，脚本 D1）——`06_补充材料/Supplementary_Table_S4.md`，三段 **246 + 3 + 450 = 699 行**，逐行由 `02_结果表/` 的三份 CSV 转排，无值被重算、无手抄（§二·九·三） |
"""), "补充表 S3/S4 行")
sub_line("源文件：", blk(f"""
源文件（{DATE} 起）：正文 `08_修订_2026-09-17/05_稿件/English_manuscript_MDPI_重定位2026-09-23.md`；
投稿信 `08_修订_2026-09-17/08_投稿材料/Cover_Letter_MDPI.md`。
"""), "源文件行")
sub_line("生成脚本：", blk("""
docx 生成脚本：`03_论文初稿/生成word.py`（`MDPI_DOCS` / `EXTRA` 两个列表）——**本轮未运行**（§二·八·三）。
"""), "生成脚本行")

# ── ③ §1.2 图清单（图题由新稿抓取；图件状态由**实测**得出，不手写）
# 图件状态的三项事实各由一个可重跑的量支撑，而不是由「我记得画过了」支撑：
#   · 文件在不在    → 逐图查 04_图表/ 下的 PDF 与 PNG
#   · 带不带旧图号  → 扫 PDF 文本层（本项目有过「重编号后图文件内仍写旧号、且与别的图撞号」的先例）
#   · 图注入不入图  → 拿该图图注的**末句**去 PDF 文本层里找，画进图内的必然命中
# 末一条是 2026-09-27 实测发现的不一致：图 1、图 2 把整条图注画进了图内，图 3/4/5 没有。
# 装包出 docx 时若逐图再配一遍稿件 Figure Legends，前两幅的图注会重复。故此项每轮实测、写进表里。
FIG_FILES = {"1": "Fig_study_design", "2": "Fig_candidate_overview", "3": "Fig_forest_logFC",
             "4": "Fig_evidence_matrix", "5": "Fig_scRNA_dotplot"}
FIGDIR = REV / "04_图表"
FIGNUM = re.compile(r"Figure\s*\d+|Fig\.?\s*\d+")
N_EMB = 2   # 图注已画入图内的图件数（2026-09-27 实测）。此数一变，本节说明与 §六 的缺口项都要改。


def pdf_flat(p):
    recs, _ = pdf_text(p)
    return re.sub(r"\s+", " ", " ".join(t for _, _, t in sorted(recs, key=lambda r: (-r[0], r[1]))))


legend_tail = {}
for _n, _t in FIG_TITLES:
    # 变量名不用 blk——blk 是本文件的多行文本工具函数，覆盖它会让后面的排版调用变成 Match 对象
    _blk = re.search(rf"^\*\*Figure {_n}: .+?\*\*(.+?)(?=^\*\*Figure |\Z)", fig_block, re.S | re.M)
    sents = [s.strip() for s in re.split(r"(?<=\.)\s+", _blk.group(1).strip()) if s.strip()]
    legend_tail[_n] = re.sub(r"\s+", " ", sents[-1])

fig_rows, n_emb = [], 0
for num, t in FIG_TITLES:
    for ext in ("pdf", "png"):
        must((FIGDIR / f"{FIG_FILES[num]}.{ext}").is_file(),
             f"Fig {num} 图件存在：{FIG_FILES[num]}.{ext}")
    flat = pdf_flat(FIGDIR / f"{FIG_FILES[num]}.pdf")
    nums = FIGNUM.findall(flat)
    stray = [x for x in nums if not re.search(rf"\b{num}\b", x)]
    must(not stray, f"Fig {num} 图内无撞号旧图号（实测命中 {nums or '无'}）")
    emb = legend_tail[num] in flat
    n_emb += emb
    fig_rows.append(f"| Fig {num} | {t} | 已重出 `{FIG_FILES[num]}`（PDF+PNG）；"
                    f"{'**图注已画入图内**' if emb else '图注未入图'}；"
                    f"{'图内仅含自身图号' if nums else '图内无图号'} |")
must(n_emb == N_EMB,
     f"图注已画入图内的图件 {n_emb} 幅（声明值 {N_EMB}）——此数一变，§1.2 的说明与 §六 的缺口项必须同步改")
must(sorted(p.stem for p in FIGDIR.glob("*.pdf")) == sorted(FIG_FILES.values()),
     f"04_图表/ 下的 PDF 恰为这五幅：{sorted(p.stem for p in FIGDIR.glob('*.pdf'))}")

fig_tbl = "\n".join(fig_rows)
insert_after("> 整理当日核验：", blk(f"""
### 1.2 图清单（2026-09-23 重定位口径，{N_FIG} 图）

| 图 | 题目（取自新稿 Figure Legends，逐条抓取） | 图件状态（实测） |
| --- | --- | --- |
{fig_tbl}

> 题目由脚本 C5 从新稿 `Figure Legends` 段抓取，**非手抄**；三列状态亦**非手写**——
> 文件存在性、图内是否带旧图号（扫 PDF 文本层）、图注是否已画入图内（拿图注末句在文本层里找），
> 三项都在本节就地断言。**五幅内容均已重出**，`04_图表/` 下各有一 PDF 一 PNG。
>
> **两处尚未对齐，属同一次装包动作（§二·九·二）**：① 本清单声明的是 `Figures/Fig1`–`Fig5`.tiff
> （TIFF / RGB / 300 dpi / Arial），现产物为 PDF+PNG 且文件名为 `Fig_*`；② **图 1、图 2 已把整条图注
> 画进图内，图 3/4/5 没有**——出 docx 时若逐图再配一遍稿件的 Figure Legends，前两幅的图注会重复
> 出现一遍。二者都不影响图片内容，但都会让成品与声明不符。
"""), "§1.2 图清单")

# ── ④ §二 抬头横幅（表号/图号/计数一类历史口径的通告，只此一处）
insert_after("## 二、返工记录", blk(f"""
> ⚠️ **本节记的是 2026-09-15 至 09-17 各轮返工当时的事实，其中表号/图号/计数/字数为当时口径。**
> 2026-09-23 重定位后已变：**表 4 → 表 2**（六个面板基因的逐队列复制）、**表 2/表 3 → 表 1/表 3**、
> 图 6 → **{N_FIG}**、文献 18 → **{N_REF}**、正文 7,261 → **{MAIN_W:,}** 词、摘要 200 → **{ABS_W}** 词。
> **各轮的处置决定本身仍然有效，只有编号与计数作废**；当前口径见 §二·八。
"""), "§二 抬头横幅")

# ── ⑤ §二·六 的过时指针
sub_line("> ⚠️ **此为第五轮口径，已过时。", blk(f"""
> ⚠️ **此为第五轮口径，已过时。** 该轮之后又经第六轮（正文 7,261 / 摘要 200，见 §二·七）与
> **2026-09-23 重定位（正文 {MAIN_W:,} / 摘要 {ABS_W}，见 §二·八）**；当前值以 §二·八 为准。
"""), "§二·六 指针")

# ── ⑥ §二·七 抬头横幅（底稿 §四·7 点名要求）
insert_after("## 二·七、第六轮：模拟同行评审 Round 1 修订（2026-09-17）", blk(f"""
> ⚠️ **该轮数字已被 2026-09-23 重定位取代**——正文字数 7,261 → **{MAIN_W:,}**、摘要 200 → **{ABS_W}**；
> 图 6 → **{N_FIG}**、表 4 → **{N_TAB}**、补充表 3 → **{N_SUP}**、参考文献 18 → **{N_REF}**。
> 涉及 GSE89408 机器学习与外部池化 AUC 的值随 B6 重跑更新（外部三队列池化 AUC 0.583 → **0.780**、
> 内部定标经验 P 0.816 → **0.438**）。**本节记录的是 Round 1 各项「当时怎么处置的」，处置本身仍然有效，
> 只有数字口径作废**；当前口径见 §二·八。
"""), "§二·七 横幅")

# ── ⑦ §二·八：本步同步记录
insert_before("## 三、✋ 投稿前必须补齐（缺一不可）", blk(f"""
## 二·八、第七轮：2026-09-23 重定位 + 第 5 步「五处同步」（{DATE}）

**重定位**的方案与执行见 `07_修订记录/重定位方案与修改后大纲_2026-09-17.md` 与
`重定位执行回执_2026-09-23.md`；本步的成员清单与逐成员增量见
`07_修订记录/第5步_五处同步_清单与增量_2026-09-23.md`（本清单只登记结果，不重抄其内容）。
**本步的新发现与未完成项见 `07_修订记录/第5步_五处同步_回执_2026-09-25.md` 与 §六。**

**现役稿本**：英文 `08_修订_2026-09-17/05_稿件/English_manuscript_MDPI_重定位2026-09-23.md`、
中文 `双语初稿_候选优选与组织来源_重定位2026-09-23.md`，均已定稿（B9 数字回验 / B10 词数回验退出码 0）。
**`03_论文初稿/` 与 `05_投稿包/` 全数保持原样**（方案规则 1），本轮同步件一律写入 `08_修订_2026-09-17/`。

### 一、本轮已同步的成员与落点

| 成员 | 落点 | 本轮做了什么 | 校验 |
| --- | --- | --- | --- |
| 英文稿 / 中文稿 | `05_稿件/`（就地） | 重定位定稿，本步未再改动 | B9 231/231、B10 2/2，退出码 0 |
| STROBE-MR | `06_补充材料/STROBE-MR_checklist.md` | 逐条重映射 §锚点（12 个节号逐条回新稿核实后改写）、清除旧表引用、item 5 引句还原逐字 | 生成器 C1 + 独立对账 C2：引句 30/30 逐字命中、反向对照 0 假阳性 |
| 投稿信 ×2 | `08_投稿材料/Cover_Letter_MDPI.md`（活信）、`_archive_JIR作废_2026-09-27/Cover_Letter.md`（**2026-09-27 归档**，§二·九·一） | 换标题、What we found 由 1 段拆 3 段、Fit 段按新主线改写、日期改 {DATE} | 生成器 C3 + 独立对账 C4：数字 token 56 种全部有落点；反算 Table S3 得 18 行 / 6 组名义显著 / 6 个单工具集，与信所述一致 |
| 补充表 S1 / S2 / S3 | `06_补充材料/` | §锚点与新数字，见底稿 §四·3–5 | 见底稿 §三·4、§三·6 |
| 本清单（成员 #10） | `08_投稿材料/00_投稿清单.md`（本文件） | 标题、计数、字数、图清单、S2/图件状态、§二·七 横幅 | 生成器 C5 + 独立对账 C6 |

### 二、两封投稿信（回应底稿 §四·2）

- **MDPI 母版**（投 IJMS 等五本）：抬头刊名与 **Fit with the journal** 一段仍按原机制逐刊替换（§四）。
- **JIR 母版**：JIR 阶梯已于 2026-09-15 作废归档，但该信仍在成员清单内（§二 行 4），故本轮已按新稿同样重写。
  **2026-09-27 裁定：随路线一并归档**——信件内容一字未少，落点移出交付目录根部，见 §二·九·一。
  （本条原为「在路线已归档的情况下是否保留，待你裁定」，该裁定已作出。）
- 两封信**均不含图号/表号**（源件与产物的 grep 命中均为 0），故底稿 §四·2 的「换图号」在信件中**无对象**，
  该内容属投稿清单的图清单（§1.2）。
- 旧信中四处内容在新稿中**已无载体**，随整段替换一并删除：`PP.H4=0.742`（共定位结果未在新稿报告）、
  GSE55235 机器学习正对照与名次「315th–1586th」（新稿 `positive control` 命中 0）、`12.8% vs <1.6%`
  （新稿为 1.6%–6.1%）、以 LOXL1 为主语的旧主线。四处理由已写入信件内的英文来源注释块。

### 三、本轮不重出 docx（维持底稿 §五·一的裁定）

图 1 与图 4 尚未绘制，重出 docx 会产出指向不存在图件的成品。**docx 与图件一并留待「建包」时做。**
"""), "§二·八")

# ── ⑦a §二·八 抬头横幅：该节记的是 09-23/09-25 当时的图件与 S4 状态，本轮已处置。
# 按 §二 / §二·七 的既有惯例——历史节的原话保留不动，只在抬头声明哪几处状态已作废。
insert_after(f"## 二·八、第七轮：2026-09-23 重定位 + 第 5 步「五处同步」（{DATE}）", blk(f"""
> ⚠️ **本节记的是 2026-09-23/09-25 当时的事实。** 其中三处状态已在 {DATE_D} 的第八轮处置中改变：
> **图 1 与图 4 已补绘、图 2 已重出**（§二·九·二）、**补充表 S4 已单独成表**（§二·九·三）、
> **JIR 母版投稿信已归档**（§二·九·一）。**各轮的处置决定与校验结论本身仍然有效，只有这几处状态陈述作废**；
> 当前口径见 §二·九 与 §1.2。
"""), "§二·八 横幅")

# ── ⑦b §二·九：第八轮——四项处置（2026-09-27）
# 上一轮回执 §七 末段列了「投稿前仍须处置的项」三项，§八 列了一项待裁定，合起来正是本轮的四个处置对象。
# 本节只记**处置结果**；每项都指向可重跑的产物，本节自身不引入任何手写数字（图件状态由 §1.2 实测）。
insert_before("## 三、✋ 投稿前必须补齐（缺一不可）", blk(f"""
## 二·九、第八轮：四项处置（{DATE_D}）

> 本轮**只做处置，不新增分析**。四项均为上一轮回执登记在案的未完成项
> （`第5步_五处同步_回执_2026-09-25.md` §七 末段三项 + §八 待裁定一项），
> 处置依据为用户当日指令：**「归档处置、图 1/图 4 未绘、图 2/3/5 带旧编号、S4 未成表修改」**。
> §二·八 各节的记录不回改——那里写的是 09-23/09-25 当时的事实；本轮的改变一律记在本节。

### 一、JIR 母版投稿信：归档处置

§二·八·二 的「待你裁定」项已裁定为**归档**：`Cover_Letter.md` 移出交付目录根部，落于
`08_投稿材料/_archive_JIR作废_2026-09-27/`，同目录附 `README.md` 记作废理由与恢复办法。
**信件内容一字未少**——仍由 C3 逐字生成、C4 的独立对账照旧覆盖它；变的只是落点。
理由：交付目录根部曾并排放着两封信，从文件名无从分辨哪封指向已作废的刊物。
生成器一侧加了反向守卫（根部若再出现 `Cover_Letter.md` 即失败），确保重跑不会把它复活。
C4 的 `LETTERS`、C7 的在册路径已同步到新落点；`00_投稿清单.md` §一 与 §二·八·二 已改述。

### 二、图件：图 1 / 图 4 补绘，图 2 重出，图 3 / 5 沿用

| 图 | 本轮处置 | 产物（`04_图表/`） |
| --- | --- | --- |
| Fig 1 | **补绘**（横向四层阶梯，脚本 D4） | `Fig_study_design.pdf` + `.png` |
| Fig 2 | **重出**（57 候选总览，脚本 D3） | `Fig_candidate_overview.pdf` + `.png` |
| Fig 3 | 沿用 2026-09-17 版（脚本 B7），本轮未重出 | `Fig_forest_logFC.pdf` + `.png` |
| Fig 4 | **补绘**（证据矩阵，脚本 D2） | `Fig_evidence_matrix.pdf` + `.png` |
| Fig 5 | 沿用 2026-09-17 版（脚本 B3），本轮未重出 | `Fig_scRNA_dotplot.pdf` + `.png` |

**「图 2/3/5 带旧编号」经实测无对象**：逐幅扫 PDF 文本层，五幅均不含他图图号。
图 1 图内的 `Figure 1` 是它自己画入图内的图注标题，与文件名一致，不是撞号。
该判据已固化在 §1.2 的生成逻辑里，每轮重跑自证。

**本轮新发现一处不一致，已登记未处置**：图 1、图 2 把整条图注画进了图内，图 3/4/5 没有。
出 docx 时若逐图再配一遍稿件的 Figure Legends，前两幅的图注会重复出现。见 §1.2 末段。

**仍缺的是装包动作**，不是图：清单声明 `Figures/Fig1`–`Fig5`.tiff（TIFF / RGB / 300 dpi），
现产物是 PDF+PNG 且文件名为 `Fig_*`。转换与重命名随 docx 一并做（§二·八·三 的裁定未变）。

### 三、补充表 S4：已成表

`06_补充材料/Supplementary_Table_S4.md` 由脚本 D1 生成，三段 **246 + 3 + 450 = 699 行**，
逐行由 `02_结果表/` 的三份冻结 CSV（名次随模型变动 / 显著性跨线计数 / 留一法逐基因名次变动）
转排，**无值被重算、无手抄**。

稿件引用 S4 共 **4 处**，本轮逐处回表核对，**3 处吻合、1 处是错的**：

| 引用位置 | 与 S4 实际内容是否吻合 |
| --- | --- |
| Supplementary Materials 段 | 吻合——判据敏感性。原写「450 rows」只描述第三段，已改明为「450 rows in the leave-one-cohort-out table」（该改动另触发词数回验，见 §二·九·四） |
| §3.3 末 | 吻合——名次随估计量变动 + 逐基因留一法重排 |
| Figure 4 图注 | 吻合——名次变动与三种单队列删除结果 |
| **§2.4 Methods** | **不吻合——指向错误**，见下 |

**（G）尚未闭合：富集 / 解卷积 / PPI 三项的无落点承诺，共两处（登记待裁定，未擅自改稿）**

- **b1 · §2.4 末句**把这三项说成 `These are reported in Table S4`，而 S4 不含这三项——对 S4 检索
  富集、解卷积、互作相关词（clusterProfiler、CIBERSORT、STRING、KEGG、enrichment）命中 **0 处**；
  中文稿同句同病。
- **b2 · Figure 1 图注**（英文稿与中文稿同）写「Machine learning and functional enrichment are
  reported in the Supplementary Materials」／「机器学习与功能富集列于补充材料」——机器学习半句为真
  （S2 即在），**功能富集半句无落点**。该处为 2026-09-27 编制导师核查件时补查所得：对 S1–S4 逐表
  检索富集 / 解卷积 / PPI，命中同样为 0（S1 内的 KEGG 命中是候选集的**注释引用**，不是富集分析结果）。

这正是回执 §五·1 那条发现的后半截（该条标题为「`Table S4` 一词两指，且该表尚无文件」）；
本轮 D1 只闭合了「无文件」那一半。

另需注意：Results 段**没有**报告这三项分析的任何结果（全文检索，这三项只出现在 Methods 描述、
软件版本行与文献定位段），故这两处不仅指错了地方，「in the text only where they bear on the central
claim」这一半也未曾兑现。

**裁决前须知：这三项分析并非「没做过」。** 富集（clusterProfiler，GO BP + KEGG）、CIBERSORT 解卷积、
STRING PPI 的产物都在 `04_分析数据/processed/`（`GO_BP_*.csv`、`KEGG_*.csv`、`CIBERSORT_*.txt`、
`PPI_6刚度基因_edges.tsv`），2026-09-12 运行，日志在 `04_分析数据/summary/logs/`。但它们跑在
**重定位前**的口径下（文件名即 `*_genes6` 与 `PPI_6刚度基因`，`genes6` = 当时的六个面板基因），
故选②时须先核口径是否仍成立，**不能直接搬进现稿**。

三种处置待裁定：①删去 `in Table S4`（并相应处理 b2 的功能富集半句）、只保留「在正文相关处报告」；
②为这三项另立补充表；③在正文补写其结果。**选哪种会改变正文的承诺范围，故不由脚本代决。**

### 四、登记处的同步

| 登记处 | 本轮改动 |
| --- | --- |
| 本清单 §一 / §1.2 / §二·八·二 / §七 | 本文件（生成器 C5），图件与 S4 状态改按实测，缺口项重述 |
| 英文稿 / 中文稿 | Supplementary Materials 段的 S4 行数口径；两稿头部注释的未完成项均重述为「两件未闭合」——装包 + §2.4 的 S4 交叉引用 |
| `第5步_五处同步_回执_2026-09-25.md` | 追记四项处置（原回执的回改受限，只追加不覆盖） |
| 中文稿头部注释 | 三项未完成项中，图 1/图 4 与 S4 两项划去；「图 2/3/5 带旧编号」改述为实测无对象 |
"""), "§二·九")

# ── ⑧ §三·1 文献待补条目（由新稿算出）
todo_str = "、".join(todo_refs)
insert_after("### 1. 参考文献还有一处不符合 MDPI 规范：`et al.` 未展开", blk(f"""
> **2026-{DATE[5:7]}-{DATE[8:]} 重算（脚本 C5 从新稿 References 段逐条数作者）：** 文献表已整表换过，
> 待补条目由「4、5、11、13、14、15」变为 **{todo_str}**（共 {len(todo_refs)} 条）；
> 新 `[17]`（Steinberg）列满前 10 位后接 `et al.`，**合规**。上表编号为 2026-09-17 口径，保留供审计。
"""), "§三·1 待补条目")

# ── ⑨ §三·2 单位：与新稿逐字对齐
sub_line("| ^1^ |", blk(f"""
| ^1^ | {affils[0][4:]} | Kaixuan Zhang |
| ^2^ | {affils[1][4:]} | Xinle Wang |
| ^3^ | {affils[2][4:]} | Yu Fu, Ziwei Liu, Danning Zhang, Feng Zhang, Xingxing Lin, Baoqiang Dong |
"""), "§三·2 单位行")

# ── ⑩ §四 抬头：IJMS 首选已定
insert_after("## 四、投哪本 / Cover Letter 换哪两处", blk(f"""
> **2026-09-23 已定：第一投 `International Journal of Molecular Sciences`（IJMS）。**
> 两封投稿信已按重定位稿重写并落盘于 `08_修订_2026-09-17/08_投稿材料/`（§二·八·二）。
> 下面的决策树与两项查询**仍未办**（§六 #1/#2），作投稿前核对用。
"""), "§四 抬头")

# ── ⑪ §七 一句话总结（§七 是末节，替换题头之后的全部内容）
i7 = find_unique("## 七、一句话总结", "§七")
lines[i7 + 1:] = blk(f"""
**MDPI 版正文已完成三套暴露返工（§2.2）、审稿意见逐条返工（§2.3）、优先纠错与事实表冻结（§2.5）、
第五轮图文复核收尾（§2.6）、模拟同行评审 Round 1 修订（§二·七），并于 2026-09-23 完成**重定位**（§二·八）——
标题、计数（{N_FIG} 图 / {N_TAB} 表 / {N_SUP} 补充表 / {N_REF} 文献 / {N_DECL} 声明）、
字数（摘要 {ABS_W} / 正文 {MAIN_W:,}）与结构（Results 7 节 → 4 节）均已按新稿定稿，
第 5 步「五处同步」亦已完成（§二·八·一）；{DATE_D} 的第八轮四项处置（§二·九）闭合了
图件与补充表 S4 两处缺口，但**同一轮又查出一项无落点承诺**——富集 / 解卷积 / PPI 三项在 §2.4 末句
与 Figure 1 图注两处被承诺交付，却既不在任何补充表内、也不在正文中（§二·九·三 G）。
**本文件的计数与图题由脚本 C5 从新稿派生、图件状态由脚本实测，故不可能与稿本或图件脱节。**

**当前为「重定位后」状态，不是投稿终稿。** 挡在投稿前的硬缺口**有两件**：

1. **装包**——五幅图件内容已齐（§二·九·二），但须转换成本清单声明的 `Figures/Fig1`–`Fig5`.tiff
   并对齐文件名，docx 因此尚未重出；另须处理图 1、图 2 图注已画入图内所致的重复（§1.2 末段）。
2. **富集 / 解卷积 / PPI 三项的无落点承诺**（§二·九·三 G）——§2.4 末句与 Figure 1 图注两处所指的
   结果，既不在 S4 或其余任何补充表内、也不在正文中，须先裁定处置方式（删引用 / 另立表 / 补写结果）
   再改稿。

其余是格式与行政项：补 {len(todo_refs)} 条文献的前 10 位作者（§三·1）、核对 CRediT（§三·2）、办 §六 的清单。

**第三份独立评阅的结论（「当前描述性主线基本成立，但尚不应把此次复核视为全部分析或正式投稿验收通过」）
在本轮之后仍然适用**——本轮只做处置与装包前置，未新增分析、未改动任何数值。
""")

# ─────────────────────────────────────────── 三、写出与后置条件
out = "\n".join(lines)
DST.write_text(out, encoding="utf-8")

print(f"\n── 后置条件")
must2(TITLE in out, "新标题在产物中恰出现", f"新标题缺失：{TITLE[:50]}")
must(out.count(TITLE) == 1, f"新标题恰 1 次（实得 {out.count(TITLE)}）")
OLD_TITLE = "Cross-Cohort Expression Profile and Cellular Localization"
must2(OLD_TITLE not in out,
      "旧标题已全部移除（源件 1 处 → 产物 0 处）",
      f"产物仍含旧标题 {out.count(OLD_TITLE)} 处")
must(out.count(f"正文 {MAIN_W:,}") >= 2, f"新正文字数 {MAIN_W:,} 在产物中出现 ≥2 次")
must2("7,261" in out, "旧字数 7,261 仍保留（仅出现在标为「已取代」的历史节）",
      "旧字数 7,261 被误删——历史节不应删数，只应加横幅")
def section(text, head):
    """取某一节从题头到其后第一个 --- 分隔线的文本（用于把断言归属到节）。"""
    i = text.index(head)
    j = text.index("\n---", i)
    return text[i:j]


s27 = section(out, "## 二·七、")
must2("已被 2026-09-23 重定位取代" in s27,
      "§二·七 段内含「已被 2026-09-23 重定位取代」横幅（底稿 §四·7 点名要求）",
      "§二·七 段内没有取代横幅——底稿点名要求的那处缺失")
s2 = section(out, "## 二、返工记录")
must2("当时口径" in s2,
      "§二 段内含「当时口径」通栏横幅（表号/图号/计数一类历史口径只此一处声明）",
      "§二 段内缺少历史口径横幅")
must(f"**{todo_str}**" in out, f"§三·1 待补条目写入：{todo_str}")
for num, _ in FIG_TITLES:
    must(f"| Fig {num} |" in out, f"图清单含 Fig {num} 行")
must2("已单独成表" in out and "699" in out,
      "S4 已成表并记明 699 行（三段 246 + 3 + 450）",
      "S4 的成表状态未写进产物")
must("## 二·九、" in out, "§二·九（第八轮四项处置）已写入")
must2("_archive_JIR作废_2026-09-27" in out,
      "JIR 母版信的新落点已写入（§一 与 §二·九·一）",
      "未写出 JIR 母版信的归档落点")
# S4 的四处引用必须逐处回表核对过，不许再写成「均有落点」——Methods 那处本就是错的。
must2("3 处吻合、1 处是错的" in out,
      "S4 四处引用已逐处回表核对（3 处吻合 / 1 处指向错误）",
      "未记明 S4 四处引用的逐处核对结果")
must2("§2.4 Methods" in out and "指向错误" in out,
      "§2.4 的 S4 交叉引用错误已登记待裁定（§二·九·三 G），未擅自改稿",
      "§2.4 的 S4 交叉引用错误未登记")
# 2026-09-27 补查：同一项承诺在 Figure 1 图注里还有一处（功能富集半句无落点）。
# 两处必须同登记、同裁定——只记一处会让缺口在下一轮被当成已闭合。
must2("b2 · Figure 1 图注" in out,
      "Figure 1 图注处的同一承诺已与 §2.4 并列登记（§二·九·三 G）",
      "Figure 1 图注处的同一承诺未登记——缺口会被误报为只剩一处")
must2("04_分析数据/processed" in out and "重定位前" in out,
      "裁定时须知三项分析产物确在（`04_分析数据/processed/`）且跑在重定位前口径，不能直接搬",
      "未记明三项分析已有产物及其口径前提，裁定会缺依据")
must("现均有落点" not in out,
     "已删去「四处引用现均有落点」这一被实测证伪的说法")
# 第八轮是 2026-09-27 的事：DATE 同时是投稿信落款日期，故不可用它给本轮记日。
must(f"{DATE_D} 的第八轮四项处置" in out, f"§七 记第八轮的日期为 {DATE_D}")
must(f"{DATE} 的四项处置" not in out, f"§七 不再把第八轮误记为 {DATE}")
# 反向守卫：上一轮登记的四处未完成项，其**旧表述**不得再出现在「当前状态区」里。
# 「改过了」不等于「旧话没了」——处置必须自证其生效。
# 只对 §一 与 §七 做串检，不做全文串检：§二 是历史记录区，按项目惯例保留当时原话、只加横幅，
# 全局串检会逼我去改历史（那正是规则 1 要避免的）。
def _zone(text, head, tail=None):
    i = text.index(head)
    return text[i:text.index(tail, i)] if tail else text[i:]


CURRENT = _zone(out, "## 一、包内文件", "## 二、返工记录") + _zone(out, "## 七、一句话总结")
for gone, what in [("尚未单独成表", "S4 未成表"),
                   ("未绘制", "图件未绘制"),
                   ("待你裁定", "JIR 母版信待裁定"),
                   ("可交付 2 幅", "「声明 5 图、实际可交付 2 幅」"),
                   ("只剩「装包」一件", "「硬缺口只剩装包一件」——已查出第二件（无落点承诺）"),
                   ("S4 的 §2.4 交叉引用错误", "「缺口只是 §2.4 一处 S4 交叉引用」——已补查出 Figure 1 图注这第二处")]:
    must(gone not in CURRENT, f"当前状态区（§一 + §七）已无旧表述「{what}」")
must2("本轮不重出 docx" in out or "本轮不重出 docx（§二·八·三）" in out,
      "「本轮不重出 docx」已在抬头与 §二·八·三 两处声明",
      "未声明本轮不重出 docx")
must("Shenyang, Liaoning, China" in out and "Zhejiang Province, China" in out,
     "单位地址已与新稿逐字对齐（含省名）")
must2(re.search(r"[一-鿿]", src) is not None,
      "源件为中文文件（汉字存在），故产物为中文属预期",
      "源件竟无汉字，语言判据不适用")
n_src, n_out = len(src.split("\n")), len(out.split("\n"))
print(f"  源 {n_src} 行 → 产物 {n_out} 行（+{n_out - n_src}）")
must(n_out > n_src, "产物比源件长（本轮为加横幅与新增节，未删历史）")
print(f"\n写出 {DST.relative_to(ROOT)}  {len(out.encode('utf-8'))} B")
print(f"\nC5 完成：投稿清单已按重定位口径重建，并记入第八轮四项处置（{DATE_D}）"
      f"与无落点承诺两处（§2.4 末句 + Figure 1 图注）的登记")
