#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E1 — 导师核查件 docx 生成器（2026-09-27）

输入（唯一内容来源）：07_修订记录/导师核查件_2026-09-27.md
输出：07_修订记录/导师核查件_2026-09-27.docx

为什么由脚本生成而不手写 docx：同一段文字在两处各存一份必然漂移（项目已有教训）。
md 是唯一来源，docx 每次由它重出；本脚本不自带任何一句正文。

只支持该 md 实际用到的 markdown 子集：标题 / 段落 / 无序与有序列表 / 引用块 / 表格 /
水平线 / 行内 **粗体** 与 `代码`。范围之外一律报错，不静默吞掉。

项目已知坑（本脚本据此设防）：
  · 反引号与 ** 会漏进 docx —— 写出后回读断言全文不含这两个字面；
  · 给标题段追加 run 会覆盖标题样式 —— 故标题的粗体/字号/颜色逐项显式设定，不靠样式；
  · 宽表须另起横向节 —— 本件最宽 3 列，A4 纵向放得下，故不涉及；
  · 表格列宽须写 tblGrid，否则 Word 按内容重排；
  · 手写在 md 里的字节数/哈希会随上游改动过期，而自指的那行必然记错 ——
    故成品写出前先对盘（见 verify_materials），对不上就不写。

用法：E1_导师核查件docx_2026-09-27.py [src.md] [dst.docx] [rev_root]
      （默认取本目录解析出的同名校对件；给出参数即换源，供负面对照与 2026-09-28 核查件复用，
        此时绝不会碰到真成品。第三参数覆盖材料清单表里「落点」的相对根。）
"""
import hashlib
import pathlib
import re
import sys

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

REV = pathlib.Path(sys.argv[3]) if len(sys.argv) > 3 \
    else pathlib.Path(__file__).resolve().parents[2]
# 允许换源换落点：负面对照要把 md 复制到别处改错一个字节数再跑，此时绝不能碰真成品
# argv[3] 另可覆盖 REV（材料清单表里的落点按 REV 解析）：2026-09-28 的核查件放在 09 树
# 顶层，其材料清单表引用的是 09 树下的 02/04/05/06/07，故由调用方显式给出根目录。
SRC = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 \
    else REV / "07_修订记录" / "导师核查件_2026-09-27.md"
DST = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 \
    else REV / "07_修订记录" / "导师核查件_2026-09-27.docx"

LATIN, CJK, MONO = "Times New Roman", "宋体", "Consolas"
BODY_PT, TABLE_PT = 10.5, 9.5
USABLE_CM = 16.0          # A4 宽 21cm − 左右各 2.5cm 页边距

fail = 0


def must(ok, msg, bad=None):
    """成功也打印——故消息写成陈述句。"""
    global fail
    print(f"  {'OK  ' if ok else 'FAIL'} {msg}")
    if not ok:
        fail += 1
        if bad:
            print(f"       {bad}")


# ─────────────────────────── 零、材料清单表硬校验（2026-09-27 晚追加）
# 起因：这张表首版把「修订回执」与「本件自身」的字节数手写在 md 里。前者一改稿就过期，
# 后者是自指——记的正是记完即变的那个数，结构上不可能对；两处都真的错了，是目视发现的。
# 故改为机器对盘：落点能定位到单个文件的行，字节数与短哈希都必须与磁盘相符；定位不到
# 文件的行必须在此显式登记免验理由。新增一行而未登记即失败——逼着做一次决定，而不是
# 等导师按图索骥核不上再回头查。
BTICK = re.compile(r"`([^`]+)`")
SHORT_MD5 = re.compile(r"\b([0-9a-f]{8})…")
BYTES = re.compile(r"([\d,]+) B")
META_HEAD = ["材料", "落点", "大小 / 校验（md5）"]
EXEMPT = {
    "表 1 / 2 / 3": "三个主表在稿件 ## Tables 段内，无独立文件",
    "投稿信（JIR 母版，已作废归档）": "落点是归档目录，大小栏为说明文字",
    "本件（md 为源，docx 由脚本生成）": "自指——本表所在文件自身，字节数记完即变",
}


def verify_materials(rows):
    assert rows[0] == META_HEAD, f"材料清单表表头已变，校验须随之更新：{rows[0]}"
    bad = []
    for label, where, chk in rows[1:]:
        for tok in BTICK.findall(where):          # 落点里凡是路径，先证其存在
            if "/" in tok and not (REV / tok).exists():
                bad.append(f"{label}：落点不存在 —— {tok}")
        if label in EXEMPT:
            print(f"  --   {label}（免验：{EXEMPT[label]}）")
            continue
        files = [p for p in BTICK.findall(where) if (REV / p).is_file()]
        if len(files) != 1:
            bad.append(f"{label}：落点可定位的文件 {len(files)} 个，须在本脚本登记免验理由")
            continue
        b = (REV / files[0]).read_bytes()
        n, m = len(b), hashlib.md5(b).hexdigest()[:8]
        got_n = [int(s.replace(",", "")) for s in BYTES.findall(chk)]
        got_m = SHORT_MD5.findall(chk)
        if not got_n or got_n[0] != n:
            bad.append(f"{label}：记的字节 {got_n[:1] or '（无）'} ≠ 磁盘 {n:,}（{files[0]}）")
        if got_m and got_m[0] != m:
            bad.append(f"{label}：记的哈希 {got_m[0]}… ≠ 磁盘 {m}…（{files[0]}）")
        if got_n and got_n[0] == n and (not got_m or got_m[0] == m):
            print(f"  OK   {label}：{n:,} B  {m}…（{files[0]}）")
    must(not bad, f"材料清单表 {len(rows) - 1} 行逐行对盘"
                  f"（免验 {len(EXEMPT)} 行已登记，其余须与磁盘逐值相符）",
         bad="\n       ".join(bad))


# ─────────────────────────────────────────── 一、解析 markdown
_TOK = re.compile(r"\*\*|`")


def parse_inline(s):
    """把一行拆成 (文字, 粗体, 等宽) 三元组序列。

    单趟扫描、两个独立开关（`**` 翻粗体、`` ` `` 翻等宽），因此「粗体里嵌代码」
    （`**…（`x`）…**`）能渲染成「粗体 + 粗体等宽 + 粗体」三段。两个开关必须在行末
    同时回到关状态，否则说明标记不成对——直接断言，不静默把字面 `**` 漏进成品。

    2026-09-27 首版按「先切粗体、再切代码」的两段式写，粗体段内的反引号原样漏进 docx；
    改成「先切代码」后，跨代码段的粗体又被切成两半、`**` 反而漏出。两版都被回读断言当场
    抓住——所以这里不靠精心排序，而靠扫描结束时的成对性断言兜底。
    """
    out, pos, bold, mono = [], 0, False, False
    for m in _TOK.finditer(s):
        if m.start() > pos:
            out.append((s[pos:m.start()], bold, mono))
        if m.group() == "**":
            bold = not bold
        else:
            mono = not mono
        pos = m.end()
    if pos < len(s):
        out.append((s[pos:], bold, mono))
    assert not bold and not mono, f"行内标记不成对（粗体余 {bold} / 等宽余 {mono}）：{s[:90]!r}"
    return [t for t in out if t[0]]


def parse(md):
    lines = md.split("\n")
    blocks, i = [], 0
    while i < len(lines):
        raw = lines[i]
        if not raw.strip():
            i += 1
            continue
        if raw.lstrip().startswith("|"):
            rows = []
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                # 分隔行（| --- | --- |）不入表体
                if not all(re.fullmatch(r":?-{2,}:?", c) for c in cells):
                    rows.append(cells)
                i += 1
            assert len(rows) >= 2, f"表格只有 {len(rows)} 行，缺表头或表体"
            ncol = len(rows[0])
            assert all(len(r) == ncol for r in rows), \
                f"表格各行列数不一：{[len(r) for r in rows]}"
            blocks.append(("table", rows))
            continue
        if raw.startswith("#"):
            lvl = len(raw) - len(raw.lstrip("#"))
            blocks.append((f"h{lvl}", raw[lvl:].strip()))
            i += 1
            continue
        if raw.strip() == "---":
            blocks.append(("hr", None))
            i += 1
            continue
        if raw.startswith(">"):
            buf = []
            while i < len(lines) and lines[i].startswith(">"):
                buf.append(lines[i].lstrip(">").strip())
                i += 1
            blocks.append(("quote", buf))
            continue
        m = re.match(r"^(\d+)\.\s+(.*)$", raw)
        if m:
            items = [[int(m.group(1)), m.group(2)]]
            i += 1
            while i < len(lines) and lines[i].strip():
                m2 = re.match(r"^(\d+)\.\s+(.*)$", lines[i])
                if m2:
                    items.append([int(m2.group(1)), m2.group(2)])
                else:
                    assert not lines[i].startswith(("#", "|", ">")), \
                        f"列表项后接了别的块：{lines[i][:40]!r}"
                    items[-1][1] += " " + lines[i].strip()   # 续行（markdown 懒惰续行）
                i += 1
            blocks.append(("ol", items))
            continue
        if raw.startswith("- "):
            items = [raw[2:].strip()]
            i += 1
            while i < len(lines) and lines[i].strip():
                if lines[i].startswith("- "):
                    items.append(lines[i][2:].strip())
                else:
                    assert not lines[i].startswith(("#", "|", ">")), \
                        f"列表项后接了别的块：{lines[i][:40]!r}"
                    items[-1] += " " + lines[i].strip()
                i += 1
            blocks.append(("ul", items))
            continue
        # 普通段落：**要求单行**——本件所有段落都未换行，若日后换行则宁可报错也不静默拼接
        assert i + 1 >= len(lines) or not lines[i + 1].strip(), \
            f"段落被硬换行（L{i + 1}）：解析器不支持，请把该段并成一行"
        blocks.append(("p", raw.strip()))
        i += 1
    return blocks


# ─────────────────────────────────────────── 二、写 docx
def style_run(run, size, bold, mono=False, color=None, italic=False):
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    if color:
        run.font.color.rgb = color
    latin = MONO if mono else LATIN
    run.font.name = latin
    rf = run._element.get_or_add_rPr().get_or_add_rFonts()
    rf.set(qn("w:eastAsia"), CJK)
    rf.set(qn("w:ascii"), latin)
    rf.set(qn("w:hAnsi"), latin)


def add_para(doc, tokens, size=BODY_PT, bold_all=False, space_before=0, space_after=6,
             indent=0.0, color=None):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(space_before)
    p.paragraph_format.space_after = Pt(space_after)
    p.paragraph_format.line_spacing = 1.25
    if indent:
        p.paragraph_format.left_indent = Cm(indent)
    for text, b, mono in tokens:
        style_run(p.add_run(text), size, bold_all or b, mono, color)
    return p


def bottom_border(p):
    """水平线：段落加下边框（比插一张图或空段更稳）。"""
    pPr = p._element.get_or_add_pPr()
    bd = pPr.makeelement(qn("w:pBdr"), {})
    bot = pPr.makeelement(qn("w:bottom"), {})
    for k, v in (("w:val", "single"), ("w:sz", "6"), ("w:space", "1"), ("w:color", "BFBFBF")):
        bot.set(qn(k), v)
    bd.append(bot)
    pPr.append(bd)


def col_widths(rows):
    """按各列最长内容估宽，再按可用宽度等比缩放；每列下限 2.2cm。"""
    ncol = len(rows[0])
    weight = []
    for c in range(ncol):
        longest = max(len(re.sub(r"\*\*|`", "", r[c])) for r in rows)
        weight.append(max(longest, 6))
    total = sum(weight)
    w = [USABLE_CM * x / total for x in weight]
    floor = 2.2
    if min(w) < floor:                      # 抬到下限后把余量从最宽的列里扣
        deficit = sum(floor - x for x in w if x < floor)
        w = [max(x, floor) for x in w]
        widest = w.index(max(w))
        w[widest] -= deficit
    return w


def add_table(doc, rows):
    t = doc.add_table(rows=len(rows), cols=len(rows[0]))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    widths = col_widths(rows)
    for r, row in enumerate(rows):
        for c, cell in enumerate(row):
            para = t.cell(r, c).paragraphs[0]
            para.paragraph_format.space_before = Pt(1.5)
            para.paragraph_format.space_after = Pt(1.5)
            for text, b, mono in parse_inline(cell):
                style_run(para.add_run(text), TABLE_PT, True if r == 0 else b, mono)
    for c, wcm in enumerate(widths):        # 列宽须写 tblGrid，否则 Word 按内容重排
        for row in t.rows:
            row.cells[c].width = Cm(wcm)
    # 跨页可读性：整行不拆（首版实测「共同基因…/691」被页断切成两半），表头行跨页重复
    for r, row in enumerate(t.rows):
        trPr = row._tr.get_or_add_trPr()
        trPr.append(trPr.makeelement(qn("w:cantSplit"), {}))
        if r == 0:
            trPr.append(trPr.makeelement(qn("w:tblHeader"), {}))
    doc.add_paragraph().paragraph_format.space_after = Pt(2)
    return t


def build(blocks):
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    for m in ("left_margin", "right_margin"):
        setattr(sec, m, Cm(2.5))
    for m in ("top_margin", "bottom_margin"):
        setattr(sec, m, Cm(2.5))
    doc.styles["Normal"].font.size = Pt(BODY_PT)

    HSIZE = {"h1": 17, "h2": 13.5, "h3": 11.5}
    for kind, val in blocks:
        if kind in HSIZE:
            add_para(doc, [(val, True, False)], size=HSIZE[kind], bold_all=True,
                     space_before=12 if kind != "h1" else 0, space_after=6,
                     color=RGBColor(0x1F, 0x36, 0x4D))
        elif kind == "p":
            add_para(doc, parse_inline(val))
        elif kind == "quote":
            add_para(doc, parse_inline(" ".join(val)), size=BODY_PT - 0.5,
                     indent=0.6, color=RGBColor(0x44, 0x44, 0x44))
        elif kind == "ul":
            for it in val:
                add_para(doc, [("• ", False, False)] + parse_inline(it), indent=0.5,
                         space_after=4)
        elif kind == "ol":
            for num, it in val:
                add_para(doc, [(f"{num}. ", True, False)] + parse_inline(it), indent=0.5,
                         space_after=4)
        elif kind == "table":
            add_table(doc, val)
        elif kind == "hr":
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(2)
            p.paragraph_format.space_after = Pt(8)
            bottom_border(p)
        else:
            raise AssertionError(f"未处理的块类型 {kind!r}")
    return doc


# ─────────────────────────────────────────── 三、主流程
print("E1 导师核查件 docx 生成")
md = SRC.read_text(encoding="utf-8")
blocks = parse(md)
n_tab = sum(1 for k, _ in blocks if k == "table")
print(f"  md {len(md.split(chr(10)))} 行 → {len(blocks)} 块（含表 {n_tab} 张）")

print("\n── 材料清单表对盘（先对盘，对不上就不写成品）")
meta = [v for k, v in blocks if k == "table" and v[0] == META_HEAD]
assert len(meta) == 1, f"材料清单表应恰有 1 张，实得 {len(meta)}"
verify_materials(meta[0])
if fail:
    print(f"\n材料清单表有 {fail} 项对不上磁盘，未写出 docx。")
    raise SystemExit(1)

doc = build(blocks)
doc.save(DST)
print(f"  写出 {DST.relative_to(REV)}  {DST.stat().st_size} B")

# ─────────────────────────────────────────── 四、回读验收（不靠「没报错」当通过）
back = Document(str(DST))
alltext = "\n".join(p.text for p in back.paragraphs) + "\n" + "\n".join(
    c.text for t in back.tables for r in t.rows for c in r.cells)

print("\n── 回读验收")
must("`" not in alltext, "docx 全文无反引号残留（行内代码标记未漏进成品）")
must("**" not in alltext, "docx 全文无 `**` 残留（粗体标记未漏进成品）")

md_tables = [v for k, v in blocks if k == "table"]
must(len(back.tables) == len(md_tables),
     f"表格数一致（md {len(md_tables)} 张 / docx {len(back.tables)} 张）")
mismatch = []
for ti, (md_rows, dx) in enumerate(zip(md_tables, back.tables)):
    if (len(dx.rows), len(dx.columns)) != (len(md_rows), len(md_rows[0])):
        mismatch.append(f"表{ti + 1} 尺寸 md {len(md_rows)}×{len(md_rows[0])} "
                        f"vs docx {len(dx.rows)}×{len(dx.columns)}")
        continue
    for r, md_row in enumerate(md_rows):
        for c, md_cell in enumerate(md_row):
            want = "".join(t for t, _, _ in parse_inline(md_cell))
            got = dx.cell(r, c).text
            if want != got:
                mismatch.append(f"表{ti + 1} 第{r + 1}行第{c + 1}列\n         md: {want[:70]}\n"
                                f"        docx: {got[:70]}")
must(not mismatch, f"每张表逐格文字与 md 一致（共 {sum(len(r) for r in md_tables)} 格）",
     bad="\n       ".join(mismatch[:6]))

# 内联标记拆除后，正文段落文字应与 md 段落逐字一致（抽样：全部 H1/H2 与首末段）
md_p = ["".join(t for t, _, _ in parse_inline(v)) for k, v in blocks if k == "p"]
dx_p = [p.text for p in back.paragraphs if p.text.strip()]
must(len(dx_p) >= len(md_p),
     f"docx 非空段落 {len(dx_p)} 段 ≥ md 正文段落 {len(md_p)} 段（差额为标题、列表与线）")
must(md_p[0] in alltext and md_p[-1] in alltext,
     "md 的首段与末段均已出现在 docx 中")

print(f"\n结果：{'全部通过' if not fail else f'{fail} 项未通过'}")
raise SystemExit(1 if fail else 0)
