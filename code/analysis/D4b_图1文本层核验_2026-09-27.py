#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
D4b — 图 1 的 **PDF 侧**核验（D4 只查了 PDF 魔数，那只证明文件是 PDF，不证明画对了）

为什么必须有这一步：本项目有过「PNG 正常、PDF 侧静默丢字」的先例（基线 pdf 设备遇非 ASCII
触发 mbcsToSbcs conversion failure，ggsave 照常返回成功），也有过「图文件内自带旧编号而
文本层校验看不到」的先例。故画完图必须**分别**核 PDF 文本层，不能只核 PNG。

四类判据：
  ① 片段数对账：逐箱数 PDF 文本层的片段，须等于 D4 写进 图1_研究设计_箱体.csv 的 n_text。
     这是唯一能发现「某行被裁掉 / 漏画」的判据——图上少一行，肉眼多半看不出来。
  ② 归属唯一：每个片段的起始 x 必须落在**恰好一个**箱子的 x 区间内（图注与括注除外）。
     箱子是等宽等距算出来的，若排版错位，归属会落到箱外或落到两个箱里。
  ③ 阶序：四个编号层的标题在 x 上严格左→右递增，且 y 严格递减（阶梯的序）。
  ④ 越界与字符：所有片段在 MediaBox 内；不含非 ASCII / 替换符（U+FFFD）。

输入：04_图表/Fig_study_design.pdf、02_结果表/图1_研究设计_箱体.csv
"""
import csv
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from pdftext import pdf_text  # noqa: E402

REV = HERE.parents[1]
PDF = REV / "04_图表" / "Fig_study_design.pdf"
BOX = REV / "02_结果表" / "图1_研究设计_箱体.csv"

FAIL = 0


def must(cond, msg):
    global FAIL
    if cond:
        print(f"  OK   {msg}")
    else:
        FAIL += 1
        print(f"  FAIL {msg}")


print("D4b 图 1 文本层核验")
recs, raw = pdf_text(PDF)
must(len(recs) > 0, f"PDF 文本层解出 {len(recs)} 个片段（0 个即解析器或文件有问题，不是「没画字」）")

m = re.search(rb"/MediaBox\s*\[([^\]]*)\]", raw)
must(m is not None, "可读出 /MediaBox")
mb = [float(v) for v in m.group(1).split()]
PW, PH = mb[2] - mb[0], mb[3] - mb[1]
print(f"  页面 {PW:.0f} x {PH:.0f} pt（{PW / 72:.2f} x {PH / 72:.2f} in）")

boxes = list(csv.DictReader(open(BOX, encoding="utf-8")))
must(len(boxes) == 5, f"箱体审计表 {len(boxes)} 行（期望 5）")
for b in boxes:
    for k in ("x_left", "x_right", "y_top", "y_bottom", "n_text"):
        b[k] = float(b[k])

# ── ① 片段数对账
# 归属必须**同时**看 x 与 y：图注左端与箱 1 的左缘 x 完全相同（都在 MARGIN 上），
# 只看 x 会把 8 行图注算进箱 1（首轮实测即如此：箱 1 报 24，实为 16 + 8 行图注）。
TOL_X, TOL_Y = 1.0, 3.0  # pt
per = {i: [] for i in range(len(boxes))}
outside = []
for y, x, t in recs:
    hit = [i for i, b in enumerate(boxes)
           if b["x_left"] * 72 - TOL_X <= x <= b["x_right"] * 72 + TOL_X
           and b["y_bottom"] * 72 - TOL_Y <= y <= b["y_top"] * 72 + TOL_Y]
    if len(hit) == 1:
        per[hit[0]].append((y, x, t))
    elif len(hit) == 0:
        outside.append((y, x, t))
    else:
        must(False, f"片段 (x={x}, y={y}) 同时落进 {len(hit)} 个箱：{t!r}")

for i, b in enumerate(boxes):
    must(len(per[i]) == int(b["n_text"]),
         f"箱 {i + 1} [{b['head']}] PDF 片段 {len(per[i])} = 审计表 n_text {int(b['n_text'])}")

# ── ② 箱外只能有两样东西：通栏括注（D4 落盘的括注文件）与图注（D4 落盘的图注文件）。
#    两样都逐字比对——本项目的教训是「图注的结构描述会随图改而失真」，逐值核对读不到散文。
norm = lambda s: re.sub(r"\s+", " ", s).strip()
want_note = norm((REV / "02_结果表" / "图1_研究设计_括注.txt").read_text(encoding="utf-8"))
want_cap = norm((REV / "02_结果表" / "图1_研究设计_图注.txt").read_text(encoding="utf-8"))

note = [o for o in outside if norm(o[2]) == want_note]
capf = [o for o in outside if o not in note]
must(len(note) == 1, f"通栏括注在 PDF 文本层中逐字相符（{len(note)} 段）")
got_cap = norm(" ".join(t for _, _, t in sorted(capf, key=lambda r: -r[0])))
must(got_cap == want_cap,
     f"图注 {len(capf)} 行拼合后与 D4 落盘的图注逐字相符（{len(got_cap)} 字）")
if got_cap != want_cap:
    for i, (a, b) in enumerate(zip(got_cap.split(), want_cap.split())):
        if a != b:
            print(f"       首个差异在第 {i + 1} 个词：PDF={a!r} 期望={b!r}")
            break

# 箱框线上不该有文字起点（左边距对齐的正文起点在 x_left + PAD，不是 x_left）
edge = [(y, x, t) for i in range(len(boxes)) for (y, x, t) in per[i]
        if abs(x - boxes[i]["x_left"] * 72) < 0.5 or abs(x - boxes[i]["x_right"] * 72) < 0.5]
must(not edge, f"没有片段的起始 x 压在箱框线上（{len(edge)} 个）")

# ── ③ 阶序
heads = {}
for y, x, t in recs:
    mm = re.match(r"^([1-4])\s\s", t)
    if mm:
        heads[int(mm.group(1))] = (y, x, t)
must(sorted(heads) == [1, 2, 3, 4], f"四个编号层标题齐备：{sorted(heads)}")
xs = [heads[k][1] for k in sorted(heads)]
ys = [heads[k][0] for k in sorted(heads)]
must(all(b > a for a, b in zip(xs, xs[1:])), f"层标题 x 严格递增：{xs}")
must(all(b < a for a, b in zip(ys, ys[1:])), f"层标题 y 严格递减（阶梯的序）：{ys}")
must(ys[0] < boxes[0]["y_top"] * 72, "数据层箱顶在层 1 之上（数据层是最高一级）")

# ── ④ 越界与字符
bad = [(y, x, t) for y, x, t in recs
       if x < -TOL_X or x > PW + TOL_X or y < -TOL_Y or y > PH + TOL_Y]
must(not bad, f"没有片段越出页面（{len(bad)} 个）")
na = [(y, x, t) for y, x, t in recs if any(ord(c) > 126 or ord(c) < 32 for c in t)]
must(not na, f"文本层无非 ASCII / 控制符（{len(na)} 个）——非 ASCII 在基线 pdf 设备上会静默丢字")
rep = [t for _, _, t in recs if "�" in t]
must(not rep, f"文本层无替换符 U+FFFD（{len(rep)} 个）")

# 全图文本拼合后，图注末句必须在位（这句是阶梯纵向位移的声明，缺了就等于加了未声明的通道）
flat = " ".join(t for _, _, t in sorted(recs, key=lambda r: (-r[0], r[1])))
flat = re.sub(r"\s+", " ", flat)
tail = "The four layers are ordered by how strongly each bounds the priority claim, not by the order in which the analyses were run."
must(tail in flat, "图注末句（阶梯序的声明）在 PDF 文本层中完整在位")

print(f"\n结果：{'全部通过' if not FAIL else f'{FAIL} 项未通过'}")
raise SystemExit(1 if FAIL else 0)
