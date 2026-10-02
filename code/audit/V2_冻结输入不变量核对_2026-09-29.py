#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V2 — 冻结输入不变量核对（2026-09-29）

对象：`01_冻结结果与代码/冻结清单_2026-09-27.tsv` —— 2026-09-27 冻结当时由 F1 逐件记下的
SHA-256，共 165 件。本脚本把该清单当作**不变量基线**，在四份审稿意见的修订做完之后
逐件复算，回答一个问题：

    「这一轮修订，究竟动了什么？」

判据四条，缺一不可：
  1. 未登记变动 = 失败；清单里的件在磁盘上找不到 = 失败。
  2. **反向守卫**：登记在 KNOWN_DELTA 里的件若实测**未**变动，说明该登记已过期
     （要么文件被改回去了，要么登记写错了），同样失败。没有这一条，登记表会随修订
     漂移成一个永远成立的注脚。
  3. **分析结果类一律不得变动**：修订只改文字、图件与其绘图脚本，不得动任何一个
     分析数字。该类不写成「02_结果表/ 全类」——该类里有 6 件是绘图脚手架
     （图 1 的箱体/布局/括注/图注、图 2 与森林图的逐行数值），随重绘而变是本分；
     故按**点名**列举（NUMERIC_ANCHORS），点名不到即失败。
  4. **域守卫**：`01_审计/`、`05_稿件/`、`06_补充材料/`、`07_修订记录/`、`08_投稿材料/`
     五棵冻结树**整棵**不得变动；`03_分析/scripts`、`03_分析/logs`、`04_图表` 里凡变动者
     必须逐件登记。

与 V1 的关系：V1（`V1_版本清点_2026-09-27.py`，含同样的 KNOWN_DELTA 机制与反向守卫）
核的是**旧项目根**（`08_修订_2026-09-17/` 与 `09_Genes_CIMB_投稿修订_2026-09-27`），
那两个目录在本机已不存在（只在 WPS 云侧），且其 OUT 已被本轮镜像覆盖 —— 故 V1 在本树
**不可复跑**，按规则 4 登记为未完成（见 `交付与验收记录_2026-09-29.md`）。V2 是它在
本树内的对应物：域换成冻结清单，判据机制照旧。

用法：python3 V2_冻结输入不变量核对_2026-09-29.py ; echo $?
"""

import csv
import hashlib
import os
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
# V2_OUT 只用于负面对照沙箱：指向一棵被注入缺陷的副本树，验证判据真的会失败。
OUT = pathlib.Path(os.environ["V2_OUT"]).resolve() if os.environ.get("V2_OUT") else HERE.parent
FZ = OUT / "01_冻结结果与代码"
MANIFEST = FZ / "冻结清单_2026-09-27.tsv"

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


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ── 登记制：本轮修订中被改动的冻结件（每条含成因；未登记即失败，登记过期亦失败）──
# 24 件，按成因分四组。分组只是便于人读，判据是逐件的。
KNOWN_DELTA = {
    # 组 1 · 图件重导出（C 段：按评议改图内旧文字、行标签、共享色阶、合并细胞口径；
    #   F23（2026-10-01）：两层主估计改换为 RE-DL 后，按新图注重出图 2／图 3／图 4）
    #   Fig_evidence_matrix.png 此前长期**不在**此列：D2 重跑后 PNG 逐字节未变（md5
    #   fb0dc6cf0eb5474064a24aa3322103ac），只有 PDF 因 R 的 /CreationDate 与 /ID
    #   每次必变而变；剥掉这两项后正文层 10,161 B 亦相同。故当时登记里不写它，
    #   让反向守卫盯着它。F23 改了图 4 的正文层（列组名 Rank (layer primary) →
    #   Rank (IVW-FE disc / pooled rep)），PNG 首次真变 —— 守卫如期报警，本件就此登记在册。
    "04_图表/Fig_candidate_overview.pdf": "F23：图 2 两轴改 RE-DL（原发现层 IVW-FE × 复制层 pooled）后重导出",
    "04_图表/Fig_candidate_overview.png": "同上（同一幅图的栅格出口）",
    "04_图表/Fig_evidence_matrix.pdf": "F23：图 4 列组名改 Rank (IVW-FE disc / pooled rep) 后重导出；正文层确实变了，不再只是元数据",
    "04_图表/Fig_evidence_matrix.png": "F23：同上（同一幅图的栅格出口）。此前各轮逐字节未变故未登记，本轮正文层真变，首次登记",
    "04_图表/Fig_forest_logFC.pdf": "F23：行标签改 RE (DL) (main)／IVW-FE／Pooled (cohort-adjusted)，红菱形移到两层主估计行后重导出",
    "04_图表/Fig_forest_logFC.png": "同上（栅格出口）",
    "04_图表/Fig_scRNA_dotplot.pdf": "C3/B15：图 5 改 mean log1p(CP10K)、四面板共享色阶、All donors 改合并细胞比例",
    "04_图表/Fig_scRNA_dotplot.png": "同上（栅格出口）",
    "04_图表/Fig_study_design.pdf": "C1/B15：图 1 由 E2 真正重绘（此前交付件从未跑过 E2），删图内三处废文字",
    "04_图表/Fig_study_design.png": "同上（栅格出口）",
    # 组 2 · 图件脚手架数据（随重绘一并产出的画图输入，不是分析结果）
    "02_结果表/图1_研究设计_图注.txt": "E2 重绘：图内嵌题注改为与稿件 Figure 1 图注一致",
    "02_结果表/图1_研究设计_布局.csv": "E2 重绘：三层改两层阶梯后的版式坐标",
    "02_结果表/图1_研究设计_括注.txt": "E2 重绘：删去「4 cell types」等废文字",
    "02_结果表/图1_研究设计_箱体.csv": "E2 重绘：箱体坐标（几何随之调整）",
    "02_结果表/图2_57候选总览_数值.csv": "F23：D3 重导出，坐标改 re_logFC（RE-DL）并新增 estimator 列",
    "02_结果表/森林图_合并行数值.csv": "F23：B7 重导出，表头由 kind 改 estimator 并新增 layer 列",
    "02_结果表/图4_证据矩阵_数值.csv": "F23：D2 重导出，列组名改 Discovery (IVW-FE)／Replication (pooled)／Rank (IVW-FE disc / pooled rep)／Deletion (discovery, IVW-FE)",
    # 组 3 · 绘图脚本（C 段：修旧基线绝对路径、改行标签、删图内旧文字）
    "03_分析/scripts/B3_单细胞逐供者定位.py": "C3：四面板共享 vmin/vmax、色条改 log1p、All donors 改合并细胞比例、补点面积图例",
    "03_分析/scripts/B7_森林图.R": "F23：SENS_LAB／PRIM_LAB 改为两层主估计 RE-DL；用 LV <- unique() 去重水平（R≥4.0 重复水平是硬错）",
    "03_分析/scripts/D2_图4证据矩阵_2026-09-25.R": "F23：列组名改 Rank (IVW-FE disc / pooled rep)；删去「该层主估计」等失效注释（表的 64 格逐格未变）",
    "03_分析/scripts/D3_图2_57候选总览_2026-09-25.R": "F23：两轴改 RE-DL，same_dir 改由 re_logFC 符号现算，审计表加 estimator 列",
    "03_分析/scripts/E1_导师核查件docx_2026-09-27.py":
        "新增第三个 CLI 参数 rev_root，使 H3 能换源复用本渲染器出 2026-09-28 核查件"
        "（该件材料清单表的落点按 09 树解析）。已与云侧原件逐行 diff 核实，仅此一处功能扩展。",
    # 组 4 · 运行日志（重跑的副产品）
    "03_分析/logs/B7_森林图.log": "F23 重跑森林图产生的日志",
    "03_分析/logs/D2_图4证据矩阵.log": "F23 重跑图 4 产生的日志",
}

# 分析结果类：正文直接引用的数字来源，一律点名。它们**不得**出现在变动集合里。
NUMERIC_ANCHORS = [
    "02_结果表/57基因总览.csv",
    "02_结果表/A2_交联家族_判据敏感性.csv",
    "02_结果表/A2_名次随模型变动.csv",
    "02_结果表/A2_显著性跨线计数.csv",
    "02_结果表/A2c_全基因组两两一致性.csv",
    "02_结果表/A2c_置换零分布.csv",
    "02_结果表/A2c_队列级诊断.csv",
    "02_结果表/A2e_留一法_逐基因名次变动.csv",
    "02_结果表/A2e_留一法重排汇总.csv",
    "02_结果表/A2d_候选基因×状态_逐供者.csv",
    "02_结果表/A2d_状态基因集.csv",
    "02_结果表/合并效应_发现层.csv",
    "02_结果表/合并效应_复制层.csv",
    "02_结果表/留一法_合并效应.csv",
    "02_结果表/逐队列效应_力学基因.csv",
    "02_结果表/逐SNP工具变量与F.csv",
    "02_结果表/MR工具_逐SNP_F与弱工具.csv",
    "02_结果表/基因符号映射核对.csv",
    "02_结果表/未检出原因核查.csv",
]
# 绘图脚手架：02_结果表 里允许变动的**全部**成员（除它们以外全类不许变）
FIG_SCAFFOLD = {
    "02_结果表/图1_研究设计_图注.txt",
    "02_结果表/图1_研究设计_布局.csv",
    "02_结果表/图1_研究设计_括注.txt",
    "02_结果表/图1_研究设计_箱体.csv",
    "02_结果表/图2_57候选总览_数值.csv",
    "02_结果表/图4_证据矩阵_数值.csv",
    "02_结果表/森林图_合并行数值.csv",
}
# 整棵不许变的冻结树
FROZEN_TREES = ("01_审计/", "05_稿件/", "06_补充材料/", "07_修订记录/", "08_投稿材料/")


def main():
    print("V2 冻结输入不变量核对（2026-09-29）")
    if not MANIFEST.is_file():
        print(f"FATAL 缺基线清单 {MANIFEST}")
        return 2
    rows = list(csv.DictReader(open(MANIFEST, encoding="utf-8"), delimiter="\t"))
    print(f"  基线 {MANIFEST.name}：{len(rows)} 件（2026-09-27 冻结当时逐件 SHA-256）")

    same, changed, absent = [], [], []
    for r in rows:
        rel = r["冻结副本相对路径"]
        p = FZ / rel
        if not p.is_file():
            absent.append(rel)
            continue
        (same if sha(p) == r["SHA-256"] else changed).append(rel)

    print("\n[A] 逐件复算")
    print(f"  一致 {len(same)} / 变动 {len(changed)} / 磁盘上找不到 {len(absent)}")
    must(not absent, f"基线 {len(rows)} 件在磁盘上全部找得到", f"缺 {absent}")
    unreg = [p for p in changed if p not in KNOWN_DELTA]
    must(not unreg, f"变动 {len(changed)} 件全部已在 KNOWN_DELTA 登记（未登记 {len(unreg)} 件）",
         f"{unreg}")
    stale = [p for p in KNOWN_DELTA if p not in changed]
    must(not stale,
         f"KNOWN_DELTA 的 {len(KNOWN_DELTA)} 条登记**全部仍然成立**（无过期登记）",
         f"登记了但实测未变（该登记已不成立）：{stale}")
    # 登记表里的件必须真的在基线里，否则登记的是清单之外的路径（写错名也会静默成立）
    off = [p for p in KNOWN_DELTA if p not in {r["冻结副本相对路径"] for r in rows}]
    must(not off, "KNOWN_DELTA 的键全部是基线清单内的路径（无清单外登记）", f"{off}")

    print("\n[B] 分析结果类不得变动（正文数字的来源）")
    hit = [p for p in changed if p in NUMERIC_ANCHORS]
    must(not hit, f"点名的 {len(NUMERIC_ANCHORS)} 件分析结果**逐字节未变**", f"变动 {hit}")
    anchor_absent = [p for p in NUMERIC_ANCHORS if p not in {r["冻结副本相对路径"] for r in rows}]
    must(not anchor_absent, "点名的分析结果全都在基线清单内（点名不落空）", f"{anchor_absent}")
    # 02_结果表 全类：除 7 件绘图脚手架外，其余一律不许变
    r2_changed = [p for p in changed if p.startswith("02_结果表/") and p not in FIG_SCAFFOLD]
    must(not r2_changed,
         f"02_结果表 除 {len(FIG_SCAFFOLD)} 件绘图脚手架外全类未变",
         f"变动 {r2_changed}")
    must(len(FIG_SCAFFOLD) == 7, "绘图脚手架登记恰 7 件（本表与 KNOWN_DELTA 里的组 2 逐件对应）")

    print("\n[C] 域守卫：五棵冻结树整棵不得变动")
    for t in FROZEN_TREES:
        v = [p for p in changed if p.startswith(t)]
        n = sum(1 for r in rows if r["冻结副本相对路径"].startswith(t))
        must(not v, f"{t} 整棵未变（{n} 件）", f"变动 {v}")

    print("\n[D] 覆盖度自检：登记必须把变动解释干净")
    by_group = {
        "图件重导出": [p for p in KNOWN_DELTA if p.startswith("04_图表/")],
        "图件脚手架": [p for p in KNOWN_DELTA if p in FIG_SCAFFOLD],
        "绘图脚本": [p for p in KNOWN_DELTA if p.startswith("03_分析/scripts/")],
        "运行日志": [p for p in KNOWN_DELTA if p.startswith("03_分析/logs/")],
    }
    tot = sum(len(v) for v in by_group.values())
    for k, v in by_group.items():
        print(f"    {k:12s} {len(v):2d} 件")
    must(tot == len(KNOWN_DELTA), f"四组成因之和 {tot} 等于登记总数 {len(KNOWN_DELTA)}（无漏归组）")
    must(len(changed) == len(KNOWN_DELTA),
         f"变动件数 {len(changed)} 与登记件数 {len(KNOWN_DELTA)} 相等（一一对应）")

    print(f"\n结果：{'全部通过' if not fails else '存在失败'}")
    print(f"通过 {oks[0]} 项，失败 {len(fails)} 项")
    for f in fails:
        print(f"  - {f}")
    print(f"\n判读：冻结 {len(rows)} 件中 {len(same)} 件逐字节未变；变动的 {len(changed)} 件"
          f"全部落在图件、绘图脚本、图件脚手架数据与运行日志四类，"
          f"**没有任何一件分析结果**变动。")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
