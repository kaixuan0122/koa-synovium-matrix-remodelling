#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V5 — 冻结输入不变量核对（2026-10-01）

**为什么换基线。** V2（2026-09-29）的基线是 `冻结清单_2026-09-27.tsv`：它把 09-27 当时的
165 件逐件字节当作不变量。第 8—11 轮按设计改了其中 **34 件**，其中 **6 件是 `NUMERIC_ANCHORS`
点名的分析结果**、**1 件在 `01_审计/` 这棵「整棵不许变」的树里**（主估计由 `IVW-FE`／
`pooled(HC3)` 升为 `REML+HK`，B2／A2／B8 三处 R 脚本原地重跑重出）。旧基线因此**在原理上
不可能再成立**，与判据松紧无关。

在 i 轮成品树上实测跑了一遍 V2（只读），失败 5 条，未登记恰 10 件 —— 就是 i 轮那 10 件：

    FAIL 变动 34 件全部已在 KNOWN_DELTA 登记（未登记 10 件）
    FAIL 点名的 19 件分析结果**逐字节未变**（变动 6 件：合并效应_{发现,复制}层、
         57基因总览、A2_{交联家族_判据敏感性, 名次随模型变动, 显著性跨线计数}）
    FAIL 02_结果表 除 7 件绘图脚手架外全类未变（同上 6 件）
    FAIL 01_审计/ 整棵未变（13 件）｜ 变动 01_审计/A2_判据敏感性_结果.md
    FAIL 变动件数 34 与登记件数 24 相等（一一对应）

（V2 末尾那句「没有任何一件分析结果变动」是它的**固定结论句**，此刻与它自己的 FAIL 相矛盾
—— 本件不改 V2（它是 09-29 那一刻的记录），只把这一处在报告里点名。）

**V5 不是把 V2 改松。** V2 一字未动、原地留存。V5 是新文件、新基线
（`冻结清单_2026-10-01.tsv`，由 F34 只哈希不拷贝重建），`KNOWN_DELTA` **清空**：新基线
就是一个固定点。V2 的判据逐条照搬，只去掉两条**在空登记下恒真**的（「未登记变动」与
「变动为 0」在空登记下是同一句话；「过期登记」恒为空集；「四组成因之和 == 登记总数」退化成
`0 == 0`）—— 这是随基线换代自然消失的，不是放宽。

**但空登记的核对脚本会变成摆设**：在一棵干净树上「全部通过」是必然的，那与没跑没有区别。
故 V5 另加两条承重判据：

  · **[D] 过户判据**：把新旧两份基线**逐行**对上。件数、路径与**行序**必须相同；SHA-256
    不同者**恰为**登记集合；其余 131 件逐位相同；字节列不同者 33 件 —— 因为
    `02_结果表/图1_研究设计_箱体.csv` 内容变了而字节数恰好仍是 336，**单比字节数会漏一件**，
    故本件所有判据一律比 SHA-256。
  · **[E] 交叉核对**：变动的 34 件在两处**独立写下**过 —— V2 的 `KNOWN_DELTA`（24 条，09-29
    那一轮的登记）与 F34 的 `EXPECTED_CHANGED`（34 件，基线重建时的登记）。两处都从**源文件
    里读出来**（`ast` 取常量，不重打一遍），加上实测的变动集合，三者必须同解。V2 那一份
    就是上面那段实测失败输出的成因记录。

用法：python3 V5_冻结输入不变量核对_2026-10-01.py ; echo $?
      V5_OUT=<沙箱树根> python3 V5_…（负面对照用，见 V5负面对照_2026-10-01.py）
"""

import ast
import csv
import hashlib
import os
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
# V5_OUT 只用于负面对照沙箱：指向一棵被注入缺陷的副本树，验证判据真的会失败。
OUT = pathlib.Path(os.environ["V5_OUT"]).resolve() if os.environ.get("V5_OUT") else HERE.parent
FZ = OUT / "01_冻结结果与代码"
MANIFEST = FZ / "冻结清单_2026-10-01.tsv"
PRIOR = FZ / "冻结清单_2026-09-27.tsv"
V2_SRC = OUT / "00_版本审计" / "V2_冻结输入不变量核对_2026-09-29.py"
F34_SRC = FZ / "F34_冻结输入_2026-10-01.py"

# 新基线即固定点：本轮不登记任何变动。留这个空字典是为了让 [A]—[C] 的判据形态与 V2 逐条对应。
KNOWN_DELTA = {}

# i 轮（2026-10-01，主估计由 IVW-FE / pooled(HC3) 升为 REML+HK）按设计改动的 10 件。
# 因果不在此重述 24 件历史登记（那 24 条的成因全文在 V2 的 KNOWN_DELTA 里，本件 [E] 直接读它，
# 避免同一份登记存两份而漂移）。
R11 = {
    "01_审计/A2_判据敏感性_结果.md":
        "A2 重跑：判据敏感性表的估计量列随主分析同步为 REML+HK（该件在「整棵不许变」的树里，"
        "是本轮唯一必须换基线的直接原因）",
    "02_结果表/合并效应_发现层.csv":
        "B2 重跑：新增 hk_logFC/hk_SE/hk_P/hk_CI_*/hk_tau2 与 pi_low/pi_high，"
        "primary_logFC 等改指 REML+HK（原指 IVW-FE）",
    "02_结果表/合并效应_复制层.csv": "B2 重跑（同上；该层 primary 原指 pooled(HC3)）",
    "02_结果表/A2_交联家族_判据敏感性.csv":
        "A2 重跑产物：判据敏感性表的估计量列随主分析同步",
    "02_结果表/A2_名次随模型变动.csv": "A2 重跑产物：名次随模型变动表的估计量列同步",
    "02_结果表/A2_显著性跨线计数.csv":
        "A2 重跑产物：S4d 的三行由 FE-vs-RE 换成 FE-vs-primary 与 pooled-vs-primary",
    "02_结果表/57基因总览.csv":
        "B8 重跑：吃 B2 新表，disc_/rep_ 两列组随估计量改名（disc_FE_*→disc_ivw_*、"
        "disc_RE_*→disc_re_*），disc_logFC[…] 现为 REML+HK",
    "03_分析/scripts/A2_交联家族_判据敏感性.R":
        "改 EXPECT_METHOD 等常量与输出的估计量列名（主估计量升为 REML+HK）",
    "03_分析/scripts/B2_合并效应与留一法.R":
        "新增 metafor::rma(method=\"REML\", test=\"knha\") 与 predict() 预测区间；"
        "primary_* 改指它，原 primary（IVW-FE／pooled）降为描述性敏感性",
    "03_分析/scripts/B8_57基因总览.R": "列名常量随 B2 的新表头同步（含重命名后的 re_* 列）",
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

# 过户判据的冻结量（2026-10-01 实测，逐条见 [D]）
N_ROWS = 165          # 两份清单各 165 件
N_TRANSITION = 34     # 09-27 → 10-01 的 SHA-256 变动件数
N_IDENTICAL = 131     # 逐位未变件数
N_BYTE_SAME = 33      # 变动件里「字节数恰好不变」的件数（故判据不得只比字节数）
BYTE_INVARIANT = "02_结果表/图1_研究设计_箱体.csv"   # 上面那唯一一件

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


def read_manifest(p):
    """返回 {冻结副本相对路径: 行}，并保序返回路径列表。"""
    rows = list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))
    return rows, [r["冻结副本相对路径"] for r in rows]


def const_from(src, name):
    """从**源文件**里取一个模块级常量（字面量）。

    用 ast 而不是 import：F34 的模块级正文没有 `if __name__` 守卫，import 会**真的跑一遍
    清单重建**（会写盘）；V2 虽守卫了，但两处取法一致更不容易出错。取不到返回 None。
    """
    if not src.is_file():
        return None
    for node in ast.parse(src.read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == name for t in node.targets):
            try:
                return ast.literal_eval(node.value)
            except ValueError:
                return None
    return None


def main():
    print("V5 冻结输入不变量核对（2026-10-01）")
    if not MANIFEST.is_file():
        print(f"FATAL 缺基线清单 {MANIFEST}")
        return 2
    if not PRIOR.is_file():
        print(f"FATAL 缺上一代基线 {PRIOR} —— [D] 过户判据无从做起，中止（不猜）")
        return 2

    rows, order = read_manifest(MANIFEST)
    prior_rows, prior_order = read_manifest(PRIOR)
    print(f"  基线 {MANIFEST.name}：{len(rows)} 件（由 F34 只哈希不拷贝重建，2026-10-01）")
    print(f"  上代 {PRIOR.name}：{len(prior_rows)} 件（09-27 冻结当时，原地留存供追溯）")

    SHA = {r["冻结副本相对路径"]: r["SHA-256"] for r in rows}
    BYTES = {r["冻结副本相对路径"]: r["字节"] for r in rows}
    measured = {k: (sha(FZ / k) if (FZ / k).is_file() else None) for k in order}

    print("\n[A] 逐件复算（新基线是固定点）")
    absent = sorted(k for k, v in measured.items() if v is None)
    changed = sorted(k for k in order if measured[k] is not None and measured[k] != SHA[k])
    same = [k for k in order if k not in set(changed)]
    print(f"  一致 {len(same)} / 变动 {len(changed)} / 磁盘上找不到 {len(absent)}")
    must(not absent, f"基线 {len(rows)} 件在磁盘上全部找得到", f"缺 {absent}")
    # 空登记下「未登记变动」与「变动为 0」是同一句话，故只留后者（V2 的两条在此合并）
    must(not changed,
         f"新基线下变动 {len(changed)} 件（新基线是固定点：清单里的件逐字节未变）",
         f"变动 {[(p, KNOWN_DELTA.get(p)) for p in changed]}")
    off = [p for p in KNOWN_DELTA if p not in set(order)]
    must(not off, "KNOWN_DELTA 的键全部是基线清单内的路径（无清单外登记）", f"{off}")

    print("\n[B] 分析结果类不得变动（正文数字的来源）")
    hit = [p for p in changed if p in NUMERIC_ANCHORS]
    must(not hit, f"点名的 {len(NUMERIC_ANCHORS)} 件分析结果**逐字节未变**", f"变动 {hit}")
    anchor_absent = [p for p in NUMERIC_ANCHORS if p not in set(order)]
    must(not anchor_absent, "点名的分析结果全都在基线清单内（点名不落空）", f"{anchor_absent}")
    r2_changed = [p for p in changed if p.startswith("02_结果表/") and p not in FIG_SCAFFOLD]
    must(not r2_changed, f"02_结果表 除 {len(FIG_SCAFFOLD)} 件绘图脚手架外全类未变",
         f"变动 {r2_changed}")
    must(len(FIG_SCAFFOLD) == 7, "绘图脚手架登记恰 7 件（与 V2 的该表逐件对应）")

    print("\n[C] 域守卫：五棵冻结树整棵不得变动")
    for t in FROZEN_TREES:
        v = [p for p in changed if p.startswith(t)]
        n = sum(1 for r in rows if r["冻结副本相对路径"].startswith(t))
        must(not v, f"{t} 整棵未变（{n} 件）", f"变动 {v}")

    # ────────────────────────────────────────────────────────────────────────
    print("\n[D] ★ 过户判据：新旧两代基线必须逐行对上（本件的承重项）")
    must(len(prior_rows) == len(rows) == N_ROWS,
         f"两代清单都是 {N_ROWS} 件", f"{len(prior_rows)} / {len(rows)}")
    must(list(prior_rows[0].keys()) == list(rows[0].keys()),
         "两代清单的五列同序同名（列换了位置会让下面的逐行比对失去意义）",
         f"{list(prior_rows[0].keys())}")
    must(prior_order == order,
         "两代清单的路径**逐行同序同值**（顺序变了，逐行比对就不可信）",
         f"首个不同的行号 {next((i for i, (a, b) in enumerate(zip(prior_order, order)) if a != b), None)}")
    prior = {r["冻结副本相对路径"]: r for r in prior_rows}
    tr_meas = sorted(k for k in order if prior[k]["SHA-256"] != measured[k])
    same_meas = [k for k in order if prior[k]["SHA-256"] == measured[k]]
    byte_meas = sorted(k for k in tr_meas if prior[k]["字节"] != BYTES[k])
    must(len(tr_meas) == N_TRANSITION and len(same_meas) == N_IDENTICAL,
         f"过户：{N_TRANSITION} 件 SHA-256 变动 / {N_IDENTICAL} 件逐位未变",
         f"实得 {len(tr_meas)} / {len(same_meas)}")
    must(len(byte_meas) == N_BYTE_SAME,
         f"变动件里**字节数**不同的只有 {N_BYTE_SAME} 件（故判据一律比 SHA-256，不比字节数）",
         f"实得 {len(byte_meas)}")
    must(BYTE_INVARIANT in tr_meas and BYTE_INVARIANT not in byte_meas,
         f"点名的那一件：{BYTE_INVARIANT} 内容变了而字节数未变 —— 只比字节数会漏它，"
         f"本件据此把「字节不变」排除在判据之外")

    print("\n[E] ★ 交叉核对：两处**独立写下**的登记 ＋ 实测，三者必须同解")
    v2k = const_from(V2_SRC, "KNOWN_DELTA")
    f34c = const_from(F34_SRC, "EXPECTED_CHANGED")
    must(isinstance(v2k, dict) and len(v2k) == 24,
         "从 V2 源文件里读到 KNOWN_DELTA 的 24 条历史登记（09-29 那一轮）",
         f"实得 {type(v2k).__name__} / {len(v2k) if isinstance(v2k, dict) else '—'}")
    must(isinstance(f34c, list) and len(f34c) == 34,
         "从 F34 源文件里读到 EXPECTED_CHANGED 的 34 件（基线重建时登记）",
         f"实得 {type(f34c).__name__} / {len(f34c) if isinstance(f34c, list) else '—'}")
    if not (isinstance(v2k, dict) and isinstance(f34c, list)):
        print("\nFATAL 两处登记读不到，[E] 无从判起（不静默跳过）")
        return 2
    v2ks, f34s, trs = set(v2k), set(f34c), set(tr_meas)
    must(not (v2ks & set(R11)),
         "V2 的 24 条与 i 轮的 10 件**无一重叠**（两轮改动不同批）", f"{sorted(v2ks & set(R11))}")
    must(v2ks <= trs,
         "V2 登记的 24 件**全部**落在实测变动集合里（历史登记与实测相符）",
         f"登记了但实测未变：{sorted(v2ks - trs)}")
    must(trs - v2ks == set(R11),
         f"实测变动扣掉 V2 那 24 件，**恰为** i 轮的 {len(R11)} 件（这就是 V2 报「未登记 10 件」的那 10 件）",
         f"差集 {sorted(trs - v2ks)}")
    must(f34s == trs,
         "F34 的 EXPECTED_CHANGED 与实测变动集合同解（基线重建的登记非事后追认）",
         f"多 {sorted(f34s - trs)} / 少 {sorted(trs - f34s)}")
    must(len(set(R11) & set(NUMERIC_ANCHORS)) == 6,
         "i 轮改动里落在 NUMERIC_ANCHORS 上的恰 6 件（V2 的 [B] 报的就是这 6 件）",
         f"{sorted(set(R11) & set(NUMERIC_ANCHORS))}")
    must("01_审计/A2_判据敏感性_结果.md" in R11,
         "i 轮改动里含 01_审计/ 那一件（V2 的 [C] 域守卫报的就是它）")

    print(f"\n结果：{'全部通过' if not fails else '存在失败'}")
    print(f"通过 {oks[0]} 项，失败 {len(fails)} 项")
    for f in fails:
        print(f"  - {f}")
    print(f"\n判读：新基线 {len(rows)} 件逐字节全部未变（变动 0 件）—— 它是本轮交付的固定点。"
          f"\n     过户：上代 165 件 → 本代 165 件的变动恰 {len(tr_meas)} 件，"
          f"与 V2 的 24 条历史登记 ∪ i 轮的 10 件逐件对应，{len(same_meas)} 件逐位未变。")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
