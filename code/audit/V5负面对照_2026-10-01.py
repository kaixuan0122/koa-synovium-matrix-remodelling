#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V5 负面对照（2026-10-01）

目的：证明 V5 的判据**真的会失败**。V5 的 `KNOWN_DELTA` 是空的，于是 [A]—[C] 在一棵干净树上
必然「全部通过」—— 恒真的判据与没跑没有区别。本件逐个把缺陷注入沙箱副本，把「哪些判据被
触发」量出来。全部在 `/private/tmp` 的沙箱副本里做，交付树全程只读（规则 1）。

五个对照，各自只测**一类**判据：

  NC-1  往一件**点名的分析结果**（合并效应_发现层.csv）尾部追加一个字节
        → 预期只触发 [A] 变动为 0 与 [B] 两条。注意它**不会**触发 [D]／[E]：该件在 09-27
          基线下本来就已变动，再追加一字节不改变「新旧两代之间的变动集合」。这正是要在
          对照里看清的事 —— 否则会误以为 [A]—[C] 与 [D]／[E] 互相兜底。
  NC-2  改一件**在两代基线之间逐位未变**的件（从清单里现取一棵冻结树里的文件）
        → 预期触发 [A]＋[C]＋[D]＋[E]。这一条才测到 [D]／[E] —— 它们比的是**两代之间**。
  NC-3  **基线对照**：沙箱原样 → 预期 rc=0、0 FAIL（证明前两条的失败不是沙箱本身带来的）
  NC-4  把沙箱里 F34 的 `EXPECTED_CHANGED` 删掉一件（磁盘一字未动）
        → 预期触发 [E] 的两条（读到的条目数、与实测同解）。这一条测的是「两处独立登记必须
          同解」本身：没有它，[E] 会退化成把 R11 重打一遍。
  NC-5  把新清单里 `02_结果表/图1_研究设计_箱体.csv` 的「字节」格 336 改成 337（磁盘不动）
        → 预期触发 [D] 的两条字节判据。NC-2 测不到它们 —— 那两条比的是**两代清单记的字节列**，
          改磁盘不会动到清单，故须单独测一次。

注：NC-1 与 NC-2 都**不**触发 [D] 的字节判据，这本身就是被量出来的事实（见 NC-1／NC-2 的
「说明」断言）—— 判据之间不是互相兜底的。

用法：python3 V5负面对照_2026-10-01.py ; echo $?
"""

import csv
import os
import pathlib
import shutil
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parent
V5 = HERE / "V5_冻结输入不变量核对_2026-10-01.py"
SANDBOX = pathlib.Path("/private/tmp") / "v5_negctl_2026-10-01"
FZ_REL = "01_冻结结果与代码/"
MANIFEST_REL = FZ_REL + "冻结清单_2026-10-01.tsv"
F34_REL = FZ_REL + "F34_冻结输入_2026-10-01.py"
ANCHOR_MAN = "02_结果表/合并效应_发现层.csv"      # 清单里的键（相对冻结根）
ANCHOR = FZ_REL + ANCHOR_MAN                     # 磁盘路径（相对树根）

fails = []
oks = [0]


def must(ok, msg, detail=None):
    if ok:
        oks[0] += 1
        print(f"  OK   {msg}")
    else:
        fails.append(msg)
        print(f"  FAIL {msg}" + (f"  | {detail}" if detail else ""))
    return ok


def fail_msgs(out):
    """取 V5 的 FAIL **消息**（不含 '  | detail' 段）——预期集合按消息比，detail 里带数字。"""
    r = []
    for ln in out.splitlines():
        s = ln.strip()
        if s.startswith("FAIL "):
            r.append(s[5:].split("  | ")[0].strip())
    return r


def run_v5(root):
    env = {"PATH": "/usr/bin:/bin:/usr/local/bin", "V5_OUT": str(root)}
    p = subprocess.run([sys.executable, str(root / "00_版本审计" / V5.name)],
                       capture_output=True, text=True, env=env)
    return p.returncode, p.stdout + p.stderr


def fresh():
    if SANDBOX.exists():
        shutil.rmtree(SANDBOX)
    SANDBOX.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(OUT, SANDBOX)


def expect(got, want, tag):
    """want = 消息子串列表；按**子串**匹配（消息里带计数，写死整句会在件数一变时假失败）。"""
    for w in want:
        must(any(w in g for g in got), f"{tag} 预期 FAIL：{w[:46]}…", f"实得 FAIL 集合 {got}")
    must(len(got) == len(want), f"{tag} FAIL 集合恰 {len(want)} 条（无连带误伤）", f"实得 {got}")


print("V5 负面对照：判据必须能失败")
if not V5.is_file():
    print(f"FATAL 缺 {V5}")
    sys.exit(2)

# 从清单里现取一件「在两代基线之间逐位未变」的冻结树文件（NC-2 用）。
# 写死路径会在清单换代时指向别处；这里按前缀现取并自证它确实未变。
_rows = list(csv.DictReader(open(OUT / FZ_REL / "冻结清单_2026-10-01.tsv",
                                 encoding="utf-8"), delimiter="\t"))
_prior = {r["冻结副本相对路径"]: r["SHA-256"]
          for r in csv.DictReader(open(OUT / FZ_REL / "冻结清单_2026-09-27.tsv",
                                       encoding="utf-8"), delimiter="\t")}
FROZEN_MAN = next(r["冻结副本相对路径"] for r in _rows
                  if r["冻结副本相对路径"].startswith("05_稿件/"))
FROZEN_IN = FZ_REL + FROZEN_MAN
must(FROZEN_MAN in _prior and _prior[FROZEN_MAN] ==
     next(r["SHA-256"] for r in _rows if r["冻结副本相对路径"] == FROZEN_MAN),
     f"NC-2 前提：现取的 {FROZEN_MAN} 在两代基线之间逐位未变（改它才会动到 [D]／[E]）")

# ── NC-3 基线：沙箱原样 ──────────────────────────────────────────────────
print("\n[NC-3] 基线（沙箱原样，不改任何东西）")
fresh()
rc, out = run_v5(SANDBOX)
must(rc == 0, "基线 rc=0", f"实得 {rc}")
must("全部通过" in out and "失败 0 项" in out, "基线 0 FAIL",
     out.strip().splitlines()[-1] if out.strip() else "")

# ── NC-1 篡改一件点名的分析结果 ──────────────────────────────────────────
print("\n[NC-1] 往 合并效应_发现层.csv 追加一个字节（测 [A]／[B]）")
fresh()
tgt = SANDBOX / ANCHOR
before = tgt.stat().st_size
with open(tgt, "ab") as fh:
    fh.write(b"\n")
must(tgt.stat().st_size == before + 1,
     "NC-1 前提：改写确已生效（字节数 +1）", f"{before} → {tgt.stat().st_size}")
rc, out = run_v5(SANDBOX)
must(rc == 1, "NC-1 rc=1", f"实得 {rc}")
expect(fail_msgs(out), [
    "新基线下变动",
    "点名的 19 件分析结果**逐字节未变**",
    "02_结果表 除 7 件绘图脚手架外全类未变",
], "NC-1")
must("过户：34 件 SHA-256 变动" in out,
     "NC-1 说明：该件本就在两代基线之间变动，故 [D]／[E] 仍全过 —— 两套判据不是互相兜底的",
     out[-400:])

# ── NC-2 改一件两代之间未变的冻结树文件 → 测 [C]／[D]／[E] ────────────────
print(f"\n[NC-2] 改一件两代之间未变的冻结树文件（{FROZEN_MAN}）")
fresh()
tgt2 = SANDBOX / FROZEN_IN
before2 = tgt2.stat().st_size
with open(tgt2, "ab") as fh:
    fh.write(b"\n")
must(tgt2.stat().st_size == before2 + 1,
     "NC-2 前提：改写确已生效（字节数 +1）", f"{before2} → {tgt2.stat().st_size}")
rc, out = run_v5(SANDBOX)
must(rc == 1, "NC-2 rc=1", f"实得 {rc}")
expect(fail_msgs(out), [
    "新基线下变动",
    "05_稿件/ 整棵未变",
    "过户：34 件 SHA-256 变动",
    "实测变动扣掉 V2 那 24 件",
    "F34 的 EXPECTED_CHANGED 与实测变动集合同解",
], "NC-2")   # 恰 5 条 = [A]1 ＋ [C]1 ＋ [D]1 ＋ [E]2；集合里没有 [B] 的消息 —— 判据各管一段
# 注意 [D] 的**字节数**那条在此不响：它比的是**两代清单记的字节列**，不是磁盘 ——
# 改磁盘不会动到清单。它要由 NC-5 单独测。

# ── NC-4 删掉沙箱里 F34 的一条登记（磁盘一字未动）→ 测 [E] 的同解判据 ──────
print("\n[NC-4] 把沙箱里 F34 的 EXPECTED_CHANGED 删掉一件（磁盘不动，测 [E] 的同解判据）")
fresh()
f34 = SANDBOX / F34_REL
src = f34.read_text(encoding="utf-8")
LINE = '    "04_图表/Fig_study_design.png",\n'
must(src.count(LINE) == 1, "NC-4 前提：待删的那一行在 F34 源里恰一处", f"实得 {src.count(LINE)}")
f34.write_text(src.replace(LINE, "", 1), encoding="utf-8")
must('"04_图表/Fig_study_design.png",' not in f34.read_text(encoding="utf-8"),
     "NC-4 前提：改写确已生效（该行已从沙箱副本里删去）")
rc, out = run_v5(SANDBOX)
must(rc == 1, "NC-4 rc=1", f"实得 {rc}")
# 两条都预期会响：[E] 里「读到的条目数」与「与实测同解」是**两条**判据，删一件两条同时触发。
# 这正是本对照要量出来的事 —— 计数守卫不替同解守卫兜底，反之亦然。
expect(fail_msgs(out), [
    "从 F34 源文件里读到 EXPECTED_CHANGED 的 34 件",
    "F34 的 EXPECTED_CHANGED 与实测变动集合同解",
], "NC-4")

# ── NC-5 改新清单里箱体那一件的「字节」格 → 测 [D] 的两条字节判据 ──────────
print("\n[NC-5] 把新清单里 图1_研究设计_箱体.csv 的字节格 336 改成 337（磁盘不动）")
fresh()
man = SANDBOX / MANIFEST_REL
lines = man.read_text(encoding="utf-8").splitlines(keepends=True)
# 须先 rstrip("\n")：keepends=True 下最后一个字段带着换行，直接比路径会一条也命不中。
hit_i = [i for i, l in enumerate(lines)
         if l.rstrip("\n").split("\t")[4:5] == ["02_结果表/图1_研究设计_箱体.csv"]]
must(len(hit_i) == 1, "NC-5 前提：该件在新清单里恰一行", f"实得 {len(hit_i)}")
i = hit_i[0]
cells = lines[i].rstrip("\n").split("\t")
must(cells[1] == "336", "NC-5 前提：该行的字节格现为 336（与前代相同）", f"实得 {cells[1]}")
cells[1] = "337"
lines[i] = "\t".join(cells) + "\n"
man.write_text("".join(lines), encoding="utf-8")
must(man.read_text(encoding="utf-8").splitlines()[i].split("\t")[1] == "337",
     "NC-5 前提：改写确已生效")
rc, out = run_v5(SANDBOX)
must(rc == 1, "NC-5 rc=1", f"实得 {rc}")
expect(fail_msgs(out), [
    "变动件里**字节数**不同的只有 33 件",
    "点名的那一件",
], "NC-5")

shutil.rmtree(SANDBOX, ignore_errors=True)
print(f"\n结果：{'全部通过' if not fails else '存在失败'}")
print(f"通过 {oks[0]} 项，失败 {len(fails)} 项")
for f in fails:
    print(f"  - {f}")
print("\n判读：NC-3 证明沙箱本身不产生 FAIL；NC-1 证明 [A]／[B] 会失败、且它们与 [D]／[E]"
      "\n     不是互相兜底的（该件本就在两代之间变动，[D]／[E] 照样全过）；NC-2 证明 [C]／[D]／[E]"
      "\n     会失败而 [B] 不会；NC-4 证明「两处独立登记必须同解」会失败；"
      "\n     NC-5 证明 [D] 的字节判据会失败（NC-2 测不到它 —— 它比的是清单的字节列，不是磁盘）。")
sys.exit(0 if not fails else 1)
