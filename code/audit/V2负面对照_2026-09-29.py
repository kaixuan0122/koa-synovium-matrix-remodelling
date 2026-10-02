#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V2 负面对照（2026-09-29）

目的：证明 V2 的四条判据**真的会失败**，而不是恒真的摆设——判据恒真时报告一样是
      「全部通过」。每个对照都在 /private/tmp 的沙箱副本里做，交付树全程只读（规则 1）。

三个对照：
  NC-1  往一件**点名的分析结果**（57基因总览.csv）尾部追加一个字节
        → 预期 FAIL：未登记 1 件 ／ 点名分析结果变动 ／ 02_结果表 全类变动
  NC-2  往沙箱的 V2 副本里给一件**实测未变**的件补一条登记
        （01_审计/A2_判据敏感性_结果.md）
        → 预期 FAIL：反向守卫（该登记已不成立）
        这一条测的是反向守卫本身——没有它，登记表会漂移成永远成立的注脚。
  NC-3  **基线对照**：沙箱原样 → 预期 rc=0、0 FAIL（证明上面两条的失败不是沙箱本身带来的）

用法：python3 V2负面对照_2026-09-29.py ; echo $?
"""

import hashlib
import pathlib
import shutil
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parent
V2 = HERE / "V2_冻结输入不变量核对_2026-09-29.py"
SANDBOX = pathlib.Path("/private/tmp") / "v2_negctl_2026-09-29"
FZ_REL = "01_冻结结果与代码/"
ANCHOR_MAN = "02_结果表/57基因总览.csv"           # 清单里的键（相对冻结根）
ANCHOR = FZ_REL + ANCHOR_MAN                      # 磁盘路径（相对树根）
UNCHANGED_MAN = "01_审计/A2_判据敏感性_结果.md"    # 实测未变，用它的 sha 反查基线值
UNCHANGED_IN = FZ_REL + UNCHANGED_MAN             # 磁盘路径


def fail_msgs(out):
    """取出 V2 的 FAIL **消息**（不含 '  | detail' 段）——预期集合按消息比，detail 里带数字。"""
    r = []
    for ln in out.splitlines():
        s = ln.strip()
        if s.startswith("FAIL "):
            r.append(s[5:].split("  | ")[0].strip())
    return r

fails = []
oks = [0]


def must(ok, msg, detail=None):
    if ok:
        oks[0] += 1
        print(f"  OK   {msg}")
    else:
        fails.append(msg)
        print(f"  FAIL {msg}" + (f"  | {detail}" if detail else ""))


def run_v2(root):
    """在 root 上跑 V2，返回 (rc, stdout)。"""
    env = {"PATH": "/usr/bin:/bin:/usr/local/bin", "V2_OUT": str(root)}
    p = subprocess.run([sys.executable, str(root / "00_版本审计" / V2.name)],
                       capture_output=True, text=True, env=env)
    return p.returncode, p.stdout + p.stderr


def fresh():
    if SANDBOX.exists():
        shutil.rmtree(SANDBOX)
    SANDBOX.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(OUT, SANDBOX)


print("V2 负面对照：判据必须能失败")
if not V2.is_file():
    print(f"FATAL 缺 {V2}")
    sys.exit(2)

# ── NC-3 基线：沙箱原样 ──────────────────────────────────────────────────
print("\n[NC-3] 基线（沙箱原样，不改任何东西）")
fresh()
rc, out = run_v2(SANDBOX)
must(rc == 0, "基线 rc=0", f"实得 {rc}")
must("全部通过" in out and "失败 0 项" in out, "基线 0 FAIL",
     out.strip().splitlines()[-1] if out.strip() else "")

# ── NC-1 篡改一件点名的分析结果 ──────────────────────────────────────────
print("\n[NC-1] 往 57基因总览.csv 追加一个字节")
fresh()
tgt = SANDBOX / ANCHOR
before = tgt.stat().st_size
with open(tgt, "ab") as fh:
    fh.write(b"\n")
must(tgt.stat().st_size == before + 1,
     "负面对照前提：改写确已生效（字节数 +1）", f"{before} → {tgt.stat().st_size}")
rc, out = run_v2(SANDBOX)
must(rc == 1, "NC-1 rc=1", f"实得 {rc}")
WANT1 = [
    # 该消息里带计数（"变动 N 件…未登记 M 件"），故按不含计数的子串匹配：
    # 写死旧计数会在件数一变时把「判据正常工作」误报成对照失败。
    "全部已在 KNOWN_DELTA 登记",
    "点名的 19 件分析结果**逐字节未变**",
    "02_结果表 除 6 件绘图脚手架外全类未变",
    "变动件数 23 与登记件数 22 相等（一一对应）",
]
got = fail_msgs(out)
for w in WANT1:
    # 注意是**子串**匹配：`w in got` 在列表上是元素相等，消息里带计数就永远比不中。
    must(any(w in g for g in got), f"NC-1 预期 FAIL：{w[:44]}…", f"实得 FAIL 集合 {got}")
must(len(got) == len(WANT1), f"NC-1 FAIL 集合恰 {len(WANT1)} 条（无连带误伤）", f"实得 {got}")

# ── NC-2 给一件实测未变的件补一条登记 → 反向守卫必须失败 ──────────────────
print("\n[NC-2] 给一件实测未变的件补登记（测反向守卫）")
fresh()
v2s = SANDBOX / "00_版本审计" / V2.name
src = v2s.read_text(encoding="utf-8")
# 该件在基线里的 SHA-256，用于自证「它确实没变」（对照不是空转）
import csv
rows = {r["冻结副本相对路径"]: r["SHA-256"]
        for r in csv.DictReader(open(SANDBOX / "01_冻结结果与代码" / "冻结清单_2026-09-27.tsv",
                                     encoding="utf-8"), delimiter="\t")}
h = hashlib.sha256((SANDBOX / UNCHANGED_IN).read_bytes()).hexdigest()
must(h == rows[UNCHANGED_MAN], "NC-2 前提：该件的实测哈希等于基线哈希（确未变动）")
OLD = '    # 组 4 · 运行日志（重跑的副产品）'
NEW = ('    "%s": "NC-2 注入的假登记（该件实测未变）",\n' % UNCHANGED_MAN) + OLD
must(src.count(OLD) == 1, "NC-2 前提：注入锚点在 V2 源里恰一处", f"实得 {src.count(OLD)}")
v2s.write_text(src.replace(OLD, NEW, 1), encoding="utf-8")
# 回读自证改写生效（空转是对照的典型失败形态）
must(UNCHANGED_MAN in v2s.read_text(encoding="utf-8"),
     "NC-2 前提：改写确已生效（假登记已写入沙箱的 V2 副本）")
rc, out = run_v2(SANDBOX)
must(rc == 1, "NC-2 rc=1", f"实得 {rc}")
got2 = fail_msgs(out)
must(any("全部仍然成立" in g for g in got2),
     "NC-2 预期 FAIL：反向守卫报「登记仍全部成立」（即假登记未被识破？）", f"实得 {got2}")
must("该登记已不成立" in out,
     "NC-2 预期：反向守卫的 detail 点名了那件已不成立的登记", f"实得 {got2}")
must(len(got2) == 3, f"NC-2 FAIL 集合恰 3 条（反向守卫 ＋ 两条计数一致性）", f"实得 {got2}")

shutil.rmtree(SANDBOX, ignore_errors=True)
print(f"\n结果：{'全部通过' if not fails else '存在失败'}")
print(f"通过 {oks[0]} 项，失败 {len(fails)} 项")
for f in fails:
    print(f"  - {f}")
print("\n判读：NC-1 证明「分析结果不得变动」这条会失败；NC-2 证明反向守卫会失败"
      "（假登记当场被抓住）；NC-3 证明沙箱本身不产生 FAIL。三者合起来说明 V2 的四条判据"
      "都不是恒真的。")
sys.exit(0 if not fails else 1)
