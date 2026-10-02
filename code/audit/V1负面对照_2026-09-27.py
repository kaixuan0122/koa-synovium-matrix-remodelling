#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V1 负面对照 —— 证明「稿件正文引用 → 文件定位」这道门不是空转

做法：用软链把项目根整个映射到一个沙箱（/private/tmp/v1_negctl），
其中**只有一份文件是真副本**：主稿，且把正文里第一处 `Table S4` 改成 `Table S5`。
S4 实体在、S5 实体不在，故若这道门有效，V1 必在「稿件正文引用 → 文件定位」段报 FAIL。

前提：V1 的路径全部由 __file__ 推导（HERE → OUT → BASE → REV），
故把 V1 原样复制到沙箱的 09_.../00_版本审计/ 下，它读的就是沙箱。
真文件一个字节都不动；对照跑完沙箱即弃。

自证（缺一不可，否则「对照通过」是空转的签名）：
  1. 沙箱主稿的 SHA-256 ≠ 真身主稿的 SHA-256
  2. 沙箱主稿含 `Table S5`，且 `Table S4` 出现次数比真身少 1
  3. 真身主稿的 SHA-256 == 跑之前的记录值（证明没被碰）

注：V1 的 must() 只累加失败计数、不中断（与 E1 的「失败即不出成品」不同），
故此处只断言退出码非零，不断言产物缺失。

已知沙箱伪影（与本对照要证的门无关，勿误读）：pathlib 的 rglob 不递归进符号链接目录，
而沙箱里 04_分析数据/ 等是软链，故沙箱日志里另有两道 FAIL——
「基线 206 项在磁盘上全部找得到」与「结果溯源」。这两条在真身运行时为 OK。
本对照只认目标门那行：`FAIL 正文引用的附件与文件名 … 待补齐 …`。
"""

import hashlib
import os
import pathlib
import shutil
import subprocess
import sys

BASE = pathlib.Path(
    "/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/"
    "Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/"
    "01.材料/01.博士材料/10.毕业/01.小论文/04.立项/滑膜基质刚度生成与刚度感知 0912")
SB = pathlib.Path("/private/tmp/v1_negctl")
REV_NAME = "08_修订_2026-09-17"
MS_NAME = "English_manuscript_MDPI_重定位2026-09-23.md"
V1_REL = pathlib.Path("09_Genes_CIMB_投稿修订_2026-09-27/00_版本审计/V1_版本清点_2026-09-27.py")

FAILS = []


def must(ok, msg, detail=None):
    print(f"  {'OK  ' if ok else 'FAIL'} {msg}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILS.append(msg)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


ms_real = BASE / REV_NAME / "05_稿件" / MS_NAME
must(ms_real.is_file(), f"真身主稿在位：{MS_NAME}")
sha_before = sha(ms_real)
ms_lines_real = ms_real.read_text(encoding="utf-8").split("\n")
# V1 判定的是 ms_lines[33:]（见 V1 的 `body = "\n".join(ms_lines[33:])`）。
# 首版把替换施加在**全文件**第一处 `Table S4` 上，那处落在扫描域之外，
# 于是门没响、而退出码非零是别的门造成的——对照假通过。改动必须落在扫描域内。
BODY_START = 33
body_before = "\n".join(ms_lines_real[BODY_START:])
n_s4_body = body_before.count("Table S4")
n_s4_all = "\n".join(ms_lines_real).count("Table S4")
must(n_s4_body >= 1, f"真身主稿的判定域（L{BODY_START + 1} 起）内 `Table S4` 出现 {n_s4_body} 次"
                     f"（全文件 {n_s4_all} 次）")

sb_lines = list(ms_lines_real)
for _i in range(BODY_START, len(sb_lines)):
    if "Table S4" in sb_lines[_i]:
        sb_lines[_i] = sb_lines[_i].replace("Table S4", "Table S5", 1)
        break
else:
    sys.exit("判定域内找不到 `Table S4`，负面对照无从施加")
body_after = "\n".join(sb_lines[BODY_START:])
must("Table S5" in body_after, "改写落在判定域内（判定域已含 `Table S5`）")
must(body_after.count("Table S4") == n_s4_body - 1,
     f"判定域内 `Table S4` 比真身少 1（{n_s4_body} → {body_after.count('Table S4')}）")
t_sb_body = "\n".join(sb_lines)

# ── 搭沙箱 ────────────────────────────────────────────────────────────
if SB.exists():
    shutil.rmtree(SB)
SB.mkdir(parents=True)
for p in sorted(BASE.iterdir()):
    if p.name.startswith("09_") or p.name == ".DS_Store" or p.name == REV_NAME:
        continue
    os.symlink(p, SB / p.name)
srev = SB / REV_NAME
srev.mkdir()
for p in sorted((BASE / REV_NAME).iterdir()):
    if p.name != "05_稿件":
        os.symlink(p, srev / p.name)
s05 = srev / "05_稿件"
s05.mkdir()
for p in sorted((BASE / REV_NAME / "05_稿件").iterdir()):
    if p.name == MS_NAME:
        (s05 / p.name).write_text(t_sb_body, encoding="utf-8")
    else:
        os.symlink(p, s05 / p.name)
sout = SB / V1_REL.parent
sout.mkdir(parents=True)
shutil.copy2(BASE / V1_REL, sout / V1_REL.name)

# ── 自证改写生效 ──────────────────────────────────────────────────────
ms_sb = s05 / MS_NAME
must(sha(ms_real) == sha_before, "真身主稿逐字节未动（沙箱未回写）", f"{sha_before[:12]}…")
must(sha(ms_sb) != sha_before, "沙箱主稿确实是改过的那份（哈希与真身不同）",
     f"沙箱 {sha(ms_sb)[:12]}…")
t_sb = ms_sb.read_text(encoding="utf-8")
must("Table S5" in t_sb, "沙箱主稿含改写后的 `Table S5`")
must(t_sb.count("Table S4") == n_s4_all - 1,
     f"沙箱主稿全文件的 `Table S4` 比真身少 1（{n_s4_all} → {t_sb.count('Table S4')}）")
must(not ms_sb.is_symlink(), "沙箱主稿是真文件而非软链（改写才有落点）")

# ── 跑沙箱 V1 ─────────────────────────────────────────────────────────
print("\n── 在沙箱内运行 V1")
r = subprocess.run([sys.executable, str(sout / V1_REL.name)],
                   capture_output=True, text=True, cwd=str(sout))
log = r.stdout + r.stderr
(sout / "负面对照_沙箱运行.log").write_text(log, encoding="utf-8")
print(f"  沙箱 V1 退出码 = {r.returncode}（期望 ≠ 0）")
must(r.returncode != 0, "沙箱 V1 以非零退出码结束（门响了）")

hit_seg = "── 稿件正文引用 → 文件定位" in log
must(hit_seg, "沙箱日志含「稿件正文引用 → 文件定位」段")
seg = log.split("── 稿件正文引用 → 文件定位", 1)[1].split("── ", 1)[0] if hit_seg else ""
must("Table S5" in seg and "不存在" in seg,
     "该段报出 `Table S5 → Supplementary_Table_S5.md ★不存在`（正是被改写的那处）")
# 必须匹配 **FAIL 行**：「待补齐」三个字在门通过时的 OK 行里也有
# （`OK 正文引用的附件与文件名 0 项待补齐：[]`），只查子串会误通过。
gate_fail = [l.strip() for l in log.split("\n")
             if l.strip().startswith("FAIL") and "待补齐" in l]
must(gate_fail, "沙箱日志含该门的 FAIL 行（`FAIL … 待补齐 …`）",
     detail=f"实测 {gate_fail or '无'}")
for line in log.split("\n"):
    if line.strip().startswith("FAIL") or "Table S5" in line:
        print(f"    │ {line.strip()}")

print(f"\n{'='*72}")
if FAILS:
    print(f"负面对照未达预期：{len(FAILS)} 项 —— {FAILS}")
    sys.exit(1)
print("负面对照全部通过：门在输入被改动时确实会响。")
print("真身项目未被改动一个字节（已用 SHA-256 对账）。")
shutil.rmtree(SB)
print("沙箱已清除。")
