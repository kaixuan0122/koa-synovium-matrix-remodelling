# -*- coding: utf-8 -*-
"""C0b —— C0 口径回源核对的负面对照。

C0 报了 1439 条断言全过。但一份「全过」的审计脚本，其价值取决于它**真的会失败**：
断言写歪、取值写错、或者破坏根本没打到脚本读的那份文件上，都会让它照样全过。
本脚本逐条把 C0 的判据弄坏，看它是否真的非零退出、且报出对得上的原因。

三处设计要点（都是踩过的坑）：
  1. **破坏必须真的到达脚本**。C0 的 BASE 是多行括号拼接（`BASE = ("…" "…")`），
     只用正则会改写失败，脚本仍读真实目录，破坏全打在副本上——8 条对照会全部"通过"。
     故这里用 sed 压成单行，并**断言**改写后脚本里 BASE 指向沙箱、且不含真实路径字面量。
  2. **对照要选真会变的单元格**。给一个本来就是 3 的值改成 3 是空操作，不是对照。
     故每条破坏都先断言确实改了、且改后的值与原值不同。
  3. 全部作业在 mktemp 目录里做，真实产物只读不写。

用法：python3 03_分析/scripts/C0b_口径核对_负面对照.py
退出码 0 表示每条对照都按预期中止。
"""
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

BASE = ("/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/"
        "Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/"
        "01.材料/01.博士材料/10.毕业/01.小论文/04.立项/"
        "滑膜基质刚度生成与刚度感知 0912/08_修订_2026-09-17")
RES = os.path.join(BASE, "02_结果表")
SCRIPT = os.path.join(BASE, "03_分析/scripts/C0_口径回源核对.py")

lines = []


def say(s=""):
    print(s)
    lines.append(s)


def rewrite_base(text: str, sandbox: str) -> str:
    """把 C0 的多行 BASE 拼接压成指向沙箱的单行。

    C0 里写的是 `BASE = ("…"` 换行续写、末行以 `")` 收尾（隐式拼接），
    不是 `BASE = (` 单独一行——按行找起止才稳，正则容易差一个字符就静默不改。
    """
    out, i, done = [], 0, False
    src = text.splitlines(keepends=True)
    while i < len(src):
        if not done and src[i].lstrip().startswith("BASE = ("):
            out.append(f'BASE = "{sandbox}"\n')
            while i < len(src) and not src[i].rstrip("\n").endswith('")'):
                i += 1
            i += 1                      # 跳过收尾那一行
            done = True
            continue
        out.append(src[i])
        i += 1
    if not done:
        raise SystemExit("改写失败：没找到 `BASE = (` 起头、`\")` 收尾的拼接——"
                         "C0 的 BASE 写法变了吗？")
    return "".join(out)


def make_sandbox():
    t = Path(tempfile.mkdtemp(prefix="c0_control_"))
    (t / "02_结果表").mkdir()
    (t / "01_审计").mkdir()
    (t / "03_分析/scripts").mkdir(parents=True)
    for f in Path(RES).glob("*.csv"):
        shutil.copy2(f, t / "02_结果表" / f.name)
    src = Path(SCRIPT).read_text(encoding="utf-8")
    rew = rewrite_base(src, str(t))
    # 自证：脚本必须只指向沙箱，否则破坏到不了脚本读的文件上
    assert f'BASE = "{t}"' in rew, "改写后 BASE 未指向沙箱"
    assert "WPS Cloud Files" not in rew, "改写后仍残留真实路径字面量"
    (t / "03_分析/scripts/C0.py").write_text(rew, encoding="utf-8")
    return t


def edit(pr: Path, fname, gene, col, new):
    """改一格。断言：确实命中一行，且新值与原值不同（否则是空操作，不是对照）。"""
    p = pr / fname
    rows = list(csv.DictReader(open(p, newline="", encoding="utf-8")))
    cols = list(rows[0])
    hit, old = 0, None
    for r in rows:
        if r["gene"] == gene:
            old = r[col]
            r[col] = new
            hit += 1
    assert hit == 1, f"对照没造对：{fname} {gene} 命中 {hit} 行"
    assert old != new, f"对照没造对：{fname} {gene}.{col} 原值就是 {new}，改了等于没改"
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    return old


def edit_const(pr: Path, fname, col, new):
    p = pr / fname
    rows = list(csv.DictReader(open(p, newline="", encoding="utf-8")))
    cols = list(rows[0])
    old = rows[0][col]
    assert old != new, f"对照没造对：{fname}.{col} 原值就是 {new}"
    for r in rows:
        r[col] = new
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


CASES = [
    ("① 发现层主估计被换成 RE-DL", "RE-DL",
     lambda pr: edit_const(pr, "合并效应_发现层.csv", "primary_method", "RE-DL")),
    ("② 总览 disc_logFC 与发现层不符", "disc_logFC",
     lambda pr: edit(pr, "57基因总览.csv", "LOXL1", "disc_logFC", "0.999")),
    ("③ 总览 rep_logFC 被填成发现层的值", "rep_logFC",
     lambda pr: edit(pr, "57基因总览.csv", "LOXL1", "rep_logFC", "0.638154383075031")),
    ("④ disc_agree_k 多记一个队列", "disc_agree_k",
     lambda pr: edit(pr, "57基因总览.csv", "ITGA2", "disc_agree_k", "2")),
    ("⑤ disc_agree_n 少记一个队列", "disc_agree_n",
     lambda pr: edit(pr, "57基因总览.csv", "ITGA2", "disc_agree_n", "2")),
    ("⑥ 该 NA 的 rep_k 被填成 3", "rep_k",
     lambda pr: edit(pr, "57基因总览.csv", "LOXL3", "rep_k", "3")),
    ("⑦ 方向列与系数符号相反", "disc_dir",
     lambda pr: edit(pr, "57基因总览.csv", "LOXL1", "disc_dir", "下调")),
    ("⑧ 复制层 k 与逐队列不符", "k=",
     lambda pr: edit(pr, "合并效应_复制层.csv", "VGLL4", "k", "2")),
    ("⑨ 复制层 primary_logFC 与逐队列不符", "VGLL4",
     lambda pr: edit(pr, "合并效应_复制层.csv", "VGLL4", "primary_logFC", "0.1")),
]


def main() -> int:
    t = make_sandbox()
    pr = t / "02_结果表"
    bro = t / "_pristine"
    shutil.copytree(pr, bro)
    say("# C0 负面对照（C0b）")
    say()
    say(f"- 沙箱：`{t}`（跑完即删；真实产物全程只读）")
    say(f"- 被检验脚本：`03_分析/scripts/C0_口径回源核对.py`（BASE 已改写指向沙箱并自证）")
    say()
    say("| 对照 | 期望中止原因 | 实测 exit | C0 报出的第一条 ❌ | 结果 |")
    say("|---|---|---|---|---|")
    bad = []
    try:
        for name, want, perturb in CASES:
            shutil.rmtree(pr)
            shutil.copytree(bro, pr)
            perturb(pr)
            r = subprocess.run([sys.executable, str(t / "03_分析/scripts/C0.py")],
                               capture_output=True, text=True)
            log = r.stdout + r.stderr
            first = next((l.strip() for l in log.splitlines() if "❌" in l), "（无 ❌ 行）")
            ok = r.returncode == 1 and want in log
            if not ok:
                bad.append(name)
            say(f"| {name} | `{want}` | {r.returncode} | {first} | {'✅' if ok else '❌'} |")
        say()
        say(f"**{len(CASES) - len(bad)}/{len(CASES)} 条对照按预期中止。**")
        if bad:
            say()
            say("未按预期中止：" + "、".join(bad) + "——C0 的对应判据形同虚设，须先修 C0。")
    finally:
        shutil.rmtree(t, ignore_errors=True)

    out = Path(BASE) / "01_审计/口径核对_负面对照_结果.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\n已写出 {out}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
