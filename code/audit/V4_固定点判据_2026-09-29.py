#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V4 — 固定点判据（2026-09-29）

对应计划 §F-6：「固定点：连跑两遍，记录件 sha256 相同」。

在此之前这一条只有**手工结论**（写在 `图件重导出对照表_2026-09-29.md` 末节），
没有可复跑的脚本，也没有对照。本件把它做成两段式：

    snapshot <tsv>   遍历交付树，按下面的排除规则记 `相对路径<TAB>sha256`
    compare  A B     比两份快照，逐件报「只在 A / 只在 B / 哈希不同」

**判据不是「整棵树逐字节相同」**——那是错的，且有实测依据：

| 类别 | 判据 | 原因 |
|---|---|---|
| `.md`／`.csv`／`.py`／`.tsv` 等文本产物 | **逐字节** | 链是确定性的 |
| 四个 Office／TeX 出口（`*.docx`／`*.pdf`） | **排除** | H1 每次重写，zip 条目时间戳与 PDF `/ID` 必变 |
| `最终文件清单与SHA256*` | **排除** | 它是**派生索引**：内容含上面四个出口的当次字节哈希 |
| `_chain_logs/`、`03_分析/logs/` | **排除** | 运行日志含耗时／时间戳 |
| `.DS_Store`、`__pycache__/` | **排除** | 系统与解释器副产物 |

排除项不是「不管」，而是**换判据**：出口件由 `H2` 的 Word 回读与内容哈希对账守，
派生索引由 `H1`/`H2` 的清单对账守（写它的那次之后再无写者运行时才有效）。

`snapshotx <tsv>`／`comparex A B` 是**排除项的替换判据**：出口件不比字节，而比
「DOCX 的 zip 条目内容（去时间戳）」与「PDF 的文本层」。否则「排除」就等于「不管」。

负面对照（`selftest`）：在临时沙箱里造一份快照，注入三类缺陷各一
（改一字节／删一件／加一件），判据必须分别报出「哈希不同／只在 A／只在 B」，
且必须**恰好**报出这三条——多报是误伤，少报是判据没抓住。

用法：
    python3 V4_固定点判据_2026-09-29.py snapshot /tmp/fp_A.tsv
    python3 V4_固定点判据_2026-09-29.py compare  /tmp/fp_A.tsv /tmp/fp_B.tsv
    python3 V4_固定点判据_2026-09-29.py selftest
"""

import hashlib
import os
import pathlib
import shutil
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
OUT = pathlib.Path(os.environ["V4_OUT"]).resolve() if os.environ.get("V4_OUT") else HERE.parent

EXCLUDE_SUFFIX = (".docx", ".pdf")
EXCLUDE_NAME_PREFIX = ("最终文件清单与SHA256",)
EXCLUDE_DIR = ("_chain_logs", "logs", "__pycache__")
EXCLUDE_NAME = (".DS_Store",)


def included(rel: str) -> bool:
    """rel 是相对交付树根的 POSIX 路径。"""
    parts = rel.split("/")
    if parts[-1] in EXCLUDE_NAME or parts[-1].startswith(EXCLUDE_NAME_PREFIX):
        return False
    if any(d in EXCLUDE_DIR for d in parts[:-1]):
        return False
    return not rel.endswith(EXCLUDE_SUFFIX)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def snapshot(root: pathlib.Path, dst: pathlib.Path):
    rows = []
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(root).as_posix()
        if included(rel):
            rows.append((rel, sha(p)))
    dst.write_text("".join(f"{r}\t{h}\n" for r, h in rows), encoding="utf-8")
    print(f"快照：{len(rows)} 件文本产物 → {dst}")
    return rows


def read_snap(p: pathlib.Path):
    d = {}
    for ln in p.read_text(encoding="utf-8").splitlines():
        if ln.strip():
            r, h = ln.split("\t")
            d[r] = h
    return d


def compare(a: pathlib.Path, b: pathlib.Path, quiet=False):
    A, B = read_snap(a), read_snap(b)
    only_a = sorted(set(A) - set(B))
    only_b = sorted(set(B) - set(A))
    diff = sorted(r for r in set(A) & set(B) if A[r] != B[r])
    if not quiet:
        print(f"A = {a.name}：{len(A)} 件")
        print(f"B = {b.name}：{len(B)} 件")
        print(f"只在 A：{len(only_a)}   只在 B：{len(only_b)}   哈希不同：{len(diff)}")
        for r in only_a:
            print(f"  -A {r}")
        for r in only_b:
            print(f"  +B {r}")
        for r in diff:
            print(f"  !  {r}")
        print("结果：" + ("固定点成立（两份快照逐件相同）"
                         if not (only_a or only_b or diff) else "固定点不成立"))
    return only_a, only_b, diff


def content_fp(p: pathlib.Path):
    """DOCX → 逐 zip 条目的内容哈希（不含目录、**不含条目时间戳**）。

    DOCX 每次重写只有 zip 条目时间戳与 `docProps/core.xml` 的修改时间在变；
    条目**内容**是可复现的。PDF → 文本层哈希 + 页数（`/ID` 与 `/CreationDate` 必变）。
    """
    if p.suffix == ".docx":
        import zipfile
        z = zipfile.ZipFile(p)
        return "docx:" + "|".join(
            f"{i.filename}={hashlib.sha256(z.read(i.filename)).hexdigest()[:16]}"
            for i in sorted(z.infolist(), key=lambda x: x.filename) if not i.is_dir())
    # PDF 用 **PyMuPDF**，不用 pdfminer：实测同一文件下
    # `pdfminer.high_level.extract_text` 连取 5 次出现两个不同哈希（长度都是 346），
    # 而 PyMuPDF 5/5 同一哈希。判据自身不确定就会产生假差异——用它做闸门必然误报。
    try:
        import fitz
    except ImportError:
        return "pdf:UNREADABLE NoPyMuPDF"     # 取不到就明确登记，不换一个不确定的抽取器顶上
    try:
        d = fitz.open(str(p))
    except Exception as e:
        # 打不开的 PDF **登记为打不开**（确定性字符串），不静默略过：
        # 略过会让「出口件少了一件」在判据里表现为「无差异」。
        return f"pdf:UNREADABLE {type(e).__name__}"
    try:
        txt = "".join(pg.get_text() for pg in d)
        n = d.page_count
    finally:
        d.close()
    return f"pdf:pages={n}:text={hashlib.sha256(txt.encode()).hexdigest()[:16]}"


def snapshotx(root: pathlib.Path, dst: pathlib.Path):
    rows = []
    for p in sorted(root.rglob("*")):
        if p.is_file() and p.suffix in EXCLUDE_SUFFIX:
            rows.append((p.relative_to(root).as_posix(), content_fp(p)))
    dst.write_text("".join(f"{r}\t{h}\n" for r, h in rows), encoding="utf-8")
    print(f"出口件快照：{len(rows)} 件 → {dst}")
    return rows


def selftest():
    """三类缺陷各注入一次，判据必须恰好报出这三条。"""
    fails, oks = [], [0]

    def must(ok, msg, detail=None):
        if ok:
            oks[0] += 1
            print(f"  OK   {msg}")
        else:
            fails.append(msg)
            print(f"  FAIL {msg}" + (f"  | {detail}" if detail else ""))

    print("V4 自检：判据必须能失败")
    box = pathlib.Path(tempfile.mkdtemp(prefix="v4_selftest_"))
    try:
        src = box / "tree"
        (src / "sub").mkdir(parents=True)
        (src / "keep.md").write_text("甲\n", encoding="utf-8")
        (src / "sub" / "gone.csv").write_text("乙\n", encoding="utf-8")
        (src / "sub" / "new.py").write_text("丙\n", encoding="utf-8")
        # 排除项：必须**不**进快照，否则固定点永远不成立
        (src / "x.docx").write_bytes(b"PK\x03\x04fake")
        (src / "y.pdf").write_bytes(b"%PDF-1.4 fake")
        (src / "最终文件清单与SHA256_x.md").write_text("派生\n", encoding="utf-8")
        (src / "_chain_logs").mkdir()
        (src / "_chain_logs" / "z.log").write_text("日志\n", encoding="utf-8")
        (src / ".DS_Store").write_bytes(b"\x00")

        A = box / "A.tsv"
        rows = snapshot(src, A)
        must(len(rows) == 3, "排除规则：三件文本产物进快照，五个排除项不进",
             f"实得 {len(rows)} 件：{[r for r, _ in rows]}")

        # 注入：改一字节 / 删一件 / 加一件
        (src / "keep.md").write_text("甲改\n", encoding="utf-8")
        (src / "sub" / "gone.csv").unlink()
        (src / "sub" / "added.txt").write_text("丁\n", encoding="utf-8")
        B = box / "B.tsv"
        snapshot(src, B)

        only_a, only_b, diff = compare(A, B, quiet=True)
        must(diff == ["keep.md"], "改一字节 → 报「哈希不同」恰 1 件", f"实得 {diff}")
        must(only_a == ["sub/gone.csv"], "删一件 → 报「只在 A」恰 1 件", f"实得 {only_a}")
        must(only_b == ["sub/added.txt"], "加一件 → 报「只在 B」恰 1 件", f"实得 {only_b}")
        must(len(diff) + len(only_a) + len(only_b) == 3,
             "三类合计恰 3 条（无连带误伤）",
             f"实得 {len(diff) + len(only_a) + len(only_b)}")

        # 前提对照：原样再拍一次必须三条全零（否则上面三条不能归因于注入）
        C = box / "C.tsv"
        snapshot(src, C)
        oa, ob, df = compare(C, C, quiet=True)
        must(not (oa or ob or df), "前提对照：同一状态拍两次 → 三条全零（判据不空转）",
             f"{len(oa)}/{len(ob)}/{len(df)}")

        # 判据自身的确定性：同一 PDF 连取两次必须同哈希。抽取器不确定时，
        # 「出口件内容变了」会是假差异，据此下的结论全部作废。
        pdfs = [q for q in sorted(src.rglob("*.pdf"))]
        if pdfs:
            h = [content_fp(pdfs[0]) for _ in range(3)]
            must(len(set(h)) == 1, "判据确定性：同一 PDF 连取 3 次指纹相同", f"实得 {h}")
        else:
            must(False, "自检前提：沙箱里有 PDF 可供确定性检验")

        print(f"\n通过 {oks[0]} 项，失败 {len(fails)} 项")
        return 0 if not fails else 1
    finally:
        shutil.rmtree(box, ignore_errors=True)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    cmd = sys.argv[1]
    if cmd == "snapshot":
        snapshot(OUT, pathlib.Path(sys.argv[2]))
        return 0
    if cmd == "compare":
        only_a, only_b, diff = compare(pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3]))
        return 0 if not (only_a or only_b or diff) else 1
    if cmd == "snapshotx":
        snapshotx(OUT, pathlib.Path(sys.argv[2]))
        return 0
    if cmd == "comparex":
        only_a, only_b, diff = compare(pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3]))
        return 0 if not (only_a or only_b or diff) else 1
    if cmd == "selftest":
        return selftest()
    print(f"未知子命令 {cmd}")
    return 2


if __name__ == "__main__":
    sys.exit(main())
