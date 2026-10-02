#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""B10 修订稿词数（R1.2）

用途：重算修订稿 Title Page 自报的正文与摘要词数。

**口径不在此脚本内重写**——直接 import 原脚本 `04_分析数据/scripts/17_word_count.py`
并只替换其 DOCS 路径，故分词、去引用标记、去 markdown 的规则与原脚本逐字相同。
（若复制一份逻辑过来，两处会各自漂移，词数就不再可比。）

用法：python3 scripts/B10_修订稿词数.py
输出：stdout；退出码 0 表示与稿件 Title Page 自报值一致
"""
import importlib.util
import re
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]          # …/08_修订_2026-09-17
PROJ = BASE.parent                                   # …/滑膜基质刚度生成与刚度感知 0912
SRC = PROJ / "04_分析数据" / "scripts" / "17_word_count.py"
EN = BASE / "05_稿件" / "English_manuscript_MDPI_重定位2026-09-23.md"
CN = BASE / "05_稿件" / "双语初稿_候选优选与组织来源_重定位2026-09-23.md"

spec = importlib.util.spec_from_file_location("wc17", SRC)
wc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wc)

res = []


def chk(name, got, want, tol=0):
    ok = abs(int(got) - int(want)) <= tol
    res.append((ok, name, got, want))


def main():
    en = EN.read_text(encoding="utf-8")
    cn = CN.read_text(encoding="utf-8")

    ab_en = wc.count(wc.abstract_section(en))
    mn_en = wc.count(wc.main_section(en))
    ab_cn = wc.count(cn.split("## 摘要（Abstract）", 1)[1]
                     .split("**关键词**", 1)[0])

    # 稿件 Title Page 自报值：`**Word count:** abstract 200 words; main text … 7261 words.`
    line = [l for l in en.splitlines() if l.startswith("**Word count:**")]
    assert len(line) == 1, f"Title Page 的 Word count 行数异常：{len(line)}"
    m = re.search(r"abstract (\d+) words; main text .*?(\d+) words", line[0])
    chk("Title Page 自报摘要词数", ab_en, int(m.group(1)))
    chk("Title Page 自报正文词数", mn_en, int(m.group(2)))

    print(f"英文稿  摘要 {ab_en} 词；正文 {mn_en} 词")
    print(f"中文稿  摘要 {ab_cn} 词（原脚本口径：中文正文不参与计数）")

    npass = sum(1 for r in res if r[0])
    for ok, name, got, want in res:
        if not ok:
            print(f"[不通过] {name}: 实测={got} 稿件自报={want}")
    print(f"\n共 {len(res)} 项，通过 {npass} 项，不通过 {len(res) - npass} 项")
    return 0 if npass == len(res) else 1


if __name__ == "__main__":
    sys.exit(main())
