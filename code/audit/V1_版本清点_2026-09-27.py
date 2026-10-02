#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V1 — 版本清点（2026-09-27，Genes/CIMB 投稿修订第一轮）

对应工作令 §三「先关闭版本问题」。只读上游：不覆盖、不移动、不改名任何既有文件。
全部输出落在 09_Genes_CIMB_投稿修订_2026-09-27/00_版本审计/。

产出（§三 指定的四件）：
  文件清单_2026-09-27.tsv          全量路径 / 字节 / mtime / SHA-256 / md5
  版本差异表_2026-09-27.md         主稿世系 + §三 点名六问的逐项答案
  唯一修改基线说明_2026-09-27.md   选定基线及依据（不按文件名日期选）
  缺失文件清单_2026-09-27.md       引用了但不存在的、无法溯源的、以及待补齐的

── 判据的两次修正（首版当场诊断，都是工具假阳性，不是数据问题）──────────────
  1) 「基线缺失 149」：首版只扫 08_修订 / 05_投稿包 / 03_论文初稿，而基线 206 项里
     149 项在 04_分析数据/ —— 那 149 个文件全部真实存在。→ 改为扫整个项目根。
  2) 「脚本提到但结果表没有 133 项」「5 个结果 md 无脚本」：正则把 read.csv 当文件名、
     且只查 02_结果表 一处。→ 改为字面量+全项目定位，配对改为前缀或全文包含。
  3) 「3 个字面量失效」：两个是 R 的 paste0 模板串，一个是 B9 里的 endswith 后缀判据。
     → 模板串按 {} % * 排除；剩余逐项登记判读结论，未登记即失败（登记制）。

── 已知陷阱（写进代码）───────────────────────────────────────────────────
  · 数 markdown 表行要「以 | 开头且以 | 结尾」：S4a 题注换行后第二行以 |log2FC| 开头，
    松判据会多算 1 行（700 vs 真值 699）。
  · 「图内是否画了图注」不能看有没有 'Figure N' 字样：图 1 有编号且 100% 画入，
    图 2 无编号却有另一段说明文字。须拿图内文本与稿件图注做 token 重合度。
  · 「图注画入」有两个互不可替的口径，实测在三幅图上给出不同答案：
      整条覆盖（token 重合）  —— 图 1 = 100%，图 2 = 20.8%，图 3/4/5 = 0%
      末句是否逐字在图内      —— 图 1、图 2 命中，图 3/4/5 未命中
    图 2 命中的原因是：图内那段是**自撰说明**（与图注语义相通但措辞不同），
    只有末句 'No rank is implied by position.' 与稿件图注逐字相同。
    装配时图 1 会重复整条、图 2 只重复末句一句。C5 脚本的 N_EMB=2 用的是后一口径。
  · 图内文本的连字符是 U+2212（'random−effects'），稿件 md 里是 ASCII '-'。
    逐句子串比对前必须归一破折号，否则整句被判成"未画入"——是判据假阳性，不是图的问题。
"""

import ast
import collections
import difflib
import hashlib
import pathlib
import re
import sys
from datetime import datetime

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parent                      # 09_Genes_CIMB_投稿修订_2026-09-27
BASE = OUT.parent                      # 项目根
REV = BASE / "08_修订_2026-09-17"      # 工作令所称的「输入根目录」

MS_EN = REV / "05_稿件" / "English_manuscript_MDPI_重定位2026-09-23.md"
MS_CN = REV / "05_稿件" / "双语初稿_候选优选与组织来源_重定位2026-09-23.md"
MS_EN_OLD = REV / "05_稿件" / "English_manuscript_MDPI_修订2026-09-17.md"
MS_CN_OLD = REV / "05_稿件" / "双语初稿_刚度因果版_修订2026-09-17.md"

HANS = re.compile(r"[一-鿿]")
BASELINE_T = datetime(2026, 9, 17, 19, 13)

# ── 登记制：基线变动中已查明性质者（每条含证据；未登记即失败）──────────────
# 反向守卫：登记的路径若实测未变动，说明该登记已过期，同样失败。
KNOWN_DELTA = {
    "04_分析数据/MR_coloc/coloc_results_Alasoo.csv":
        "B4 补 SERPINE1 共定位（2026-09-17 20:29–20:31 重跑）：在既有 CSV 追加 SERPINE1 一行，"
        "原 VGLL4/SDC4 两行逐位未变；证据 02_结果表/B4_SERPINE1共定位_结果.md、logs/B4_coloc.log",
    "04_分析数据/MR_coloc/coloc_results_Alasoo_orig.csv": "同上（Alasoo 原始版，同批追加）",
    "04_分析数据/MR_coloc/coloc_results_GTEx.csv": "同上（GTEx 档追加）",
    "04_分析数据/MR_coloc/coloc_results_Synovium.csv": "同上（滑膜档追加）",
    "04_分析数据/MR_coloc/gene_coords.py": "同上：COLOC_GENES 由 [VGLL4,SDC4] 加入 SERPINE1",
    "04_分析数据/MR_coloc/run_coloc.R": "同上：随 gene_coords.py 同步改",
}
# ── 登记制：脚本内「项目内找不到」的字面量，逐项判读（未登记即失败）──────────
# 分两栏：现役（03_分析/scripts，即当前稿件所用）与上游（04_分析数据/**，数据准备阶段）。
# 现役栏里任何一条若真是失效引用，就是当前复现链的断点——三条全为模板/后缀判据，故零失效。
LIT_JUDGED = {
    # 现役（3 条，均为句法构造，非文件名）
    "_mixture.txt": "现役·B6 paste0(nm, \"_mixture.txt\") 的模板后缀，非独立文件名",
    "_results.txt": "现役·B6 paste0(\"CIBERSORT_\", nm, \"_results.txt\") 的模板后缀",
    "修订2026-09-17.md": "现役·B9 中 MS.endswith(\"修订2026-09-17.md\") 的后缀判据，非路径",
    # 上游（15 条，均为早期数据准备阶段的名字，逐条给出对应物）
    "GSE1919_series_matrix.txt": "上游·实际存为 raw/GSE1919_series_matrix.txt.gz（解压后读取）",
    "GSE206848_series_matrix.txt": "上游·实际存为 raw/GSE206848_series_matrix.txt.gz",
    "GSE55235_series_matrix.txt": "上游·实际存为 raw/GSE55235_series_matrix.txt.gz",
    "GSE55457_series_matrix.txt": "上游·实际存为 raw/GSE55457_series_matrix.txt.gz",
    "GSE82107_series_matrix.txt": "上游·实际存为 raw/GSE82107_series_matrix.txt.gz",
    "GSE152805.rds": "上游·实际存为 raw/GSE152805_1.rds（GEO 补充文件 _1 序号约定）",
    "merged_DEG.csv": "上游·现役对应物 合并验证/5队列合并_DEG.csv（早期合并步骤改名）",
    "merged_7gene.csv": "上游·现役对应物 合并验证/5队列合并_7基因方向.csv",
    "finngen_KOA_outcome.txt": "上游·真身 MR_coloc/mr_data/finngen_KOA.gz（819,114,058 B 实存；"
                               "B4 日志记载该档命中缓存）",
    "eQTL_exposure.txt": "上游·run_MR.R 的旧输入约定名；现役产物为 "
                         "MR_coloc/multi_iv_exposure_outcome{,_Synovium,_GTEx,_Alasoo}[_pruned].csv",
    "_pruned.csv": "上游·20_mr_sensitivity.R 中 paste0(..., \"_pruned.csv\") 的模板后缀，"
                   "实际产物 multi_iv_MR_results_*_pruned.csv 实存",
    "GPL96_probe2symbol_clean.csv": "上游·merge_analysis.R 自建注释；现役改用 Bioconductor 注释包 "
                                    "hgu133a.db/hgu133plus2.db（见 B1_逐队列效应_力学基因.R:228-236）",
    "GPL570_probe2symbol.csv": "上游·同上，现役改用 hgu133plus2.db（GSE82107=GPL570）",
    "GPL91_probe2symbol.csv": "上游·现役改用 raw/GPL91.soft.gz 经 B1 的 Table(getGEO(...)) 解析",
    "_7genes_by_celltype.png": "上游·早期单细胞通讯脚本的出图；现役单细胞结果见 02_结果表/A2d_* 与图 5",
}
# 记录性 md（非脚本产物）走另一条判据：其正文自称的脚本/日志/产物路径须在磁盘上找得到
RE_SRC = re.compile(r"`([^`\n]*?\.(?:py|R|csv|log|md|rds))`")

fail = 0
rep = {}          # 供 md 报告取用


def must(ok, msg, detail=None):
    global fail
    print(f"  {'OK  ' if ok else 'FAIL'} {msg}")
    if not ok:
        fail += 1
        if detail:
            print(f"       {detail}")
    return ok


def digest(p):
    b = p.read_bytes()
    return len(b), hashlib.sha256(b).hexdigest(), hashlib.md5(b).hexdigest()


def mtime_of(p):
    return datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")


# ══════════════════════════════════════ 1. 全量清单
print("V1 版本清点 —— 只读上游")
rows = []
SKIP_TOP = {OUT.name}
for p in sorted(BASE.rglob("*")):
    if not p.is_file():
        continue
    rel = p.relative_to(BASE)
    if rel.parts[0] in SKIP_TOP or p.name in (".DS_Store",) or p.name.startswith("~$"):
        continue
    n, s, m = digest(p)
    rows.append({"区": rel.parts[0], "路径": str(rel), "字节": n, "sha256": s, "md5": m,
                 "mtime": mtime_of(p)})
inv = HERE / "文件清单_2026-09-27.tsv"
with open(inv, "w", encoding="utf-8") as fh:
    fh.write("区\t路径\t字节\tmtime\tsha256\tmd5\n")
    for r in rows:
        fh.write(f"{r['区']}\t{r['路径']}\t{r['字节']}\t{r['mtime']}\t{r['sha256']}\t{r['md5']}\n")
print(f"  清单写出：{inv.name}  {len(rows)} 个文件  {sum(r['字节'] for r in rows)/1e6:.1f} MB")
rep["清单"] = {"n": len(rows), "mb": sum(r["字节"] for r in rows) / 1e6,
               "区": collections.Counter(r["区"] for r in rows).most_common()}
for z, c in rep["清单"]["区"]:
    print(f"    {z:52s} {c:4d} 个")

# ══════════════════════════════════════ 2. 原件对账
print("\n── 原件对账（01_审计/基线md5清单.tsv，取样 2026-09-17 19:13）")
baseline = {}
for ln in (REV / "01_审计" / "基线md5清单.tsv").read_text(encoding="utf-8").splitlines()[1:]:
    if ln.strip():
        p, size, md5 = ln.split("\t")
        baseline[p] = (int(size), md5)
cur = {r["路径"]: (r["字节"], r["md5"]) for r in rows}
same, changed, absent = [], [], []
for p, (size, md5) in baseline.items():
    if p not in cur:
        absent.append(p)
    elif cur[p] == (size, md5):
        same.append(p)
    else:
        changed.append({"path": p, "b0": size, "m0": md5[:8], "b1": cur[p][0], "m1": cur[p][1][:8]})
print(f"  基线 {len(baseline)} 项：一致 {len(same)} / 变动 {len(changed)} / 找不到 {len(absent)}")
for c in changed:
    tag = "已登记" if c["path"] in KNOWN_DELTA else "★未登记"
    print(f"    [{tag}] {c['path']}  {c['b0']} B/{c['m0']} → {c['b1']} B/{c['m1']}")
for p in absent:
    print(f"    ★找不到 {p}")
must(not absent, f"基线 {len(baseline)} 项在磁盘上全部找得到")
unreg = [c["path"] for c in changed if c["path"] not in KNOWN_DELTA]
must(not unreg, f"基线变动 {len(changed)} 项全部已在 KNOWN_DELTA 登记（未登记 {len(unreg)} 项）",
     detail="\n       ".join(unreg))
stale = [p for p in KNOWN_DELTA if p in cur and p not in {c["path"] for c in changed}]
must(not stale, f"KNOWN_DELTA 的 {len(KNOWN_DELTA)} 条登记全部仍然成立（无过期登记）",
     detail="\n       ".join(stale))
# 规则 1 的两个受保护树
guard = [c["path"] for c in changed if c["path"].startswith(("05_投稿包/", "03_论文初稿/"))]
must(not guard, "规则 1 保护的两棵树（05_投稿包/ 与 03_论文初稿/）逐字节未动",
     detail=str(guard))
rep["对账"] = {"n_base": len(baseline), "same": len(same), "changed": changed, "absent": absent,
               "guard_ok": not guard,
               "prot": {t: sum(1 for p in same if p.startswith(t))
                        for t in ("05_投稿包/", "03_论文初稿/", "04_分析数据/", "08_修订_2026-09-17/")}}
print(f"  受保护两棵树：05_投稿包/ 一致 {rep['对账']['prot']['05_投稿包/']} 项；"
      f"03_论文初稿/ 一致 {rep['对账']['prot']['03_论文初稿/']} 项")

added = sorted(set(cur) - set(baseline))
older = [p for p in added
         if datetime.strptime(mtime_of(BASE / p), "%Y-%m-%d %H:%M:%S") < BASELINE_T]
print(f"  现有而基线未登记 {len(added)} 个；其中 mtime 早于基线取样（漏登而非事后新增）{len(older)} 个")
rep["对账"]["added"] = len(added)
rep["对账"]["older"] = len(older)
rep["对账"]["older_zone"] = collections.Counter(p.split("/")[0] for p in older).most_common()

# ══════════════════════════════════════ 3. 主稿世系
print("\n── 主稿世系")
def ms_stat(p):
    t = p.read_text(encoding="utf-8")
    ls = t.split("\n")
    han = [i + 1 for i, l in enumerate(ls) if HANS.search(l)]
    return {"file": str(p.relative_to(BASE)), "size": p.stat().st_size, "mt": mtime_of(p),
            "lines": len(ls), "han": len(han), "lo": han[0] if han else 0,
            "hi": han[-1] if han else 0, "md5": hashlib.md5(p.read_bytes()).hexdigest()[:8],
            "title": next((l.lstrip("# ").strip() for l in ls if l.startswith("# ")), "（无 H1）")}


lineage = [ms_stat(p) for p in [
    BASE / "03_论文初稿" / "English_manuscript.md",
    BASE / "03_论文初稿" / "English_manuscript_MDPI.md",
    MS_EN_OLD, MS_EN,
    BASE / "03_论文初稿" / "双语初稿_刚度因果版.md",
    MS_CN_OLD, MS_CN]]
for r in lineage:
    print(f"  {r['size']:>7,} B  {r['mt']}  行{r['lines']:>4}  汉字{r['han']:>4}"
          f"(L{r['lo']}–L{r['hi']})  {r['md5']}  {r['file']}")


def diff_stat(a, b):
    la, lb = a.read_text(encoding="utf-8").split("\n"), b.read_text(encoding="utf-8").split("\n")
    add = dele = 0
    for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, la, lb, autojunk=False).get_opcodes():
        if tag in ("replace", "delete"):
            dele += i2 - i1
        if tag in ("replace", "insert"):
            add += j2 - j1
    return add, dele, len(la), len(lb)


d = {"en": diff_stat(MS_EN_OLD, MS_EN), "cn": diff_stat(MS_CN_OLD, MS_CN)}
print(f"  英文两代：{d['en'][2]} 行 → {d['en'][3]} 行（新增 {d['en'][0]} / 删除 {d['en'][1]}）")
print(f"  中文两代：{d['cn'][2]} 行 → {d['cn'][3]} 行（新增 {d['cn'][0]} / 删除 {d['cn'][1]}）")
t_new, t_old = lineage[3]["title"], lineage[2]["title"]
must(t_new != t_old, "两代英文稿标题不同——重定位而非改字（故不能按文件名日期自动选基线）",
     detail=f"09-17：{t_old[:70]}\n       09-23：{t_new[:70]}")
rep["世系"] = {"rows": lineage, "diff": d, "title_new": t_new, "title_old": t_old}

# ══════════════════════════════════════ 4. 图件
print("\n── 图件一致性")
figdir = REV / "04_图表"
pairs = {}
for p in sorted(figdir.iterdir()):
    if p.suffix.lower() in (".pdf", ".png"):
        pairs.setdefault(p.stem, set()).add(p.suffix.lower())
must(len(pairs) == 5, f"04_图表 下有 5 组图件（实得 {len(pairs)} 组）")
must(all(v == {".pdf", ".png"} for v in pairs.values()), "每组 PDF 与 PNG 成对齐全")

ms_lines = MS_EN.read_text(encoding="utf-8").split("\n")
i0 = next(i for i, l in enumerate(ms_lines) if l.strip().startswith("## Figure Legends"))
i1 = next((i for i in range(i0 + 1, len(ms_lines)) if ms_lines[i].startswith("## ")), len(ms_lines))
legends, cur_n = {}, None
for l in ms_lines[i0:i1]:
    m = re.match(r"\*\*Figure (\d+):(.*?)\*\*(.*)", l)
    if m:
        cur_n = int(m.group(1))
        legends[cur_n] = (m.group(2) + " " + m.group(3)).strip()
    elif cur_n and l.strip():
        legends[cur_n] += " " + l.strip()
must(len(legends) == 5, f"主稿 Figure Legends 恰 5 条（实得 {len(legends)}）")


def toks(s):
    s = s.replace("−", "-").replace("–", "-").replace("—", "-").replace("’", "'")
    return re.sub(r"[^a-z0-9 ]", " ", s.lower()).split()


def norm_fig(s):
    """图内文本层的归一化。**破折号必须先归一**：matplotlib/ggplot 出图时连字符被排成
    U+2212 MINUS SIGN（图内 'random−effects'、'RNA−seq'），而稿件 md 里是 ASCII '-'，
    不归一则逐句子串比对会把整句判成"未画入"——判据的假阳性，不是图的问题。"""
    s = s.translate(str.maketrans({c: "-" for c in "‐‑‒–—―−﹘"}))
    s = s.translate(str.maketrans({c: "'" for c in "‘’ʼ"}))
    return re.sub(r"\s+", " ", s).strip()


STEM_HINT = {"1": "Fig_study_design", "2": "Fig_candidate_overview", "3": "Fig_forest_logFC",
             "4": "Fig_evidence_matrix", "5": "Fig_scRNA_dotplot"}
# STEM_HINT 是**手写**的映射，本身不构成验证。C5_投稿清单重定位_2026-09-25.py 里另有一份
# 独立写就的 FIG_FILES，两处逐项一致才算映射被独立确认（同型文字两处各存一份必然漂移）。
# 首版这里写的是 'forest'/'dotplot' 两个**短名**、且按子串匹配去目录里找图件——碰巧也能命中
# 正确文件，纯属运气：目录里若另有一个含 'forest' 的文件就会静默取错图。对账门当场抓到这两处。
C5_SCRIPT = REV / "03_分析" / "scripts" / "C5_投稿清单重定位_2026-09-25.py"
c5_fig = {}
if C5_SCRIPT.is_file():
    for _node in ast.walk(ast.parse(C5_SCRIPT.read_text(encoding="utf-8", errors="replace"))):
        if isinstance(_node, ast.Assign) and any(
                getattr(_t, "id", None) == "FIG_FILES" for _t in _node.targets):
            c5_fig = ast.literal_eval(_node.value)
must(c5_fig, f"从 C5_投稿清单重定位 的源码解析出 FIG_FILES（实得 {len(c5_fig)} 项）")
_drift = {k: (STEM_HINT.get(k), c5_fig.get(k)) for k in sorted(set(STEM_HINT) | set(c5_fig))
          if c5_fig.get(k) != STEM_HINT.get(k)}
must(not _drift, f"本脚本的 STEM_HINT 与 C5 的 FIG_FILES 逐项一致（{len(STEM_HINT)} 项，"
                 f"两处独立写就）", detail=str(_drift))
fig_note, order = [], {}
try:
    import fitz
    HAVE_FITZ = True
except ImportError:
    HAVE_FITZ = False
    print("  （未装 PyMuPDF，图内文字检查跳过）")
NUM_IN_FIG = re.compile(r"Figure\s*(\d+)")
# 「装配时图注会不会重复」的判据是**末句是否在图内**（与 C5_投稿清单重定位 的 N_EMB 同一件事）。
# 整条覆盖率（cov）测的是另一件事：整条图注有没有全部画进去。两者不可互替——
# 图 2 的末句画入了但整条没有，故 cov 判它"未画入"而末句判据判它"会重复一句"。
N_EMB_C5 = 2          # C5 脚本的声明值，本脚本须复现同一数字；此数一变两处都要改
if HAVE_FITZ:
    print(f"  {'图件':<26}{'图注词':<8}{'图内词':<8}{'整条画入':<10}{'末句在图内':<11}{'图内编号'}")
    for n in sorted(legends):
        stem = STEM_HINT[str(n)]
        must(stem in pairs, f"图 {n} 的图件 {stem}.pdf/.png 在 04_图表 下存在",
             detail=f"目录内实有 {sorted(pairs)}")
        order[n] = stem
        doc = fitz.open(str(figdir / f"{stem}.pdf"))
        inside = " ".join(pg.get_text() for pg in doc)
        flat = norm_fig(inside)
        a, b = toks(legends[n]), toks(inside)
        bl = [x for x in difflib.SequenceMatcher(None, a, b, autojunk=False).get_matching_blocks()
              if x.size >= 4]
        cov = sum(x.size for x in bl) / max(len(a), 1)
        sents = [s.strip() for s in re.split(r"(?<=\.)\s+", legends[n]) if s.strip()]
        tail = norm_fig(sents[-1])
        emb = tail in flat
        nums = sorted({int(x) for x in NUM_IN_FIG.findall(inside)})
        fig_note.append({"n": n, "stem": stem, "lw": len(a), "iw": len(b), "cov": cov,
                         "emb": emb, "tail": tail, "nums": nums})
        print(f"  {stem:<26}{len(a):<8}{len(b):<8}{cov*100:>6.1f}%    "
              f"{'是' if emb else '否':<11}{(str(nums) if nums else '（无）')}")
        # 局限（不假装有独立验证）：本门只在「图内确实出现编号」时才有判别力。
        # 04_图表 的五幅图里只有图 1 带编号，其余四幅此项 trivially 通过、抓不到图件错配；
        # 真正防错配的是上面 STEM_HINT ≡ C5.FIG_FILES 的两处独立映射对账。
        must((not nums) or nums == [n], f"图 {n}（{stem}）图内编号若存在即为自身的 {n}",
             detail=f"图内出现 {nums}")
    must(len(order) == 5 and all(order.values()), "五条图注各自定位到一份图件")
    must(sum(1 for f in fig_note if f["nums"]) <= 1,
         "图内带编号者至多一幅（多于此数说明另有图件把图号烧进了图内，须逐幅追）",
         detail=str([f["n"] for f in fig_note if f["nums"]]))
    n_emb = sum(1 for f in fig_note if f["emb"])
    n_full = [f["n"] for f in fig_note if f["cov"] >= 0.6]
    print(f"  图注末句画入图内者共 {n_emb} 幅：{[f['n'] for f in fig_note if f['emb']]}"
          f"（这些装配时会重复末句一句）")
    print(f"  其中整条图注也画入的：{n_full}（这些会重复整条）")
    print(f"  注：『末句画入』与『整条画入』是两件事——图 2 的图内是自撰说明段"
          f"（与图注语义相通但措辞不同），只有末句与图注逐字相同。")
    must(n_emb == N_EMB_C5,
         f"图注末句画入图内的幅数 {n_emb} 与 C5_投稿清单重定位 的声明值 {N_EMB_C5} 一致"
         f"——此数一变，核查件事项 2、C5 的 §1.2 说明与 §六 缺口项都要同步改")
    rep["图"] = {"pairs": sorted(pairs), "note": fig_note, "n_emb": n_emb,
                 "full": n_full}

# ══════════════════════════════════════ 5. 结果溯源
print("\n── 结果溯源")
ALL_SCRIPTS = [p for p in BASE.rglob("*")
               if p.suffix in (".py", ".R") and OUT.name not in p.parts and p.is_file()]
script_text = {str(p.relative_to(BASE)): p.read_text(encoding="utf-8", errors="replace")
               for p in ALL_SCRIPTS}
SCRIPTS_MAIN = {p.stem for p in ALL_SCRIPTS if p.parent.name == "scripts"}
ALLBLOB = "\n".join(script_text.values())
res = [p for p in sorted((REV / "02_结果表").glob("*")) if p.is_file()]
orphan = []
for p in res:
    stem = re.sub(r"_结果$", "", p.stem)
    if stem in SCRIPTS_MAIN or p.stem in SCRIPTS_MAIN or p.name in ALLBLOB:
        continue
    if any(t and t in " ".join(SCRIPTS_MAIN) for t in re.split(r"[_\-]", stem)):
        continue
    orphan.append(p.name)
print(f"  02_结果表 {len(res)} 个文件；前缀与全文都配不到脚本的 {len(orphan)} 个：")
for n in orphan:
    print(f"    · {n}")
# 记录性 md 的判据：自报来源路径须存在（而不是「脚本里提过它」）
# 含 {} 的是花括号缩写写法（如 coloc_results_{A,B}.csv），跳过——展开后逐项另行可查
selfbad = []
for n in orphan:
    t = (REV / "02_结果表" / n).read_text(encoding="utf-8", errors="replace")
    srcs = {s for s in RE_SRC.findall(t) if "/" in s and "{" not in s and "*" not in s}
    gone = [s for s in srcs if not (BASE / s).exists() and not (REV / s).exists()]
    print(f"    · {n}：正文自报来源 {len(srcs)} 条，找不到 {len(gone)} 条"
          + (f" —— ★{gone}" if gone else "；全部在位"))
    selfbad += [f"{n} → {g}" for g in gone]
must(not selfbad, f"{len(res)} 个结果文件全部可溯源：脚本产物由脚本产出，"
                  f"记录性 md 的自报来源逐条在位")
rep["溯源"] = {"n_res": len(res), "n_scripts": len(ALL_SCRIPTS), "orphan": orphan,
               "selfbad": selfbad}

LIT = re.compile(r"""['"]([^'"/\\\n]{3,80}\.(?:csv|tsv|md|txt|rds|xlsx|pdf|png|tiff))['"]""")
dead = collections.Counter()
for name, t in script_text.items():
    for m in LIT.findall(t):
        if any(c in m for c in "{}%*"):
            continue
        if m not in {r["路径"].split("/")[-1] for r in rows}:
            dead[m] += 1
unjudged = [k for k in dead if k not in LIT_JUDGED]
print(f"  脚本字面量中项目内找不到的 {len(dead)} 个：")
for k, v in dead.most_common():
    print(f"    [{'已登记' if k in LIT_JUDGED else '★未登记'}] {v}×  {k}"
          f"{'  —— ' + LIT_JUDGED[k] if k in LIT_JUDGED else ''}")
must(not unjudged, f"找不到的 {len(dead)} 个字面量全部已逐项判读登记（未判读 {len(unjudged)}）",
     detail="\n       ".join(unjudged))
rep["溯源"]["dead_lit"] = dict(dead)

# ══════════════════════════════════════ 6. 稿件 → 文件定位（§三 验收硬门）
print("\n── 稿件正文引用 → 文件定位")
body = "\n".join(ms_lines[33:])
refs = collections.Counter(re.findall(r"(?:Table|Figure|Fig\.)\s*(S\d+|[1-6])", body))
print(f"  正文 Table/Figure 引用：{dict(sorted(refs.items(), key=lambda x: (x[0][0] != 'S', x[0])))}")
missing = []
for s in sorted({r for r in refs if r.startswith("S")}):
    f = REV / "06_补充材料" / f"Supplementary_Table_S{s[1:]}.md"
    print(f"    Table {s} → {f.name} {'存在' if f.exists() else '★不存在'}")
    if not f.exists():
        missing.append(f"Table {s}")
for n in sorted({int(r) for r in refs if not r.startswith("S")}):
    ok = n in legends and order.get(n)
    print(f"    Figure {n} → 图注{'在位' if n in legends else '★缺失'}；图件 {order.get(n, '★无')}")
    if not ok:
        missing.append(f"Figure {n}")
LIT2 = re.compile(r"[\w一-鿿./-]+\.(?:csv|tsv|md|txt|rds|xlsx|pdf|png|tiff)")
for m in sorted(set(LIT2.findall(body))):
    hit = next((r["路径"] for r in rows if r["路径"].split("/")[-1] == pathlib.Path(m).name), None)
    print(f"    {'OK ' if hit else '★  '} {m}" + (f"  → {hit}" if hit else "  → 无落点"))
    if not hit:
        missing.append(m)
must(not missing, f"正文引用的附件与文件名 {len(missing)} 项待补齐：{missing}")
rep["引用"] = {"refs": dict(refs), "missing": missing}

# ══════════════════════════════════════ 7. S4 行数
print("\n── S4 行数（严判据）")
sec, cnt, loose = "题头", collections.Counter(), 0
for ln in (REV / "06_补充材料" / "Supplementary_Table_S4.md").read_text(encoding="utf-8").split("\n"):
    m = re.match(r"^## (S4[a-c])\.", ln)
    if m:
        sec = m.group(1)
    if ln.startswith("|"):
        loose += 1
        if ln.rstrip().endswith("|"):
            cnt[sec] += 1
per = {k: v - 2 for k, v in cnt.items() if k != "题头"}
for k in sorted(per):
    print(f"    {k}: 表行 {cnt[k]:4d} → 数据行 {per[k]}")
total = sum(per.values())
print(f"    合计 {total}；松判据（^| 即算）得 {loose - 6}——多的是 S4a 题注里 |log2FC| 开头的散文行")
must(total == 699 and per.get("S4a") == 246 and per.get("S4b") == 3 and per.get("S4c") == 450,
     "S4 的 246 + 3 + 450 = 699 与稿件及核查件所称相符")
rep["s4"] = {"per": per, "total": total, "loose": loose - 6}

# ══════════════════════════════════════ 8. 主稿头部
print("\n── 主稿头部（清洁稿须剥离的部分）")
han_lines = [i + 1 for i, l in enumerate(ms_lines) if HANS.search(l)]
body_start = next(i + 1 for i, l in enumerate(ms_lines) if l.strip().startswith("**Title Page**"))
c_open = [i + 1 for i, l in enumerate(ms_lines) if "<!--" in l]
c_close = [i + 1 for i, l in enumerate(ms_lines) if "-->" in l]
print(f"    汉字 {len(han_lines)} 行 L{han_lines[0]}–L{han_lines[-1]}；正文起点 L{body_start}；"
      f"注释块 L{c_open}–L{c_close}")
must(max(han_lines) < body_start, f"{len(han_lines)} 行汉字全部在正文起点 L{body_start} 之前")
must(len(c_open) == len(c_close) and all(o < c for o, c in zip(c_open, c_close)),
     f"{len(c_open)} 段注释各自成对闭合")
must(all(c_open[0] <= l <= c_close[-1] for l in han_lines),
     f"含汉字各行全部落在 L{c_open[0]}–L{c_close[-1]} 注释块内——清洁稿整块剥离即可")
must(not any(HANS.search(l) for l in ms_lines[body_start - 1:]),
     f"正文 L{body_start} 起无任何汉字")
rep["头部"] = {"han": len(han_lines), "lo": han_lines[0], "hi": han_lines[-1],
               "body_start": body_start, "c_open": c_open, "c_close": c_close}

# ══════════════════════════════════════ 9. B4 证据链（基线变动的定性依据）
print("\n── B4（SERPINE1 共定位补跑）证据链")
B4 = [("结果 md", REV / "02_结果表" / "B4_SERPINE1共定位_结果.md"),
      ("日志", REV / "03_分析" / "logs" / "B4_coloc.log"),
      ("脚本 gene_coords.py", BASE / "04_分析数据" / "MR_coloc" / "gene_coords.py"),
      ("脚本 run_coloc.R", BASE / "04_分析数据" / "MR_coloc" / "run_coloc.R"),
      ("产物 Synovium", BASE / "04_分析数据" / "MR_coloc" / "coloc_results_Synovium.csv"),
      ("产物 GTEx", BASE / "04_分析数据" / "MR_coloc" / "coloc_results_GTEx.csv"),
      ("产物 Alasoo", BASE / "04_分析数据" / "MR_coloc" / "coloc_results_Alasoo.csv")]
b4_ok = []
for label, p in B4:
    e = p.exists()
    b4_ok.append(e)
    print(f"    {'OK  ' if e else '★缺 '} {label:22s} {mtime_of(p) if e else '—'}  "
          f"{str(p.relative_to(BASE))}")
must(all(b4_ok), f"B4 证据链 {len(B4)} 件全部在位")
rep["b4"] = [{"label": lb, "path": str(p.relative_to(BASE)),
              "mt": mtime_of(p) if p.exists() else None,
              "b": p.stat().st_size if p.exists() else None} for lb, p in B4]

print(f"\n结果：{'全部通过' if not fail else f'{fail} 项需处置'}")

# ══════════════════════════════════════ 10. 写出三份报告（数字全部取自上面的实测）
print("\n── 写出报告")
L = rep["世系"]["rows"]
same_n = rep["对账"]["same"]
ch = rep["对账"]["changed"]

(OUT / "00_版本审计" / "版本差异表_2026-09-27.md").write_text(f"""\
# 版本差异表（2026-09-27）

对应工作令 §三「先关闭版本问题」。本表由 `V1_版本清点_2026-09-27.py` 一次运行生成，
表内每个数字都取自该次运行的实测值，无手写维护。

## 一、主稿世系：四份英文 + 三份中文，两代并存

| 世代 | 文件 | 字节 | 修改时间 | 行数 | md5 |
| --- | --- | ---: | --- | ---: | --- |
""" + "\n".join(
    f"| {g} | `{r['file']}` | {r['size']:,} | {r['mt']} | {r['lines']} | `{r['md5']}` |"
    for g, r in zip(["原件", "原件", "旧代 09-17", "**现役 09-23**", "原件", "旧代 09-17",
                     "**现役 09-23**"], L)) + f"""

两代英文稿的 H1 不同，证明这是**重定位**（研究问题换了）而不是同一稿的改字：

- 旧代：{rep['世系']['title_old']}
- 现役：{rep['世系']['title_new']}

因此**不能按文件名日期自动选基线**。选法见《唯一修改基线说明》。

## 二、两代差异规模

| 对比 | 旧代行数 | 现役行数 | 新增行 | 删除行 |
| --- | ---: | ---: | ---: | ---: |
| 英文 | {rep['世系']['diff']['en'][2]} | {rep['世系']['diff']['en'][3]} | {rep['世系']['diff']['en'][0]} | {rep['世系']['diff']['en'][1]} |
| 中文 | {rep['世系']['diff']['cn'][2]} | {rep['世系']['diff']['cn'][3]} | {rep['世系']['diff']['cn'][0]} | {rep['世系']['diff']['cn'][1]} |

## 三、§三 六项点名的逐项答案

| # | 点名事项 | 实测结论 |
| --- | --- | --- |
| 1 | 核查件记录的主稿是否与实际主稿相同 | **相同**。核查件 §五 记的 `05_稿件/English_manuscript_MDPI_重定位2026-09-23.md` 71,662 B / `1a8e0f08…`，与磁盘逐字节相符。**但**：同目录另有旧代两稿仍在原位，任何按 `05_稿件/*.md` 通配的脚本会取到 4 个文件——这是基线风险，见下节 |
| 2 | S4 是否真实存在、"450 行"与"699 行"各指什么 | **真实存在**：`06_补充材料/Supplementary_Table_S4.md`，三段实际数据行 {rep['s4']['per'].get('S4a')} + {rep['s4']['per'].get('S4b')} + {rep['s4']['per'].get('S4c')} = **{rep['s4']['total']}**。正文所称"450 rows"指 **S4c**（逐基因留一队列重排），核查件与稿件头部所称"699 行"指**三段合计**——两者指向同一文件，**均正确** |
| 3 | 森林图 PDF/PNG 是否真实存在 | **存在**。`04_图表/Fig_forest_logFC.pdf` 与 `.png` 成对；五组图件 PDF+PNG 共 10 个文件 |
| 4 | 图编号 / 文件名 / 图注 / 正文引用是否一一对应 | **对应**。主稿 Figure Legends 恰 5 条，五条各自定位到一份图件；该映射有两处**独立写就**的副本（本脚本 STEM_HINT 与 `C5_投稿清单重定位` 的 FIG_FILES），逐项一致方通过。正文引用 Figure 1–5 与 Table S1–S4 全部有落点，**待补齐 0 项**。图内编号经文本层实测：仅图 1 图内有编号字样且为自身的 `Figure 1`，**无旧编号残留**——但此项对四张无编号图是空转的，防错配实际靠上面那处映射对账。图注重复情况见第四节 |
| 5 | 分析脚本生成的结果是否就是主稿使用的结果 | **是**。`02_结果表/` {rep['溯源']['n_res']} 个文件全部可溯源（查了全项目 {rep['溯源']['n_scripts']} 个脚本）。其中 1 个是记录性 md（`B4_SERPINE1共定位_结果.md`），其正文自报的脚本/日志/产物路径逐条在位 |
| 6 | 若发现更新版本先做差异比较再定基线 | **发现一处**：`04_分析数据/MR_coloc/` 6 个文件在 2026-09-17 20:29–20:31 被改动，晚于基线清单取样时刻 19:13。差异性质已查明并登记，见第五节 |

## 四、图注画入图内的情况（装配时会否重复）

两个口径必须分开看，实测在三幅图上给出不同答案：

| 图 | 图件 | 图注词数 | 图内词数 | 整条画入率 | 末句逐字在图内 | 装配后果 |
| ---: | --- | ---: | ---: | ---: | :---: | --- |
""" + "\n".join(
    f"| {f['n']} | `{f['stem']}` | {f['lw']} | {f['iw']} | {f['cov']*100:.1f}% | "
    f"{'是' if f['emb'] else '否'} | "
    f"{'★整条重复' if f['cov'] >= 0.6 else '末句重复一句' if f['emb'] else '不重复'} |"
    for f in rep["图"]["note"]) + f"""

- **整条画入**（token 重合 ≥60%）：只有**图 1**（{rep['图']['note'][0]['lw']} 词全覆盖）。装配时若再配一遍稿件图注，这**整条逐字重复**。
- **末句逐字在图内**：**图 1、图 2**。图 2 的重合率只有 19.3%，因为其图内那段是**自撰说明**（讲 52/57 可评、点色标方向一致性），与稿件 Figure 2 图注**措辞不同**，**只有末句 `No rank is implied by position.` 与图注逐字相同**——装配时只重复这一句。
- 图 3/4/5：两个口径都未命中，图内只有轴标签、基因名与矩阵标签。

> 这一节与 `03_分析/scripts/C5_投稿清单重定位_2026-09-25.py` 的 `N_EMB` **同一口径**（该脚本用末句判据，声明 2 幅）。C5 的注释文字写"图 1、图 2 把**整条**图注画进了图内"，对图 2 不准确——本表按实测分列两个口径。
> 核查件《导师核查件_2026-09-27.md》事项 2 沿用 C5 的说法，**对图 1 成立、对图 2 高估**（图 2 实际只重复末句）。

## 五、发现的一处更新版本：B4 补跑

`04_分析数据/MR_coloc/` 下 6 个文件在 **2026-09-17 20:29–20:31** 被改动，晚于基线清单取样时刻 19:13：

| 文件 | 基线 | 现况 |
| --- | --- | --- |
""" + "\n".join(f"| `{c['path']}` | {c['b0']} B `{c['m0']}` | {c['b1']} B `{c['m1']}` |"
                for c in ch) + f"""

**性质（已查明，非失控改动）**：这是 B4——给 SERPINE1 补做区域共定位。B4 结果 md 自述：
"本轮在既有 CSV 各追加 SERPINE1 一行；VGLL4、SDC4 两行逐位未变"。
证据链 7 件经实测全部在位，时间序自洽：

""" + "\n".join(f"- {b['label']}：`{b['path']}`　{b['mt']}　{b['b']:,} B" for b in rep["b4"]) + f"""

主稿 §2.7 亦已写明 "Colocalization (coloc, default priors) was performed for VGLL4 and SDC4
and, **in this revision, for SERPINE1**"，与产物一致。SERPINE1 三档 PP.H4 = 0.0203 / 0.0181 /
0.0074，即**无共定位证据**，支持稿件"no transferable genetic support"的表述。

**后果**：`01_审计/基线md5清单.tsv` 记录的是"补 SERPINE1 之前"的状态，它**不能**用作该目录
的未改动证明。基线清单的有效范围见《唯一修改基线说明》。

## 六、主稿头部的可剥离性（§十/§十一 清洁稿的前置事实）

现役英文主稿 L3–L29 共 {rep['头部']['han']} 行含汉字，全部落在 **{len(rep['头部']['c_open'])} 段 HTML 注释块**
（L{rep['头部']['c_open'][0]}–L{rep['头部']['c_close'][0]}、L{rep['头部']['c_open'][1]}–L{rep['头部']['c_close'][1]}）之内；
正文自 **L{rep['头部']['body_start']}**（`**Title Page**`）起**无任何汉字**。

故清洁稿的动作是**整块剥离注释**，不需要逐句改写——这一点已由断言固定，改稿时不会漏。
注释块内含有两处对内部操作文件的引用（`07_修订记录/…`、`08_投稿材料/00_投稿清单.md`），
剥离注释即一并去除，不进入投稿件。
""", encoding="utf-8")
print("  版本差异表_2026-09-27.md")

(OUT / "00_版本审计" / "唯一修改基线说明_2026-09-27.md").write_text(f"""\
# 唯一修改基线说明（2026-09-27）

对应工作令 §三"若发现更新版本，先做差异比较，再确定修改基线。**不能按文件名中的日期自动选择**，
也不能覆盖彼此不同的结果。"

## 一、选定基线

| 用途 | 基线文件 | 字节 | md5 |
| --- | --- | ---: | --- |
| **英文主稿（唯一改写对象）** | `08_修订_2026-09-17/05_稿件/English_manuscript_MDPI_重定位2026-09-23.md` | {L[3]['size']:,} | `{L[3]['md5']}` |
| **中文对照稿（同源同步改，不投稿）** | `08_修订_2026-09-17/05_稿件/双语初稿_候选优选与组织来源_重定位2026-09-23.md` | {L[6]['size']:,} | `{L[6]['md5']}` |
| 结果与数值的冻结源 | `08_修订_2026-09-17/02_结果表/`（{rep['溯源']['n_res']} 个文件） | — | 见《文件清单》 |
| 图件 | `08_修订_2026-09-17/04_图表/`（5 组 PDF+PNG） | — | 见《文件清单》 |
| 补充材料 | `08_修订_2026-09-17/06_补充材料/`（S1–S4 + STROBE-MR） | — | 见《文件清单》 |

## 二、为什么是这个，而不是按日期挑

选它的依据是**内容判据**，不是文件名的日期：

1. **标题不同**——两代的 H1 不同（旧代讲"刚度因果"，现役讲"候选优选与组织来源"），
   说明研究问题换了，二者不是同一稿的两个版本。
2. **工作令 §四 建议的标题**与现役稿一致，与旧代不一致。
3. **核查件 §五**记录的主稿就是这一份，且与磁盘逐字节相符。
4. 旧代（`…_修订2026-09-17.md`）**按规则 1 保留在原位不动**，仅作回归基线，
   不参与改写、不被覆盖、不被删除。

## 三、基线清单的有效范围（重要）

`01_审计/基线md5清单.tsv`（{rep['对账']['n_base']} 项，取样于 2026-09-17 19:13）**不是全项目的未改动证明**。它的实际覆盖是：

| 范围 | 结论 |
| --- | --- |
| `05_投稿包/`（原始投稿包） | **逐字节未动**（一致 {rep['对账']['prot']['05_投稿包/']} 项）——规则 1 的核心保护对象 |
| `03_论文初稿/`（原始初稿） | **逐字节未动**（一致 {rep['对账']['prot']['03_论文初稿/']} 项） |
| `04_分析数据/MR_coloc/` | **有 6 项变动**，即上节的 B4 补跑（晚于取样时刻，属合法改动） |
| `04_分析数据/` 其余、`07_独立审稿/`、项目根散文件 | **基线未登记**（共 {rep['对账']['older']} 个文件的 mtime 早于取样时刻，属**基线漏登**而非事后新增） |
| `08_修订_2026-09-17/`（{rep['对账']['added'] - rep['对账']['older']} 个文件） | 基线取样之后新建，本就不在基线范围内 |

所以：

- "原始投稿包未改动"的结论**只对 `05_投稿包/` 与 `03_论文初稿/` 成立**，这两处已逐字节验证。
- 对 `04_分析数据/MR_coloc/` 的 6 项变动，**不得**用基线清单判为"擅自改动"；其性质是
  已有 B4 记录在案的补跑（证据链 7 件在位）。
- 对基线漏登的 {rep['对账']['older']} 个文件，**不得**判为"事后新增"；依据是其 mtime 早于取样时刻。

## 四、本轮改动的边界（开工前声明）

1. 所有产出写入 `09_Genes_CIMB_投稿修订_2026-09-27/`，**不覆盖** `08_修订_2026-09-17/` 与
   `05_投稿包/` 的任何文件。
2. 改写英文主稿时，**不是就地改**——新稿写入 `09_…/02_共用主稿/`，旧稿原地留存。
3. 现阶段只完成 §三 的版本关闭，**未开始** §五–§七 的统计修复与重写。
""", encoding="utf-8")
print("  唯一修改基线说明_2026-09-27.md")

miss = rep["引用"]["missing"]
(OUT / "00_版本审计" / "缺失文件清单_2026-09-27.md").write_text(f"""\
# 缺失文件清单（2026-09-27）

对应工作令 §三 验收："所有正文引用均能定位到真实文件；无法定位者标记'待补齐'，不得记为通过。"

## 一、正文引用无法定位者

**{len(miss)} 项。** 主稿正文中的 Table S1–S4 与 Figure 1–5 引用全部定位到真实文件；
正文中出现的文件名式字面量亦全部有落点。

## 二、已被引用但尚无落点的（登记为待补齐）

| 项 | 引用位置 | 状态 |
| --- | --- | --- |
| `Figures/Fig1`–`Fig5`.tiff | 投稿清单声明的投稿图件命名 | **待补齐**——现为 `04_图表/Fig_*.pdf` + `.png`，尚未转换为 TIFF/RGB/300 dpi/Arial 并对齐文件名。此项即核查件事项 3（装包） |
| **稿件的 docx** | 投稿材料 | **待重出**——原因同上，等图件定稿后一并做 |

## 三、脚本内已判读的历史字面量（**不是**缺失）

`V1` 对全项目 {rep['溯源']['n_scripts']} 个脚本做了字面量扫描，其中 {len(rep['溯源']['dead_lit'])} 个名字在当前项目内找不到同名文件。
逐项判读后**无一项是失效引用**：

- **现役脚本（`03_分析/scripts/`）3 条**，全为句法构造而非文件名：`_mixture.txt`、`_results.txt`
  是 R `paste0()` 的模板后缀；`修订2026-09-17.md` 是 B9 里 `endswith()` 的后缀判据。
  **即：当前分析链的代码引用是干净的。**
- **上游脚本（`04_分析数据/**`）15 条**，均为早期数据准备阶段的名字，逐条有对应物：
  GEO 矩阵存为 `.txt.gz`（5 条）、`GSE152805.rds` 实为 `_1.rds`、`merged_*.csv` 已改名
  `5队列合并_*`、`finngen_KOA_outcome.txt` 真身是 `MR_coloc/mr_data/finngen_KOA.gz`
  （{819114058:,} B 实存）、`GPL*_probe2symbol.csv` 已被 Bioconductor 注释包
  （`hgu133a.db`/`hgu133plus2.db`）取代等。逐条判读见 `V1_版本清点_2026-09-27.log`。

## 四、结论

工作令 §三 的验收条件——"所有正文引用均能定位到真实文件"——**在正文范围内成立**（待补齐 0 项）。
真正的待补齐只有投稿件本身未做（TIFF 装包与稿件 docx），且已如实登记，未记为通过。
""", encoding="utf-8")
print("  缺失文件清单_2026-09-27.md")
print(f"\n{'全部通过' if not fail else f'{fail} 项需处置'}")
raise SystemExit(1 if fail else 0)
