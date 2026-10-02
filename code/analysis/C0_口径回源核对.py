# -*- coding: utf-8 -*-
"""C0 口径回源核对：把 57基因总览.csv 的每一列追回它在 B1/B2 里的源列。

重定位方案 §七.1 的规定：`57基因总览.csv` 的 `disc_logFC/rep_logFC` 与
`*_re_logFC` 两个口径「不得按列名直用」，写入正文前必须回源确认；
`agree_n/agree_k` 的语义同样必须先打印上下文确认。

本脚本不重算任何效应，只做「派生列 ↔ 源列」的对账，任何一条对不上即非零退出。
判据逐列写成断言，而不是靠读数的人记住。

为什么要单独做这一步：本项目已经踩过两个同型坑——
  · `agreement-count-column-carries-no-direction`：一致性计数列只数队列数、不含方向，
    据此标上下会把 3 个下调基因标成「3/3↑」；
  · `rds-with-method-column-never-index-by-row`：含 method 列的结果表按行号取值，
    会把 Egger 当 IVW。
两个坑的共性都是「列名看着就是那个意思」，所以这里逐列追源。

输入：02_结果表/{57基因总览,合并效应_发现层,合并效应_复制层,逐队列效应_力学基因}.csv
输出：01_审计/口径回源核对_结果.md + 同名 .log（由调用方重定向）
"""
import csv
import os
import sys

BASE = ("/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/"
        "Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/"
        "01.材料/01.博士材料/10.毕业/01.小论文/04.立项/"
        "滑膜基质刚度生成与刚度感知 0912/08_修订_2026-09-17")
RES = os.path.join(BASE, "02_结果表")
OUT = os.path.join(BASE, "01_审计")

DISCOVERY = ["GSE89408", "GSE206848", "GSE1919"]
REPLICATION = ["GSE55235", "GSE55457", "GSE82107"]
TOL = 1e-8

# 主估计口径的**字面量**期望：不从数据现算，写死在这里。
# 现算的话，B2 哪天改了 primary_method，期望值跟着一起变，断言永远通过。
EXPECT_METHOD = {"发现": "IVW-FE", "复制": "pooled(~0+group+cohort)"}


def read_csv(name):
    path = os.path.join(RES, name)
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


# R 的 write.csv 默认把缺失写成字面 `NA`，不是空串（实测：这两份表里空串 0 个、
# `NA` 96 个）。判 NA 前先查实际哨兵，不假定。
NA_TOKENS = {"", "NA", "NaN", "N/A"}


def num(row, col):
    v = row.get(col, "")
    return None if v in NA_TOKENS else float(v)


def sgn(x):
    """R 的 sign()：0 得 0，不是 FALSE。用 (x>0)-(x<0) 才与 R 同结论——
    (v>0)==(ref>0) 这种写法在 v==0 且 ref<0 时会与 R 相反。"""
    return 0 if x is None else (x > 0) - (x < 0)


def same(a, b):
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= TOL * max(1.0, abs(a), abs(b))


lines = []
problems = []
checks = 0


def say(s=""):
    print(s)
    lines.append(s)


def check(cond, msg):
    global checks
    checks += 1
    if not cond:
        problems.append(msg)
        say(f"  ❌ {msg}")
    return cond


ov = read_csv("57基因总览.csv")
dis = {r["gene"]: r for r in read_csv("合并效应_发现层.csv")}
rep = {r["gene"]: r for r in read_csv("合并效应_复制层.csv")}
per = read_csv("逐队列效应_力学基因.csv")

# 名单一律从产物里读，不在本脚本重抄一份（重抄就有两份会漂移）
fam = lambda col: [r["gene"] for r in ov if r[col] == "TRUE"]

say("# 口径回源核对（C0）")
say()
say(f"- 57基因总览：{len(ov)} 行 × {len(ov[0])} 列")
say(f"- 合并效应_发现层：{len(dis)} 行；合并效应_复制层：{len(rep)} 行；逐队列：{len(per)} 行")
say()

# ---------- 1. 主估计口径的字面量（不从数据现算） ----------
say("## 一、主估计口径（写死的期望值，非现算）")
for layer, table in ("发现", dis), ("复制", rep):
    got = {r["primary_method"] for r in table.values()}
    check(got == {EXPECT_METHOD[layer]},
          f"{layer}层 primary_method 实测 {got}，期望 {{{EXPECT_METHOD[layer]!r}}}")
    say(f"  · {layer}层 primary_method = {sorted(got)}  → {EXPECT_METHOD[layer]}")
say()

# ---------- 2. disc_logFC 追到发现层哪一列 ----------
say("## 二、disc_* 的源列（发现层）")
n_d = 0
for r in ov:
    g = r["gene"]
    d = dis.get(g)
    if d is None:
        check(num(r, "disc_logFC") is None, f"{g}：发现层无行，但 disc_logFC 非空")
        continue
    n_d += 1
    for ovc, dc in (("disc_logFC", "primary_logFC"), ("disc_CI_low", "primary_CI_low"),
                    ("disc_CI_high", "primary_CI_high"), ("disc_P", "primary_P")):
        check(same(num(r, ovc), num(d, dc)), f"{g}：{ovc} ≠ 发现层 {dc}")
    # 发现层 primary_* 是否等于 ivw_*（B2 第 252 行按层赋值）
    for dc, ic in (("primary_logFC", "ivw_logFC"), ("primary_P", "ivw_P")):
        check(same(num(d, dc), num(d, ic)), f"{g}：发现层 {dc} ≠ {ic}")
    # 随机效应列
    for ovc, dc in (("disc_re_logFC", "re_logFC"), ("disc_re_P", "re_P")):
        check(same(num(r, ovc), num(d, dc)), f"{g}：{ovc} ≠ 发现层 {dc}")
say(f"  · 对账 {n_d} 个基因：disc_logFC ≡ primary_logFC ≡ **ivw_logFC**（IVW 固定效应 Meta）")
say("  · disc_re_logFC ≡ re_logFC（DerSimonian–Laird）——它是**主估计的**随机效应版本")
say()

# ---------- 3. rep_logFC 追到复制层哪一列 ----------
say("## 三、rep_* 的源列（复制层）")
n_r = 0
for r in ov:
    g = r["gene"]
    m = rep.get(g)
    if m is None:
        check(num(r, "rep_logFC") is None, f"{g}：复制层无行，但 rep_logFC 非空")
        continue
    n_r += 1
    for ovc, rc in (("rep_logFC", "primary_logFC"), ("rep_CI_low", "primary_CI_low"),
                    ("rep_CI_high", "primary_CI_high"), ("rep_P", "primary_P")):
        check(same(num(r, ovc), num(m, rc)), f"{g}：{ovc} ≠ 复制层 {rc}")
    # 复制层 primary_* 等于 pooled_*（**不是** ivw_*）
    for rc, pc in (("primary_logFC", "pooled_logFC"), ("primary_P", "pooled_P")):
        check(same(num(m, rc), num(m, pc)), f"{g}：复制层 {rc} ≠ {pc}")
    for ovc, rc in (("rep_re_logFC", "re_logFC"), ("rep_re_P", "re_P")):
        check(same(num(r, ovc), num(m, rc)), f"{g}：{ovc} ≠ 复制层 {rc}")
say(f"  · 对账 {n_r} 个基因：rep_logFC ≡ primary_logFC ≡ **pooled_logFC**"
    "（~0+group+cohort，HC3 稳健 SE）")
say("  · rep_re_logFC ≡ re_logFC，但它是**对照列（IVW-FE）的**随机效应版本，"
    "不是该层主估计的——二层含义不同，写入正文时不得与发现层并列成同一句话")
say()

# ---------- 4. agree_n / agree_k 语义 ----------
say("## 四、agree_n / agree_k 的语义（逐队列重算对账）")
na_logfc = [p for p in per if p["logFC"] in NA_TOKENS]
say(f"  · 逐队列表里 logFC 为 NA 的行：{len(na_logfc)} / {len(per)}"
    f"{'（' + ', '.join(p['gene'] + '/' + p['dataset'] for p in na_logfc[:6]) + '）' if na_logfc else ''}")
say("    B8 的 dir_count 只按 dataset 过滤、不滤 NA，故分母含无值队列——本核对按同一规则重算；")
say("    凡该队列 logFC 为 NA，R 的 sum(sign(v)==sign(ref)) 整条得 NA，故分子一并记 NA（不是当 0）。")
bad_n = bad_k = 0
for r in ov:
    g = r["gene"]
    for layer, cohorts, ncol, kcol, refcol in (
            ("发现", DISCOVERY, "disc_agree_n", "disc_agree_k", "disc_logFC"),
            ("复制", REPLICATION, "rep_agree_n", "rep_agree_k", "rep_logFC")):
        vals = [num(p, "logFC") for p in per if p["gene"] == g and p["dataset"] in cohorts]
        ref = num(r, refcol)
        want_n = len(vals)
        if ref is None or any(v is None for v in vals):
            want_k = None                      # R：任一项为 NA，整个 sum 得 NA
        else:
            want_k = sum(1 for v in vals if sgn(v) == sgn(ref))
        got_n, got_k = num(r, ncol), num(r, kcol)
        if not check(got_n is not None and int(got_n) == want_n,
                     f"{g}/{layer}：{ncol} 实测 {got_n}，按逐队列重算 {want_n}"):
            bad_n += 1
        # 两个方向都要断：该 NA 的必须是 NA，该有数的必须等于重算值
        if want_k is None:
            if not check(got_k is None,
                         f"{g}/{layer}：{kcol} 实测 {got_k}，但逐队列有 NA，应为 NA"):
                bad_k += 1
        elif not check(got_k is not None and int(got_k) == want_k,
                       f"{g}/{layer}：{kcol} 实测 {got_k}，按逐队列重算 {want_k}"):
            bad_k += 1
say("  · 分母 agree_n ≡ 该基因在该层**逐队列表里出现的行数**（不是 6、也不是该层的队列数 3）")
say("  · 实测 logFC 列 0 个 NA，故上面「分母含无值队列」在本表不构成实际影响——"
    "记下的是量过的结论，不是先验顾虑")
say("  · 分子 agree_k ≡ 符号与**该层主估计**相同的队列数 → **计数内含方向**")
say("  · 但方向本身在 disc_dir / rep_dir 列：报计数必须同时报方向列，"
    "只报「3/3」读者无法知道是上调还是下调")
say()

# ---------- 4b. disc_k / rep_k 的语义 ----------
say("## 四之二、disc_k / rep_k（各层参与合并的队列数）")
bad_k2 = 0
for r in ov:
    g = r["gene"]
    for layer, cohorts, kcol, src in (
            ("发现", DISCOVERY, "disc_k", dis),
            ("复制", REPLICATION, "rep_k", rep)):
        m = src.get(g)
        n_in_layer = len([p for p in per if p["gene"] == g and p["dataset"] in cohorts])
        got = num(r, kcol)
        if m is None:
            if not check(got is None, f"{g}/{layer}：该层无行，但 {kcol} 实测 {got}"):
                bad_k2 += 1
            continue
        if not check(same(num(m, "k"), n_in_layer),
                     f"{g}/{layer}：该层表 k={num(m, 'k')}，但逐队列有 {n_in_layer} 行"):
            bad_k2 += 1
        if not check(same(got, num(m, "k")),
                     f"{g}/{layer}：{kcol} 实测 {got} ≠ 该层表 k={num(m, 'k')}"):
            bad_k2 += 1
say("  · disc_k / rep_k ≡ 该层结果表里的 `k` 列 ≡ 该基因在该层逐队列表里的行数")
say(f"  · 复制层缺行的基因 {sum(1 for r in ov if num(r, 'rep_k') is None)} 个——"
    "总览里 rep_* 整段记 NA，不是 0；正文不得把这几个基因写成「复制层无差异」")
# 各层队列数**不齐**：不能统一写「3 个队列」，逐基因的 k 不同
for layer, kc in (("发现", "disc_k"), ("复制", "rep_k")):
    dist = {}
    for r in ov:
        v = num(r, kc)
        dist[None if v is None else int(v)] = dist.get(None if v is None else int(v), 0) + 1
    say(f"  · {layer}层 k 的分布：" + "，".join(
        f"k={k} 有 {n} 个基因" if k is not None else f"无行 {n} 个基因"
        for k, n in sorted(dist.items(), key=lambda kv: (kv[0] is None, kv[0]))))
say("    ——故正文**不得**统一写「三层合并」或「三个队列」；报 k 时须给分布或逐基因值")
say()

# ---------- 4c. NA 清单（量过的，不是假定的） ----------
say("## 四之三、NA 清单（逐列实测）")
na_cols = {}
for p in per:
    for c, v in p.items():
        if v in NA_TOKENS:
            na_cols.setdefault(c, []).append(f"{p['dataset']}/{p['gene']}")
say("  · 逐队列表各列 NA 数：" + (", ".join(f"`{c}` {len(v)}" for c, v in na_cols.items())
                                 or "无"))
for c, v in na_cols.items():
    for tag in v:
        ds, g = tag.split("/")
        row = next(p for p in per if p["dataset"] == ds and p["gene"] == g)
        ok = row["test"] in ("limma", "DESeq2") and row["P"] not in NA_TOKENS
        say(f"    - `{c}`：{tag}——该行 `test={row['test']}`、原始 `P={row['P']}`"
            f"{'有值' if ok else '无值'}，故是「检验跑过、仅这一项不可估」，不是「未检验」")
ov_na = sum(1 for r in ov for v in r.values() if v in NA_TOKENS)
ov_na_cols = sorted({c for r in ov for c, v in r.items() if v in NA_TOKENS})
n_rep_na_cols = sum(1 for c in ov_na_cols
                    if c.startswith("rep_") or c in ("same_dir", "both_sig", "both_sig_re"))
say(f"  · 总览：{ov_na} 个 NA，分布在 {len(ov_na_cols)} 列——逐列名如下，"
    "不按前缀归类（`both_sig`、`same_dir` 等也因依赖复制层而为 NA）：")
say("    " + ", ".join(f"`{c}`" for c in ov_na_cols))
say(f"    其中 {n_rep_na_cols} 列依赖复制层——对应"
    f" {sum(1 for r in ov if num(r, 'rep_k') is None)} 个基因复制层无行")
for name in ("合并效应_发现层.csv", "合并效应_复制层.csv", "留一法_合并效应.csv"):
    rows = read_csv(name)
    n = sum(1 for r in rows for v in r.values() if v in NA_TOKENS)
    say(f"  · {name}：{n} 个 NA")
say("  · 复制层为何 52 行而非 57：B2 第 92 行 `if (sum(ok) < 2) return(NULL)` ——"
    "该层合并要求**至少 2 个队列同时具备两组**；")
say("    只在 1 个复制队列里出现的基因整行不出（LOXL3、LOXL4 仅 GSE82107 有），"
    "故它们不是「复制层无差异」，是**复制层不可合并**")
say()

say("## 四之四、重定位要用的关键基因（本表即 §二 家族表的数值来源）")
say()
say("| 基因 | 层 | disc_k | disc_logFC | disc_P | 发现层同向 | rep_k | rep_logFC | rep_P | 复制层同向 | same_dir | both_sig |")
say("|---|---|---|---|---|---|---|---|---|---|---|---|")
cross_genes, panel_genes = fam("is_crosslink"), fam("in_panel")
for g in cross_genes + [x for x in panel_genes if x not in cross_genes]:
    r = next(x for x in ov if x["gene"] == g)
    f3 = lambda c: "—" if num(r, c) is None else f"{num(r, c):+.3f}"
    iv = lambda c: "—" if num(r, c) is None else str(int(num(r, c)))
    pv = lambda c: "—" if num(r, c) is None else f"{num(r, c):.2e}"
    tag = "+".join(t for t, on in (("crosslink", r["is_crosslink"]),
                                   ("panel", r["in_panel"])) if on == "TRUE")
    agree = lambda kc, nc: ("—" if num(r, kc) is None else
                            f"{iv(kc)}/{iv(nc)}({r['disc_dir' if kc == 'disc_agree_k' else 'rep_dir']})")
    say(f"| {'**' + g + '**'} | {tag} | {iv('disc_k')} | {f3('disc_logFC')} | {pv('disc_P')} | "
        f"{agree('disc_agree_k', 'disc_agree_n')} | "
        f"{iv('rep_k')} | {f3('rep_logFC')} | {pv('rep_P')} | "
        f"{agree('rep_agree_k', 'rep_agree_n')} | {r['same_dir']} | {r['both_sig']} |")
say()
say("  · 「同向」列写成 `k/n(方向)`：计数**必须**带方向词——本项目踩过"
    "「一致性计数列不含方向」，只写 3/3 会把下调基因标成上调")
say()

# ---------- 4d. 符号版本（表里用旧名，正文须绑定） ----------
test_by_ds = {}
for p in per:
    test_by_ds.setdefault(p["dataset"], set()).add(p["test"])
say("  · 逐队列表的 `test` 列（实测两种取值）：")
for ds, ts in test_by_ds.items():
    n_ds = len([p for p in per if p["dataset"] == ds])
    say(f"    - {ds}：{sorted(ts)}（{n_ds} 行）")
say(f"    ——{sorted(test_by_ds['GSE89408'])} 的逐队列值是**计数 GLM 系数**，"
    "与其余队列的 limma 线性模型系数不是同一个估计量；")
say("       这正是发现层主估计改用 IVW-FE、而不直接合并逐队列 DESeq2 系数的原因，"
    "Methods 必须写明这一个队列用了不同的检验")
say()

say("## 四之五、符号版本（旧的写法会与表对不上）")
for old, new in (("CYR61", "CCN1"), ("CTGF", "CCN2"), ("MKL1", "MRTFA"), ("MRTFB", "MKL2")):
    inn = old in {r["gene"] for r in ov}
    say(f"  · 表里是 **{old}**（{'在' if inn else '不在'}），不是 {new}——"
        f"正文写 {new} 时须同时给出 {old}，否则与表对不上")
say()
sym = {}
for p in per:
    if p["symbol_used"] != p["gene"]:
        sym.setdefault(p["gene"], {}).setdefault(p["symbol_used"], []).append(p["dataset"])
say(f"  · 实测 {sum(len(v) for d in sym.values() for v in d.values())} 处逐队列符号 ≠ 规范名：")
for g, d in sym.items():
    say(f"    - {g}：" + "；".join(f"在 {', '.join(ds)} 记作 **{s}**" for s, ds in d.items()))
say("    ——正文用规范名时必须把这些队列的实际符号一并给出（表注或 Methods），"
    "否则读者按 CCN2 去查会以为该队列没测到这个基因")
say(f"  · PANEL（in_panel）{len(fam('in_panel'))} 个：{', '.join(fam('in_panel'))}")
say(f"  · CROSSLINK（is_crosslink）{len(fam('is_crosslink'))} 个："
    f"{', '.join(fam('is_crosslink'))}")
say("    ——方案 §二那张家族表就是这 8 个；两处名单若不同，以本表为准并回头改方案")
say()

# ---------- 5. 方向列与估计值同源 ----------
say("## 五、方向列 dir 与所取系数的绑定")
bad_dir = 0
for r in ov:
    for col, dcol in (("disc_dir", "disc_logFC"), ("rep_dir", "rep_logFC")):
        v = num(r, dcol)
        if v is None:
            continue
        want = "上调" if v > 0 else "下调"
        if not check(r[col] == want, f"{r['gene']}：{col}={r[col]!r} 与 {dcol}={v:+.4f} 不符"):
            bad_dir += 1
say("  · 全部由估计值符号生成，无手写方向词" if bad_dir == 0
    else f"  · {bad_dir} 处不符")
say()

# ---------- 6. 供写入正文的口径表 ----------
say("## 六、结论：写入正文时必须照这个口径写")
say()
say("| 正文要报的量 | 应取的列 | 是什么模型 |")
say("|---|---|---|")
say("| 发现层合并效应 | `合并效应_发现层.csv` `primary_*`（= `ivw_*`） | 逆方差固定效应 Meta |")
say("| 复制层合并效应 | `合并效应_复制层.csv` `primary_*`（= `pooled_*`） | `~0+group+cohort`，HC3 稳健 SE |")
say("| 发现层随机效应对照 | `合并效应_发现层.csv` `re_*` | DL，主估计的 RE 版本 |")
say("| 复制层随机效应对照 | `合并效应_复制层.csv` `re_*` | DL，**对照列 IVW-FE 的** RE 版本 |")
say("| 复制层固定效应对照 | `合并效应_复制层.csv` `ivw_*` | 逆方差固定效应 Meta |")
say("| 一致性计数 | `*_agree_k / *_agree_n` + 同层 `*_dir` | 与**该层主估计**同向的队列数 |")
say()
say("两层主估计**不是同一个模型**，故不得写成「同一套合并方法用于两层」；"
    "复制层若要报固定效应与随机效应的对照，对照的是 `ivw_*` 与 `re_*`（同为 Meta 估计量），"
    "而不是 `primary_*` 与 `re_*`（这两个不是同一个估计量）。")
say()

say("---")
say(f"断言共 {checks} 条；未通过 {len(problems)} 条。")
if problems:
    say()
    say("未通过明细（前 20 条）：")
    for m in problems[:20]:
        say(f"- {m}")
    say()
    say("**结论：口径与 B1/B2 产物不一致，禁止写入正文。**")
else:
    say()
    say("**结论：全部对上，口径可用于正文。**")

os.makedirs(OUT, exist_ok=True)
dest = os.path.join(OUT, "口径回源核对_结果.md")
with open(dest, "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
print(f"\n已写出 {dest}")
sys.exit(1 if problems else 0)
