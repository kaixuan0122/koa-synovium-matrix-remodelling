#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
C3 — 两封投稿信按 2026-09-23 重定位改写（第 5 步成员 #4/#5）

做法：整段替换，不做子串微调。两封信的正文段落各为 markdown 的一行，
故按「行首标记」唯一匹配后整行替换；签名块、称呼、抬头与两个机制
（MDPI 五本刊共用 / JIR Article Collection）保持源件不变。

为什么整段换而不是改数字：
  旧信有四处内容在新稿中**已无载体**，逐个改数字会留下自相矛盾的句子——
    · 「PP.H4=0.742」——新稿 PP.H4 出现 0 次（共定位结果未报告）
    · 「positive control / 315th–1586th by |log2FC|」——新稿均无
    · 「12.8% vs <1.6% in every other cell type」——新稿为 1.6%–6.1%
    · 旧「What we found」整段以 LOXL1 为主语陈述，与重定位后的
      「判据依赖性」主线相冲突
  故四段（problem / what we did / what we found / two findings）+ Fit 段重写。

输出：08_投稿材料/Cover_Letter_MDPI.md（活信，交付目录根部）
      08_投稿材料/_archive_JIR作废_2026-09-27/Cover_Letter.md + README.md（归档，见下）

JIR 母版信的处置（2026-09-27 裁定）：
  JIR 阶梯已于 2026-09-15 作废并整体归档（`05_投稿包/_archive_JIR作废_2026-09-15/`），
  但底稿 §二 行 4 仍把该信列为成员，本轮遂照做重写。两者张力由用户裁定为**归档处置**：
  该信仍逐字产出（规则 1 要求保留产物），但落点移出交付目录根部——否则打开
  `08_投稿材料/` 会看到一封指向作废刊物的信与 MDPI 活信并列，且从文件名无从分辨。
"""
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[3]
SRC_DIR = ROOT / "03_论文初稿"
DST_DIR = ROOT / "08_修订_2026-09-17" / "08_投稿材料"
JIR_ARCHIVE = DST_DIR / "_archive_JIR作废_2026-09-27"

ARCHIVE_NOTE = """# JIR 母版投稿信 — 随 JIR 路线一并归档

**本目录是归档件，不是交付件。** 交付目录根部只保留 `Cover_Letter_MDPI.md`。

## 为什么在这里

1. **路线已作废。** JIR（*Journal of Inflammation Research*）阶梯于 **2026-09-15** 作废并整体归档，
   原始投稿包见 `05_投稿包/_archive_JIR作废_2026-09-15/`（按规则 1 原样保留，未改动）。
   首选期刊为 **IJMS**（*International Journal of Molecular Sciences*）。
2. **底稿仍把它列为成员。** `07_修订记录/第5步_五处同步_清单与增量_2026-09-23.md` §二 行 4
   把该信列为「待重写」，§五·二 只裁定「JIR 稿与 JIR 投稿清单不同步」，两者存在张力。
   故 **2026-09-25** 的 C3 照底稿执行，重写了该信并把张力登记在
   `第5步_五处同步_回执_2026-09-25.md` §八「待用户裁定」。
3. **2026-09-27 裁定：归档。** 用户处置指令「归档处置」即指此项。

## 处置方式

- 信件**内容一字未少**——仍由 C3 逐字生成，与 MDPI 母版同源同校（C4 独立对账覆盖两封）。
  归档不等于删除：规则 1 要求保留产物，且日后若恢复 JIR 路线可直接取用。
- 只是**落点改变**：不再写在 `08_投稿材料/` 根部。C3 末尾有反向守卫
  （交付目录根部若再出现 `Cover_Letter.md` 即失败），确保重跑不会把它复活。
- 引用本信时请用**本目录下的路径**；`00_投稿清单.md` 与 C7 成员表已同步。
"""


TITLE = ("Criterion-Dependent Priority of Matrix-Remodelling Candidates in Knee "
         "Osteoarthritis Synovium: Cross-Cohort Replication and Donor-Resolved Localisation")
DATE = "25 September 2026"

SRC_NOTE = """<!-- Source: 03_论文初稿/{name}. Rewritten for the 2026-09-23 repositioning (step-5 sync, member #4/#5).
     Whole paragraphs were replaced, not patched number by number: the previous letter stated four things
     that this submission no longer reports anywhere. (1) It quoted a colocalization posterior
     (PP.H4=0.742); the revised manuscript reports no colocalization result. (2) It described a
     machine-learning positive control in GSE55235 with panel ranks "315th-1586th by |log2FC|"; that
     paragraph is not in the revised manuscript. (3) It gave LOXL1's fibroblast localisation as
     "12.8% vs <1.6% in every other cell type"; the revised range is 1.6%-6.1%. (4) Its "What we found"
     paragraph led on LOXL1 as the finding, which the repositioning replaced with criterion-dependence.
     Every number and every quoted phrase below was re-checked against
     08_修订_2026-09-17/05_稿件/English_manuscript_MDPI_重定位2026-09-23.md and its tables; the
     reconciliation is run independently by C4. Nothing was carried over on the strength of the
     previous letter alone. The journal-specific mechanism of each master is unchanged. -->"""

# ── 段落文本（EN；两封信共用，Fit 段按刊分写） ──────────────────────────────

P_OPEN_MDPI = ('On behalf of my co-authors, I am pleased to submit our original research article entitled '
               f'**"{TITLE}"** for consideration as an Article in *International Journal of Molecular Sciences*.')

P_OPEN_JIR = ('On behalf of my co-authors, I am pleased to submit our original research article entitled '
              f'**"{TITLE}"** for consideration for publication in *Journal of Inflammation Research*, for the '
              'Article Collection "Imaging, Omics and Artificial Intelligence in Musculoskeletal Inflammation '
              'and Joint Degeneration."')

P_PROBLEM = ('**The problem.** Joint stiffness is a defining symptom of knee osteoarthritis (KOA) and a direct '
             'readout of the mechanical state of intra-articular soft tissues. Synovial fibrosis — fibroblast '
             'activation, excessive extracellular matrix deposition, and increased collagen cross-linking — is a '
             'major pathological basis of that stiffness. Work in this field is largely **output-type**: a study '
             'nominates a set of matrix-remodelling or mechanotransduction genes, and the nomination is the '
             'result. Two questions have therefore remained open in the synovium. First, how much does a '
             'candidate\'s priority depend on the criterion used to establish it — on the meta-analytic '
             'estimator, and on whether candidates are ordered by effect size or by significance? Second, which '
             'synovial cells carry the candidates, given that bulk tissue mixes compartments and that a '
             'nominate-and-validate design cannot place a case–control difference in a cell type? Neither '
             'question has been asked of this candidate space.')

P_DID = ('**What we did.** We applied one uniform per-cohort evaluation to 57 annotation-driven matrix-remodelling '
         'and mechanotransduction candidates (13 functional categories) across six public synovial cohorts — one '
         'discovery layer (three cohorts, 74 samples: 34 OA / 40 normal) and one replication layer (three cohorts, '
         '57 samples: 30 OA / 27 normal). For every candidate we report per-cohort effect sizes with confidence '
         'intervals and heterogeneity, meta-analysis under a declared primary estimator and a random-effects '
         'counterpart, and leave-one-cohort-out re-ranking. On that basis we compared the eight cross-linking '
         'candidates of the curated set under one framework for the first time, localised the candidates in '
         'donor-resolved single-cell data from three OA synovial samples, and report an exploratory *cis*-eQTL '
         'Mendelian randomisation (MR) assessment against FinnGen R12 knee osteoarthritis under three exposures — '
         'tissue-matched OA synovium, GTEx fibroblast and macrophage. Machine learning and functional enrichment '
         'are reported in the Supplementary Materials, with their calibrations.')

P_FOUND_1 = ('**What we found.** **Priority was criterion-dependent rather than a property of the candidates.** In '
             'the discovery layer the six reported panel genes behaved as previously described (LOXL1 up in all '
             'three cohorts; SDC4, SERPINE1 and VGLL4 down in all three), but the conclusions move with the '
             'estimator: SDC4 falls from P=6.6×10⁻⁷ under the declared fixed-effect estimator to P=0.111 under '
             'random effects, and ITGA5 from 0.022 to 0.477. **In the replication layer, analysed cohort by cohort '
             'rather than in aggregate, only LOXL1 and SERPINE1 are direction-consistent in all three cohorts**, '
             'and four of the six panel genes lose significance under random effects — LOXL1\'s replication P '
             'rising from 0.0028 to 0.0601. Under leave-one-cohort-out, only LOXL1 and SERPINE1 keep both their '
             'direction and P<0.01 in every deletion; LOX loses significance when GSE89408 is removed, PLOD1 when '
             'GSE206848 is removed, PLOD2 when either is removed, and TGM2 changes sign. Comparing the eight '
             'cross-linking candidates under one framework, only LOXL1 and LOX are significant in both layers '
             'under each layer\'s primary estimator and **only PLOD1 is significant under random effects in both**, '
             'while **LOXL4 has the smallest discovery P of the family (3.0×10⁻⁴) and no replication data at '
             'all**, being evaluable in two cohorts and absent from the GPL96 platform. The family therefore '
             'yields a conditional ordering rather than a single winner: under a fixed-effect or pooled estimator '
             'with a |log2FC| criterion LOXL1 is first, under random effects it is not; PLOD1 is first by '
             'random-effects significance in both layers, and SERPINE1 is first by cohort robustness with the '
             'larger effects in both layers. We report LOXL1 as the case that exhibits criterion-dependence, not '
             'as the answer to it, and we give the criterion alongside every ranking.')

P_FOUND_2 = ('**The cellular source is bounded rather than assigned.** The three single-cell donors were 86.83%, '
             '87.54% and 43.89% fibroblasts — **a spread of 43.65 percentage points** — so per-donor values are '
             'reported throughout and between-donor fractions are of limited comparability; the data localise '
             'candidates within OA tissue only. LOXL1 has its highest expression fraction in fibroblasts (12.80%, '
             'against 1.6%–6.1% in the other four cell types) and SDC4 likewise (37.48%), whereas VGLL4 is broadly '
             'expressed (65.45% in endothelial cells) and SERPINE1 is highest in endothelial cells (69.10%) but '
             'present in the other three at 17.13%–39.60%. The dataset contains **no normal controls**, so the '
             'bulk case–control difference cannot be assigned to a cell type, and a fibroblast-state association '
             'that might have narrowed this was tested against two null calibrations and was not established '
             '(P=0.31).')

P_FOUND_3 = ('**No transferable genetic support.** Across 18 gene–exposure combinations six estimates reach '
             'nominal significance, and the two with the smallest P values — VGLL4 under the macrophage exposure '
             '(β=−0.036, P=6.5×10⁻⁵) and SERPINE1 under the tissue-matched synovial exposure (β=−0.039, P=0.002) '
             '— are both single-variant Wald ratios, for which the reported P is by construction that of the '
             'variant\'s own association with the outcome, so no across-combination correction applies to them and '
             'none is claimed. VGLL4\'s evidence differs between the three eQTL sources (smallest eQTL P being '
             '6.8×10⁻⁴ in synovium, 1.8×10⁻⁶ in fibroblast and 3.2×10⁻¹¹ in macrophage), and its only nominally '
             'significant multi-instrument estimate, in GTEx fibroblasts, rests on two instruments that disagree '
             'with each other and falls to P=0.639 under stricter LD pruning. An instrument that predicts VGLL4 '
             'expression in macrophages is not thereby a valid instrument in synovial fibroblasts.')

P_NOSUPPORT = ('**Findings that did not support our initial hypothesis, and how we report them.** We draw your '
               'attention to these explicitly, because how they are handled is a fair test of the paper\'s rigour. '
               'First, **the study\'s own leading candidate is its most fragile headline**: LOXL1 reaches P=0.0028 '
               'in the replication layer under its primary pooled estimator but P=0.0601 under the random-effects '
               'counterpart, and it leads the cross-linker family only under a fixed-effect or pooled criterion. We '
               'therefore report it as the case that exhibits criterion-dependence rather than as the answer to it, '
               'and its random-effects value is given wherever its replication significance is stated. Second, '
               '**the genetic layer delivered no transferable causal support in the tissue of interest**: the two '
               'smallest-P estimates are single-variant Wald ratios that carry no across-combination correction, '
               'and the one multi-instrument signal that is nominally significant rests on two instruments that '
               'disagree with each other and does not survive stricter pruning. We report the MR assessment as '
               'exploratory throughout rather than as causal evidence, and we attribute the causal warrant for '
               'VGLL4 to the independent experimental literature rather than to our instruments. Third, **the '
               'cell-type question is bounded rather than closed**: the single-cell dataset has no normal controls, '
               'and the fibroblast-state association was tested against two null calibrations and was not '
               'established (P=0.31), so the localisation results stand as localisation only. Machine learning '
               'likewise does not support a diagnostic reading of the panel — the authoritative internal '
               'cross-validated AUC is 0.647 against a random gene-set null median of 0.635 (P=0.438), and the '
               'external value of 0.780 has no null distribution of its own and is not used to support any claim. '
               'For the remaining genes we report MR estimates as *not interpretable* rather than as negative, and '
               'we state which control each significance claim rests on — the genome-wide Benjamini–Hochberg '
               'procedure, or the panel Bonferroni threshold α=0.05/7≈0.00714, which LOXL1 at Fisher P=0.008 does '
               'not pass.')

P_FIT_MDPI = ('**Fit with the journal.** *International Journal of Molecular Sciences* publishes molecular-level '
              'studies of disease mechanism that integrate multi-omics data with causal inference, and this '
              'manuscript sits in that space: multi-cohort transcriptomic integration with per-cohort interval '
              'estimation and explicit criterion sensitivity, a same-framework comparison of a collagen '
              'cross-linking enzyme family, donor-resolved single-cell localisation, and an exploratory *cis*-eQTL '
              'Mendelian randomisation assessment — all applied to the molecular basis of joint stiffness in '
              'osteoarthritis. The repositioning of this revision is that the deliverable is a criterion-explicit, '
              'reproducible priority list for synovial matrix remodelling rather than a single named target, and '
              'the machine-learning layer is reported with its random gene-set calibration rather than as a '
              'diagnostic claim.')

P_FIT_JIR = ('**Fit with the collection.** The work sits directly at the intersection the collection names. It is '
             'omics-driven (multi-cohort transcriptomic integration with meta-analysis under two estimators), it '
             'applies genetic causal inference (exploratory MR) to separate cause from consequence, it resolves the '
             'cellular source by single-cell RNA sequencing, and it reports machine learning with an explicit '
             'random gene-set calibration rather than as a diagnostic claim — all in the context of musculoskeletal '
             'joint degeneration and the inflammatory synovial milieu.')

P_REPORT = ('**Reporting and declarations.** The study is reported in accordance with STROBE-MR, and the completed '
            'checklist is enclosed; it records, for each item, whether this submission reports it, does not report '
            'it, or does not apply to the design. This manuscript is original, has not been published previously, '
            'and is not under consideration elsewhere. All authors have read and approved the submitted version and '
            'declare no conflicts of interest. The study used only publicly available data and required no ethics '
            'approval or informed consent. All analysis scripts are available from the corresponding author upon '
            'reasonable request, and every numerical value reported was checked against the primary analysis output '
            'files.')

PARAS = [("**The problem.**", [P_PROBLEM]),
         ("**What we did.**", [P_DID]),
         ("**What we found.**", [P_FOUND_1, P_FOUND_2, P_FOUND_3]),
         ("**Two findings that did not support", [P_NOSUPPORT])]


def replace_block(lines, marker, new_paras, tag):
    """行首标记唯一匹配后，整行替换为若干段（段间空行）。"""
    hits = [i for i, l in enumerate(lines) if l.startswith(marker)]
    assert len(hits) == 1, f"{tag}: 行首标记 {marker!r} 命中 {len(hits)} 处（期望 1）"
    block = []
    for j, p in enumerate(new_paras):
        if j:
            block.append("")
        block.append(p)
    lines[hits[0]:hits[0] + 1] = block
    return lines


def build(src_name, p_open, p_fit, fit_marker):
    src = (SRC_DIR / src_name).read_text(encoding="utf-8")
    lines = src.split("\n")

    # 1. 日期
    assert lines[0].startswith("# Cover Letter"), f"{src_name}: 首行不是标题"
    di = [i for i, l in enumerate(lines) if re.fullmatch(r"\d{1,2} \w+ 2026", l)]
    assert len(di) == 1, f"{src_name}: 日期行命中 {len(di)} 处"
    lines[di[0]] = DATE

    # 2. 抬头段（投稿信本体首句）
    replace_block(lines, "On behalf of my co-authors", [p_open], f"{src_name}/opening")

    # 3. 四个内容块（What we found 由 1 段拆成 3 段）
    for mk, paras in PARAS:
        replace_block(lines, mk, paras, f"{src_name}/{mk[:24]}")

    # 4. Fit 段（按刊不同，标记也不同）
    replace_block(lines, fit_marker, [p_fit], f"{src_name}/fit")

    # 5. 声明段
    replace_block(lines, "**Reporting and declarations.**", [P_REPORT], f"{src_name}/reporting")

    # 6. 来源注释块：插在原有 HTML 注释块之后（保留中文工作提示与两个机制）
    end = [i for i, l in enumerate(lines) if l.strip() == "============================================================================ -->"]
    assert len(end) == 1, f"{src_name}: 原注释块结束标记命中 {len(end)} 处"
    lines.insert(end[0] + 1, "")
    lines.insert(end[0] + 2, SRC_NOTE.format(name=src_name))
    return "\n".join(lines)


DST_DIR.mkdir(parents=True, exist_ok=True)
JIR_ARCHIVE.mkdir(parents=True, exist_ok=True)
assert DST_DIR.is_dir(), f"输出目录未建立：{DST_DIR}"
assert JIR_ARCHIVE.is_dir(), f"归档目录未建立：{JIR_ARCHIVE}"

out_mdpi = build("Cover_Letter_MDPI.md", P_OPEN_MDPI, P_FIT_MDPI, "**Fit with the journal.**")
out_jir = build("Cover_Letter.md", P_OPEN_JIR, P_FIT_JIR, "**Fit with the collection.**")

# ────────────────────────── 后置条件 ──────────────────────────
def must(cond, msg):
    assert cond, f"后置条件失败: {msg}"
    print(f"  post OK  {msg}")


def must2(cond, ok_msg, bad_msg):
    """消息分成功/失败两版：must() 成功时也打印消息，故成功版必须是陈述句。"""
    assert cond, f"后置条件失败: {bad_msg}"
    print(f"  post OK  {ok_msg}")


for tag, out in [("MDPI", out_mdpi), ("JIR", out_jir)]:
    src = (SRC_DIR / ("Cover_Letter_MDPI.md" if tag == "MDPI" else "Cover_Letter.md")).read_text(encoding="utf-8")
    body = re.sub(r"<!--.*?-->", "", out, flags=re.S)

    # 旧标题消失、新标题恰 1 次
    must("Cross-Cohort Expression Profile and Cellular Localization" not in out,
         f"[{tag}] 旧标题已移除")
    must(out.count(TITLE) == 1, f"[{tag}] 新标题恰 1 次")

    # 新稿已不报告的断言：信件**本体**必须全部清除。
    # 注释块里允许出现（来源注释要说明删的是什么），但每一处命中都须归属到注释区间。
    spans = [(m.start(), m.end()) for m in re.finditer(r"<!--.*?-->", out, flags=re.S)]
    for bad in ["PP.H4", "0.742", "positive control", "315th", "1586", "vs <1.6% in every other cell type",
                "colocalization to separate cause from consequence", "gene-set calibration against small-cohort optimism"]:
        n_body = body.count(bad)
        must(n_body == 0, f"[{tag}] 无载体断言已从信件本体清除: {bad!r}（本体命中 {n_body}）")
        for m in re.finditer(re.escape(bad), out):
            inside = any(a <= m.start() < b for a, b in spans)
            must2(inside,
                  f"[{tag}] {bad!r} 的全部命中均在注释块内（共 {len(re.findall(re.escape(bad), out))} 处）",
                  f"[{tag}] 命中 {bad!r} 出现在注释之外：…{out[max(0, m.start() - 60):m.start() + 40]!r}")

    # 头条数字与关键值在位
    for good in ["0.0028", "0.0601", "3.0×10⁻⁴", "43.65", "0.647", "0.635", "0.438",
                 "6.5×10⁻⁵", "P=0.639", "P=0.31", "12.80%", "0.00714", "LOXL4", "PLOD1", "SERPINE1"]:
        must(good in out, f"[{tag}] 关键值在位: {good}")

    # 日期
    must(out.count(DATE) == 1, f"[{tag}] 日期已更新")

    # 段落行守恒：源 7 行（开篇 + 4 内容段 + Fit + 声明）→ 出 9 行
    # （唯一增行：What we found 由 1 段拆成 3 段，其余一一对应，无段丢失）
    n_src = sum(1 for l in src.split("\n") if l.startswith(("On behalf of my co-authors", "**")))
    n_out = sum(1 for l in out.split("\n") if l.startswith(("On behalf of my co-authors", "**")))
    must(n_src == 7 and n_out == n_src + 2, f"[{tag}] 段落行数：源 {n_src} → 出 {n_out}（期望 7 → 9）")

    # 签名块逐字保留
    sig_src = [l for l in src.split("\n") if l.startswith(("Sincerely", "Baoqiang Dong", "Professor",
                                                            "College of Acupuncture", "Liaoning University",
                                                            "Shenyang", "Email:", "Telephone:", "*On behalf of all authors:"))]
    sig_out = [l for l in out.split("\n") if l.startswith(("Sincerely", "Baoqiang Dong", "Professor",
                                                           "College of Acupuncture", "Liaoning University",
                                                           "Shenyang", "Email:", "Telephone:", "*On behalf of all authors:"))]
    must(sig_src == sig_out and len(sig_src) == 9, f"[{tag}] 签名块逐字保留（{len(sig_src)} 行）")

    # 两个机制保留
    if tag == "MDPI":
        must("*International Journal of Molecular Sciences*" in out and "Dear Editors," in out,
             "[MDPI] 抬头刊名与 Dear Editors, 机制保留")
        must(out.count("五本刊（按投稿顺序）") == 1, "[MDPI] 五本刊共用机制注释保留")
    else:
        must("*Journal of Inflammation Research*" in out and "Dear Prof. Quan," in out,
             "[JIR] 抬头刊名与具名主编机制保留")
        must("no placeholder" in out or "Suggested reviewers" in out, "[JIR] 推荐审稿人提示保留")

    # 信件本体不得含汉字（注释块除外）
    zh = set(re.findall(r"[一-鿿]", body))
    must(zh == set(), f"[{tag}] 信件本体无汉字")
    must(set(re.findall(r"[一-鿿]", out)) != set(), f"[{tag}] 中文工作提示注释保留")

    # 信件不含图号/表号（源件实测 0 处；§四.2 的「换图号」在信件中无对象）
    must(len(re.findall(r"\bFigure\b|\bTable[^ ]", src)) == 0,
         f"[{tag}] 源件无图/表编号（故本轮无需换图号）")
    must(len(re.findall(r"\bFigure\b|\bTable S?\d", out)) == 0, f"[{tag}] 产物不含图/表编号")

    print(f"  —— {tag} 共 {len(out.split(chr(10)))} 行 / {len(out.encode('utf-8'))} B")

OUTS = [(DST_DIR / "Cover_Letter_MDPI.md", out_mdpi, 4000),
        (JIR_ARCHIVE / "Cover_Letter.md", out_jir, 4000),
        (JIR_ARCHIVE / "README.md", ARCHIVE_NOTE, 400)]
for p, txt, floor in OUTS:
    p.write_text(txt, encoding="utf-8")
for p, txt, floor in OUTS:
    assert p.is_file() and p.stat().st_size > floor, f"产物异常：{p}"
    assert p.read_text(encoding="utf-8") == txt, f"回读不一致：{p}"
    print(f"  写出 {p.relative_to(ROOT)}  {p.stat().st_size} B")

# 反向守卫：既然裁定归档，交付目录根部就不得再有 JIR 母版信——否则重跑会把它复活，
# 「归档」等于没做。负面处置必须自证其生效，不能只靠「我改过脚本了」。
must(not (DST_DIR / "Cover_Letter.md").exists(),
     "交付目录根部无 Cover_Letter.md（JIR 母版未复活）")
print("\nC3 完成：MDPI 活信已落盘，JIR 母版已归档并附 README")
