#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
C1 — STROBE-MR checklist 重映射（第 5 步成员 #9）

输入: 03_论文初稿/STROBE-MR_checklist.md                        （源，未改动）
参照: 08_修订_2026-09-17/05_稿件/English_manuscript_MDPI_重定位2026-09-23.md
输出: 08_修订_2026-09-17/06_补充材料/STROBE-MR_checklist.md

做法：逐条锚点替换，每处断言在源中恰好出现 1 次；行级替换按行首前缀唯一匹配。
不按偏移推编号（旧 §3.x 七节 → 新 §3.x 四节），每一处 § 引用均回新稿核实后重写。
缺失项按规则 4 记为 Partial / No 并明写「本稿未报告」，不在文字中模拟。

后置条件在文件末尾统一断言。
"""
import re
import pathlib
import hashlib

ROOT = pathlib.Path(__file__).resolve().parents[3]
SRC = ROOT / "03_论文初稿" / "STROBE-MR_checklist.md"
DST = ROOT / "08_修订_2026-09-17" / "06_补充材料" / "STROBE-MR_checklist.md"

src = SRC.read_text(encoding="utf-8")
out = src
log = []


def sub(old, new, tag):
    """唯一性断言的子串替换。"""
    global out
    n = out.count(old)
    assert n == 1, f"[{tag}] 锚点出现 {n} 次（期望 1 次）: {old[:80]!r}"
    out = out.replace(old, new)
    log.append(f"  sub  {tag}: 1 处")


def line(prefix, new_line, tag):
    """按行首前缀唯一匹配的整行替换。"""
    global out
    lines = out.split("\n")
    idx = [i for i, l in enumerate(lines) if l.startswith(prefix)]
    assert len(idx) == 1, f"[{tag}] 前缀 {prefix!r} 匹配 {len(idx)} 行（期望 1 行）"
    old_line = lines[idx[0]]
    lines[idx[0]] = new_line
    out = "\n".join(lines)
    log.append(f"  line {tag}: 替换 1 行（旧长 {len(old_line)} → 新长 {len(new_line)}）")


# ── 1. 首页：标题与设计段 ────────────────────────────────────────────────
sub(
    "**Manuscript:** Cross-Cohort Expression Profile and Cellular Localization of "
    "Matrix-Remodelling and Mechanotransduction-Related Genes in Knee Osteoarthritis Synovium",
    "**Manuscript:** Criterion-Dependent Priority of Matrix-Remodelling Candidates in "
    "Knee Osteoarthritis Synovium: Cross-Cohort Replication and Donor-Resolved Localisation",
    "title",
)

line(
    "**Study design:**",
    "**Study design:** Two-sample summary-data Mendelian randomization, performed in parallel "
    "under three *cis*-eQTL exposures relevant to different degrees (eQTL Catalogue release 7): "
    "tissue-matched OA synovium (QTD000519, n=72), GTEx fibroblast (QTD000216, n=483) and "
    "Alasoo 2018 macrophage (QTD000001, n=84). Outcome: FinnGen R12 knee osteoarthritis "
    "(M13_ARTHROSIS_KNEE; 61,356 cases / 315,115 controls). Six genes were interrogated, giving "
    "18 gene–exposure combinations. The MR analyses are reported as exploratory and as a bounded "
    "appendix to the localisation result: the genetic association evidence for VGLL4 differs "
    "between the three eQTL sources used, and the two smallest *P* values in the MR set come from "
    "rows resting on a single variant — for which the reported *P* is by construction that of the "
    "variant's own association with the outcome and to which no across-combination correction "
    "applies (§2.8) — so no panel gene carries transferable genetic support in synovium, and each "
    "result is bounded by the cell type in which its instrument was derived.",
    "study-design",
)

# ── 2. 首页注释：加英文来源说明（保留原中文工作提示） ────────────────────
sub(
    "     ============================================================================ -->\n",
    "     ============================================================================ -->\n\n"
    "<!-- Source: 03_论文初稿/STROBE-MR_checklist.md. Synced for the 2026-09-23 repositioning.\n"
    "     Every section anchor and every quotation below was re-checked item by item against\n"
    "     08_修订_2026-09-17/05_稿件/English_manuscript_MDPI_重定位2026-09-23.md and its tables;\n"
    "     no anchor was remapped by offset. Items marked Partial or No record what the revised\n"
    "     manuscript no longer reports, and nothing was simulated in place of a missing result.\n"
    "     The \"Item\" column remains a paraphrase. -->\n",
    "provenance-comment",
)

# ── 3. item 1：摘要引句与题名描述 ────────────────────────────────────────
sub(
    '("Exploratory Mendelian randomization, using three cis-eQTL exposures (synovium, '
    'fibroblast, macrophage) against FinnGen R12 KOA, gave no transferable causal support …")',
    '("Exploratory Mendelian randomisation gave no transferable genetic support …")',
    "item1-abstract-quote",
)
sub(
    "it describes a cross-cohort expression and localization study and names neither a method nor a gene",
    "it describes a cross-cohort replication and donor-resolved localisation study and names "
    "neither a method nor a gene",
    "item1-title-description",
)

# ── 4. item 2：MR 理由未单独展开 ─────────────────────────────────────────
line(
    "| 2 |",
    "| 2 | Background/rationale — exposure, plausibility, why MR | Partial | §1, paragraphs 1–2 "
    "give the biological rationale (stiffness generation by the LOX/LOXL family; stiffness sensing "
    "through integrins and syndecans and Hippo-mediated response), and §1, paragraph 3 sets out the "
    "published work this study is positioned against — of the five named studies, three used "
    "synovium ([7], [11], [12]) and two of those supported their nominations by Mendelian "
    "randomisation [11,12]. The rationale **for MR as "
    "a method is not developed separately**: §1, paragraph 4 carries the genetic module as "
    "\"exploratory genetic evidence as a bounded appendix\" rather than as a stated causal hypothesis. |",
    "item2",
)

# ── 5. item 3：问句在 §1 第 4 段（原第 3 段） ────────────────────────────
sub(
    '§1, paragraph 3 states the questions asked ("… and which synovial cell type carries it")',
    '§1, paragraph 4 states the questions asked ("We therefore asked which matrix-remodelling and '
    'mechanotransduction-related candidates change reproducibly in KOA synovium, how far their '
    'priority depends on the criterion used to establish it, and which synovial cells carry them.")',
    "item3-paragraph-and-quote",
)

# ── 6. item 4a ───────────────────────────────────────────────────────────
line(
    "| 4a |",
    "| 4a | Study design and setting | Yes | §2.7; the design is also identified in the Abstract "
    "(\"Exploratory Mendelian randomisation gave no transferable genetic support …\"). The study "
    "setting and the four evidence layers are shown in Figure 1. |",
    "item4a",
)

# ── 7. item 4b：MR 无功效计算；MDE 属转录组层 ────────────────────────────
line(
    "| 4b |",
    "| 4b | Participants — eligibility, selection, **sample size / power** | Partial | §2.7 gives "
    "the sample size of each of the three eQTL exposures (n=72 / 483 / 84) and of the outcome GWAS "
    "(61,356 cases / 315,115 controls). **No power calculation is reported for the MR module, and "
    "the revised manuscript reports no minimum detectable effect for the MR estimates.** The "
    "minimum detectable effect defined in §2.3 (2.802×SE at α=0.05 and 80% power) belongs to the "
    "transcriptomic layers, where it is used to classify uninformative genes; the cohort-level "
    "values in the Table 2 footnote (genome-wide median 0.952 log2FC for GSE82107 against 0.642 and "
    "0.643; 1.38 against 1.25 on the common d scale) are properties of the replication cohorts, not "
    "of the MR estimates. |",
    "item4b",
)

# ── 8. item 4c：变体表与其脚注已不在本稿 ─────────────────────────────────
line(
    "| 4c |",
    "| 4c | Genetic variant measurement, QC and selection | Partial | §2.7: the top three "
    "*cis*-eQTLs by *P* within ±500 kb were taken as candidate instruments, with no P-value "
    "threshold applied so that the selection rule was identical across exposures; **instrument "
    "independence was established by measured linkage disequilibrium rather than by physical "
    "distance** — "
    "pairwise r² computed from the 1000 Genomes phase 3 GRCh38 phased release (n=2,548; 503 "
    "European-ancestry samples), queried by genomic position over remote tabix, then greedy pruning "
    "to r²≤0.1 by ascending eQTL *P*, with variants absent from the reference panel discarded rather "
    "than assumed independent. Measured LD exceeded r²=0.1 for six pairs and the seven affected "
    "variants were removed; the synovial exposure was unaffected (maximum pairwise r²=0.0065). §2.7 "
    "also states that the ±500 kb window is narrower than the *cis*-region as defined by the eQTL "
    "Catalogue, so every gene's extracted variants reached the window boundary and the window is a "
    "stated methodological constraint rather than a complete *cis*-region capture. **Two QC records "
    "of the previous version are not carried:** the variant-level table listing effect/other allele, "
    "EAF, eQTL β (SE) and *P*, and KOA β (SE) for each instrument, and its footnote recording every "
    "variant removed and why. The removals are now summarised in §2.7 (six pairs, seven variants) "
    "and per gene–exposure set in the \"IV set changed\" column of Table S3; per-variant alleles and "
    "EAF are not reported anywhere in this submission. |",
    "item4c",
)

# ── 9. item 4d ───────────────────────────────────────────────────────────
line(
    "| 4d |",
    "| 4d | Assessment of exposures, outcomes and other variables | Yes | §2.7; the three exposures "
    "are named with their eQTL Catalogue identifiers and sample sizes, and the disease outcome "
    "definition (FinnGen R12, M13_ARTHROSIS_KNEE) is given. §2.7 also states why release 7 is "
    "specified — the same GTEx accession was re-processed as GTEx v10 in the release 8 pre-release "
    "(n=659), so the accession alone does not identify the data. |",
    "item4d",
)

# ── 10. item 4e：知情同意原文已不同 ──────────────────────────────────────
sub(
    "Conversely, **no participant-level consent was obtained, because none was required**; the "
    "Informed Consent Statement says the same.",
    "Conversely, **no participant-level consent was obtained, because none was required**; the "
    "Informed Consent Statement reads \"Not applicable. All data analyzed were obtained from public "
    "repositories in de-identified form.\"",
    "item4e-consent",
)

# ── 11. item 5：三项假设逐字重取 ─────────────────────────────────────────
line(
    "| 5 |",
    "| 5 | **Core instrumental-variable assumptions** (relevance, independence, exclusion "
    "restriction) | Yes | §2.7, third paragraph. Quotations below are verbatim from that paragraph. "
    "*Relevance*: each instrument \"is a *cis*-eQTL for the exposure gene, quantified by the "
    "gene-level F statistic.\" *Independence*: \"instruments are assumed unrelated to confounders "
    "of the exposure–outcome relationship — **stated, not tested**, since LD pruning addresses "
    "correlation between instruments and cannot establish freedom from confounding [15].\" Residual "
    "correlation among retained instruments is stated in §2.7, second paragraph — \"six retained "
    "pairs still exceed r²=0.01, the largest being 0.063 (GTEx fibroblast VGLL4)\" — together with "
    "the statement that no correction was applied for it and that the resulting bias in the pooled "
    "standard error was \"**not quantified** rather than assumed small\". *Exclusion restriction*: "
    "\"instruments are assumed to affect KOA only through the exposure gene, which is not directly "
    "testable; the MR-Egger intercept tested directional horizontal pleiotropy where three or more "
    "instruments survived, and Cochran's Q tested heterogeneity.\" The same paragraph states that "
    "the exclusion restriction is **not directly testable**. **The two analyses reported in place of "
    "an independence check are no longer framed as such:** the decomposition into single-variant "
    "Wald ratios survives only as the §2.8 statement that a single-instrument row is not one of 18 "
    "combinations under test, and the re-pruning at r²≤0.01 is now pointed to Table S3 as a "
    "threshold comparison rather than described as a sensitivity analysis. Per-variant and per-set "
    "details, including per-SNP F, are given in Table S3. |",
    "item5",
)

# ── 12. item 6a：单位说明只在 §2.7，Table 2 脚注无此项 ───────────────────
sub(
    "| 6a | Handling of quantitative variables | Yes | §2.5 and the Table 2 footnote. β is the effect",
    "| 6a | Handling of quantitative variables | Yes | §2.7. β is the effect",
    "item6a",
)
sub(
    "so the per-allele eQTL β in Table 3 is already in units of the SD of normalised expression",
    "so the per-allele eQTL β is already in units of the SD of normalised expression",
    "item6a-stale-table3",
)

# ── 13. item 6b：协调化未描述但从分析实际执行 ────────────────────────────
line(
    "| 6b |",
    "| 6b | Handling of genetic variants and weights | Partial | §2.7 (*cis*-eQTL used as "
    "instrument; Wald ratio where one instrument survived, IVW where two or more). **Allele "
    "harmonisation is not described anywhere in this submission**: the sentence that described it in "
    "the previous version stood in the footnote of the variant-level table, which this revision does "
    "not carry. The analysis does perform it, so this is a **reporting gap in the revised "
    "manuscript, not a missing analysis** (the executing code is named in the accompanying revision "
    "record), and the description should be restored before submission. |",
    "item6b",
)

# ── 14. item 6c：Wald/IVW 阈值逐字 ───────────────────────────────────────
line(
    "| 6c |",
    "| 6c | MR estimator and related statistics | Yes | §2.7: \"Genes retaining one instrument were "
    "analysed by Wald ratio and genes retaining two or more by IVW; MR-Egger and the weighted median "
    "were available only where three or more instruments survived.\" The revised manuscript does not "
    "name one estimator as primary **within the MR module** — the \"declared primary estimator\" "
    "language of §2.3 and Table 1 refers to the transcriptomic meta-analysis, not to MR. |",
    "item6c",
)

# ── 15. item 6d：丢弃规则保留，其余三条未承载 ────────────────────────────
line(
    "| 6d |",
    "| 6d | Missing data | Partial | §2.7: variants absent from the 1000 Genomes reference panel, "
    "whose LD was therefore unverifiable, were discarded rather than assumed independent. The "
    "previous version's further statements — that the summary statistics are complete for every "
    "selected variant, that variants missing from either dataset were discarded rather than imputed, "
    "and that harmonisation left no variant unaligned — are **not carried** in the revised "
    "manuscript; the last of these cannot be checked against it at all, because harmonisation is no "
    "longer described (item 6b). |",
    "item6d",
)

# ── 16. item 6e：§2.7 → §2.8；单工具行的落点 → §3.4 ─────────────────────
sub(
    "| 6e | Multiple testing | Yes | §2.7 sets out a **three-level hierarchy**",
    "| 6e | Multiple testing | Yes | §2.8 sets out a **three-level hierarchy**",
    "item6e-anchor",
)
sub(
    "and both rest on a single variant (§3.5; Table 2 footnote). |",
    "and both rest on a single variant (§3.4). |",
    "item6e-tail",
)

# ── 17. item 7：方法仍在，结果见 item 12 ─────────────────────────────────
line(
    "| 7 |",
    "| 7 | Methods used to assess the assumptions | Yes | §2.7 (MR-Egger intercept for directional "
    "horizontal pleiotropy; Cochran's Q for heterogeneity; gene-level F for instrument strength, "
    "\"F>10 considered strong\"). The **methods** remain described, but their **results** are no "
    "longer reported — see item 12. |",
    "item7",
)

# ── 18. item 8：加权中位数/Egger 结果与共定位结果均未报告 ────────────────
line(
    "| 8 |",
    "| 8 | Sensitivity and additional analyses | Partial | §2.7 and §3.4; the full "
    "main-versus-stricter (r²≤0.01) comparison across all 18 gene–exposure combinations, with "
    "per-variant F statistics, is given in Table S3. Two analyses named in the previous version are "
    "**not reported in this revision**: the weighted-median and MR-Egger *results* (item 12), and "
    "the colocalization *results* — §2.7 still describes colocalization as performed (for VGLL4 and "
    "SDC4, and for SERPINE1 in this revision), but no colocalization posterior probability appears "
    "in the revised manuscript or in any of its supplementary tables, whereas the previous version "
    "reported them for two of the genes. |",
    "item8",
)

# ── 19. item 9a：§2.7 → §2.8，补 metafor 与 sandwich ─────────────────────
sub(
    "| 9a | Software and packages, with versions | Yes | §2.7 (R 4.6.1; DESeq2 1.52.0; limma "
    "3.68.5; clusterProfiler 4.20.0; STRINGdb 2.24.0;",
    "| 9a | Software and packages, with versions | Yes | §2.8 (R 4.6.1; DESeq2 1.52.0; limma "
    "3.68.5; clusterProfiler 4.20.0; metafor; sandwich; STRINGdb 2.24.0;",
    "item9a",
)

# ── 20. item 9b：新稿主语是 selection 不是 the panel ─────────────────────
sub(
    '| 9b | Pre-registration / protocol | Yes | §2.2: the panel "was',
    '| 9b | Pre-registration / protocol | Yes | §2.2: selection "was',
    "item9b",
)

# ── 21. item 10a：补 Figure 1 与 §3.1 ────────────────────────────────────
sub(
    "| 10a | Descriptive data — numbers at each stage, flow | Yes | §2.1; Table 1 |",
    "| 10a | Descriptive data — numbers at each stage, flow | Yes | §2.1 and §3.1; Table 1 and "
    "Figure 1 |",
    "item10a",
)

# ── 22. item 10b：MR 汇总统计未成表 ──────────────────────────────────────
line(
    "| 10b |",
    "| 10b | Summary statistics for exposures and outcomes | Partial | **Not tabulated in the "
    "revised submission.** The gene–exposure estimates under both pruning thresholds are in "
    "Table S3, but the per-variant summary statistics — SNP, effect/other allele, EAF, eQTL β (SE) "
    "and *P*, KOA β (SE), gene-level F — that the previous version gave in its Table 3 are not "
    "carried, and the previous version's summary figure (a forest plot of the MR estimates) is not "
    "part of this revision. Tables 2 and Figure 5 of the revised manuscript are the replication-"
    "cohort and single-cell-localisation items, not MR. |",
    "item10b",
)

# ── 23. item 10c：变体级特征未报告 ───────────────────────────────────────
line(
    "| 10c |",
    "| 10c | Genetic variant characteristics | Partial | Per-variant characteristics (SNP, "
    "effect/other allele, EAF, eQTL β (SE) and *P*, KOA β (SE)) are **not reported** in the revised "
    "submission. What is reported per gene–exposure set is the number of instruments retained, the "
    "IVW or Wald β (SE) with 95% CI and *P*, whether the instrument set changed under stricter "
    "pruning, and the per-SNP F values (Table S3). |",
    "item10c",
)

# ── 24. item 10d：来源人群/祖源不再写在该句 ──────────────────────────────
line(
    "| 10d |",
    "| 10d | Sample overlap between exposure and outcome samples | Yes | §2.7, third paragraph: "
    "the three eQTL datasets and the FinnGen outcome \"draw on non-overlapping sample frames, so no "
    "exposure–outcome sample overlap is expected\"; the same sentence states that the three eQTL "
    "datasets are **not** independent of one another, which is why they are compared rather than "
    "pooled. The source populations and ancestries (UK OA patients, United States GTEx donors, "
    "Estonian Biobank donors, Finnish FinnGen participants) are **not stated** in the revised "
    "manuscript; ancestry appears only in §4 (\"The genetic data are entirely European-ancestry\") "
    "and §2.7 gives 503 European-ancestry samples for the LD reference panel. |",
    "item10d",
)

# ── 25. item 11a：变体级关联未报告 ───────────────────────────────────────
line(
    "| 11a |",
    "| 11a | Variant–exposure and variant–outcome associations | Partial | **Not reported at "
    "variant level** in the revised submission (item 10c). Per-variant F statistics are given in "
    "Table S3; the per-variant eQTL and KOA effect estimates that the previous version tabulated are "
    "not carried. |",
    "item11a",
)

# ── 26. item 11b：估计值在 Table S3，无 MR 图 ────────────────────────────
line(
    "| 11b |",
    "| 11b | MR estimates with uncertainty | Yes | Table S3 gives β, SE, 95% CI and *P* for all 18 "
    "gene–exposure combinations under the main and stricter pruning thresholds. No MR figure forms "
    "part of this revision — the previous version's forest plot has been dropped. |",
    "item11b",
)

# ── 27. item 11d：无图 ───────────────────────────────────────────────────
line(
    "| 11d |",
    "| 11d | Plots | No | **No MR plot is included in this revision.** The previous version's "
    "forest plot of the MR estimates has been dropped, and the revised manuscript's five figures are "
    "the study design (Figure 1), the 57-candidate overview (Figure 2), per-cohort panel effect "
    "sizes (Figure 3), the cross-linker evidence matrix (Figure 4) and the single-cell localisation "
    "(Figure 5). The MR estimates are given numerically in Table S3 only. |",
    "item11d",
)

# ── 28. item 12：诊断结果已无落点，工具强度在 Table S3 ───────────────────
line(
    "| 12 |",
    "| 12 | Assessment of assumptions — results | Partial | **The diagnostic results are no longer "
    "reported.** §2.7 still states the methods, but the revised manuscript and its supplements carry "
    "**no MR-Egger intercept *P* value** and **no Cochran's Q value**; the only occurrences of "
    "\"Cochran\" in the manuscript are the two methodological mentions in §2.3 and §2.7, and no Q "
    "statistic is reported anywhere. What **is** reported is instrument strength: Table S3 gives "
    "per-SNP F for every retained variant and states that for multi-instrument sets the gene-level F "
    "is the arithmetic mean of those per-SNP values, defined only where two or more instruments "
    "survive (the single variant retained for VGLL4 under the macrophage exposure has F=62.3 on its "
    "own). The previous version's statement that no multi-instrument set in the macrophage exposure "
    "reaches the conventional threshold of 10 remains computable from Table S3, but is **not stated** "
    "in the revised manuscript. A single-variant F and a gene-level F are different quantities and "
    "are not compared here; per-SNP F values for every combination are listed in Table S3. |",
    "item12",
)

# ── 29. item 13：其他分析结果大多未承载 ──────────────────────────────────
line(
    "| 13 |",
    "| 13 | Other analyses — results | Partial | §3.4, in the paragraph headed \"Exploratory genetic "
    "evidence\". The revised manuscript states that \"no panel gene obtained transferable genetic "
    "support\" and names the two smallest *P* values — VGLL4 under the macrophage exposure "
    "(β=−0.036, *P*=6.5×10⁻⁵) and SERPINE1 under the tissue-matched synovial exposure (β=−0.039, "
    "*P*=0.002) — noting that both are single-variant Wald ratios for which \"no across-combination "
    "correction applies\" (§2.8). **Three items of the previous version are not carried:** (i) the "
    "count that six of the 18 gene–exposure combinations reach nominal significance, which remains "
    "computable from the main-analysis columns of Table S3 but is not stated in the text; (ii) the "
    "colocalization results, which are not reported in the revised manuscript or in any "
    "supplementary table although §2.7 still describes the colocalization as performed; and (iii) "
    "the statement that the weighted-median and MR-Egger results an unpruned analysis would have "
    "supported no longer exist, since no weighted-median or MR-Egger result is now reported at all. |",
    "item13",
)

# ── 30. item 15：五项编号限制 → §4 散文段 ────────────────────────────────
line(
    "| 15 |",
    "| 15 | Limitations | Yes | §4, the paragraph beginning \"Several limits bound these "
    "conclusions.\" The revised manuscript states its limits as continuous prose rather than as a "
    "numbered list, and the set of limits differs from the previous version's. (1) No matrix "
    "stiffness, cross-linking activity or protein abundance was measured, so every \"stiffness\" "
    "term names an annotated role. (2) No gene's expression could be assigned to a cell type: the "
    "single-cell dataset has no normal controls, only three donors, and a 43.65-percentage-point "
    "spread in fibroblast share, and the fibroblast-state association \"was tested against two null "
    "calibrations and **was not established**\" (§2.6). (3) The discovery cohorts are modest "
    "(74 samples) and demographic covariates are inconsistently annotated across public cohorts, so "
    "effect sizes are imprecise. (4) Machine learning does not support a diagnostic reading: the "
    "authoritative internal cross-validated AUC is 0.647 against a random gene-set null median of "
    "0.635 (*P*=0.438), the external 0.780 has no null of its own, and the earlier normalisation of "
    "GSE89408 had to be corrected (§2.1, §2.5, Table S2). (5) The MR analyses are exploratory and "
    "carry no causal weight: no P-value threshold was applied at instrument selection, six "
    "gene–exposure sets retain a single instrument, and the exposures are not interchangeable, since "
    "the genetic evidence for VGLL4 differs between sources (§3.4). (6) The genetic data are "
    "entirely European-ancestry, and the LD pruning threshold was calibrated on European samples. "
    "(7) The analysis is retrospective and re-uses published cohorts; it is not pre-registered, and "
    "the candidate set was assembled on knowledge grounds. |",
    "item15",
)

# ── 31. item 17：可推广性措辞回新稿 ──────────────────────────────────────
sub(
    "the genetic association evidence for VGLL4 differs between the eQTL sources used and causal "
    "support in synovial fibroblasts is not established, so the macrophage-derived result cannot be "
    "extended to synovial fibroblasts",
    "the genetic evidence for VGLL4 differs between the eQTL sources used, so the macrophage-derived "
    "result cannot be extended to synovial fibroblasts; §3.4 states that \"An instrument that "
    "predicts VGLL4 expression in macrophages is not thereby a valid instrument for VGLL4 "
    "expression in synovial fibroblasts\"",
    "item17",
)

# ── 32. item 20：数据共享逐字回新稿 ──────────────────────────────────────
line(
    "| 20 |",
    "| 20 | Data sharing | Yes | Declarations — \"GSE89408, GSE206848, GSE1919, GSE55235, "
    "GSE55457, GSE82107, and GSE152805 are publicly available from GEO\"; cis-eQTL data from eQTL "
    "Catalogue release 7 (QTD000519 / QTD000216 / QTD000001); FinnGen R12 knee osteoarthritis GWAS "
    "(M13_ARTHROSIS_KNEE) from the FinnGen repository; scripts \"available from the corresponding "
    "author upon reasonable request\". |",
    "item20",
)

# ══════════════ 后置条件 ══════════════════════════════════════════════════
def must(cond, msg):
    assert cond, f"后置条件失败: {msg}"
    log.append(f"  post OK  {msg}")


# 旧锚点全部消失
for bad in ["§3.5", "§2.5)", "Figure 5 (forest", "five numbered limitations",
            "Table 2; Figure 5", "Table 3 (SNP", "Supplementary Table S3", "PP.H4", "P=0.096"]:
    must(bad not in out, f"旧内容已清除: {bad!r}")

# 旧标题消失、新标题在
must("Cross-Cohort Expression Profile and Cellular Localization" not in out, "旧标题已移除")
must("Criterion-Dependent Priority of Matrix-Remodelling Candidates" in out, "新标题已写入")

# 新锚点分布（逐条核实后的期望计数）
must(out.count("§2.7") >= 10, f"§2.7 出现 {out.count('§2.7')} 次（MR 方法层，期望 ≥10）")
must(out.count("§2.8") >= 4, f"§2.8 出现 {out.count('§2.8')} 次（统计层，期望 ≥4）")
must(out.count("§3.4") >= 4, f"§3.4 出现 {out.count('§3.4')} 次（MR 结果层，期望 ≥4）")
must(out.count("§2.5") == 1, f"§2.5 出现 {out.count('§2.5')} 次（仅机器学习限制一处）")
must(out.count("§2.6") == 1, f"§2.6 出现 {out.count('§2.6')} 次（成纤维状态关联一处）")
must("Table S3" in out, "Table S3 已引用")

# 摘要引句用英式拼写（逐字取自新摘要）
must("Exploratory Mendelian randomisation gave no transferable genetic support" in out,
     "新摘要引句逐字在文件中")

# 条目行守恒：STROBE-MR 共 35 条主题行，一条不多一条不少（替换不增删条目）
item_re = re.compile(r"^\| (\d+[a-e]?) \| ")
n_items_src = sum(1 for l in src.split("\n") if item_re.match(l))
n_items_out = sum(1 for l in out.split("\n") if item_re.match(l))
must(n_items_src == 35, f"源条目行 {n_items_src}（期望 35）")
must(n_items_out == n_items_src, f"条目行守恒：源 {n_items_src} → 出 {n_items_out}")

# 行数：唯一增行来自来源注释块（空行 + 6 行注释 = +7），替换本身不增删行
must(len(out.split("\n")) == len(src.split("\n")) + 7,
     f"行数：源 {len(src.split(chr(10)))} → 出 {len(out.split(chr(10)))}（期望 +7，仅注释块）")
must(out.count("\n<!-- Source:") == 1, "来源注释块恰 1 处")

# 表体结构守恒：表格行数与源一致
must(sum(1 for l in out.split("\n") if l.startswith("| ")) ==
     sum(1 for l in src.split("\n") if l.startswith("| ")), "表格行数与源一致")

# 汉字：条目区（注释之外）必须一个都没有——英文产物；注释块是中文工作提示，提交前删除
zh = lambda s: set(re.findall(r"[一-鿿]", s))
strip_comments = lambda s: re.sub(r"<!--.*?-->", "", s, flags=re.S)
zh_body_out = zh(strip_comments(out))
must(zh_body_out == set(), f"条目区含汉字 {sorted(zh_body_out)}（英文产物须全英文）")
must(zh_body_out == zh(strip_comments(src)), "条目区的汉字集合与源一致（源条目区亦为空）")
must(zh(out) >= zh(src), "注释块中文未被删除（工作提示保留）")

# 节号守卫：清单里出现的每个 § 编号必须真实存在于新稿（防按偏移推编号）
# 新稿节标题（05_稿件/English_manuscript_MDPI_重定位2026-09-23.md 实测）：
#   §1 Introduction｜§2 Methods（§2.1 数据与预处理 … §2.8 统计）｜§3 Results（§3.1–§3.4）｜
#   §4 Discussion｜§5 Conclusion。§2 是 Methods 整层（item 9b 用「reported in §2」指该整层）。
SECS_NEW = {"1", "2", "2.1", "2.2", "2.3", "2.4", "2.5", "2.6", "2.7", "2.8",
            "3.1", "3.2", "3.3", "3.4", "4", "5"}
secs_out = set(re.findall(r"§(\d+(?:\.\d+)?)", out))
must(secs_out <= SECS_NEW, f"出现新稿不存在的节号 {sorted(secs_out - SECS_NEW)}")
SECS_EXPECTED = {"1", "2", "2.1", "2.2", "2.3", "2.5", "2.6", "2.7", "2.8", "3.1", "3.4", "4"}
must(secs_out == SECS_EXPECTED,
     f"节号均存在于新稿且与逐条核实一致：{sorted(secs_out)}")

# 陈旧表引用守卫：新稿 Table 3 是交联家族表，凡引 Table 3 处必须明写是「上一版」的。
# 非空性前提：源文件确有 Table 3 引用（旧稿的变体级表），否则本守卫是空转。
must(src.count("Table 3") >= 1,
     f"守卫非空转：源文件有 {src.count('Table 3')} 处旧 Table 3 引用可被清除")
_t3 = [m.start() for m in re.finditer(r"Table 3", out)]
for st in _t3:
    ctx = out[max(0, st - 90):st]
    must("previous version" in ctx, f"Table 3 引用须标明为上一版：…{ctx[-66:]!r}")
print(f"    陈旧表引用守卫：源 {src.count('Table 3')} 处 → 出 {len(_t3)} 处（均已标明上一版或已删除）")

# Reported 列分布（列 3）
codes = {}
for l in out.split("\n"):
    m = re.match(r"^\| (\d+[a-e]?) \| .*? \| ([^|]+?) \|", l)
    if m:
        codes[m.group(2).strip()] = codes.get(m.group(2).strip(), 0) + 1
log.append(f"  Reported 列分布: {codes}")
must(codes.get("Yes", 0) >= 18, "Yes ≥18 条")
must(codes.get("No", 0) == 1, "No 恰 1 条（item 11d）")
must(codes.get("Partial", 0) >= 8, "Partial ≥8 条")
must(codes.get("Not done", 0) == 1, "Not done 恰 1 条（item 11c）")

DST.parent.mkdir(parents=True, exist_ok=True)
DST.write_text(out, encoding="utf-8")

print("C1 STROBE-MR 重映射 —— 完成")
print(f"源  : {SRC}  ({len(src.encode())} B, {len(src.split(chr(10)))} 行)")
print(f"出  : {DST}  ({len(out.encode())} B, {len(out.split(chr(10)))} 行)")
print(f"源 md5: {hashlib.md5(src.encode()).hexdigest()}")
print(f"出 md5: {hashlib.md5(out.encode()).hexdigest()}")
print()
for l in log:
    print(l)
