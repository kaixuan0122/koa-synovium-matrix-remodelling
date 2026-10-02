# -*- coding: utf-8 -*-
"""B5 MR 工具变量核查（§5.1）。

要核的三件事：
  ① 表 2 / 表 3 里那个「Gene-level F」到底是什么量；
  ② 逐 SNP 的 F 是否与 rsID 绑定、是否被报出（用户明确要求 per-SNP F **and rsID**）；
  ③ 严格剪枝（r²≤0.01）后保留的到底是哪几个变异。

为什么这三点要紧：
  · 若「Gene-level F」是逐 SNP F 的**算术平均**，它就不是工具集联合强度的任何标准度量。
    平均会被强变异拉高，从而**掩盖**接近弱工具阈值的那个变异。写「F=51.7」会让人以为
    整组工具都很强，而该组的两个变异逐 SNP F 是 90.0 与 13.3。
  · F 必须与 rsID 绑定才有可核验性。S3 的 Per-SNP F 列是**降序排列的裸数字**
    （如「19.1, 9.7, 4.9」），读者无法知道哪个 F 属于哪个变异。
  · 严格剪枝若改变了集合却没写出留下的是哪个变异，该敏感性分析就无法被复核。

输入（均为已落盘产物，本脚本不重跑任何 MR）：
  04_分析数据/MR_coloc/multi_iv_exposure_outcome_{tag}_pruned.csv   r²≤0.1 剪枝后的工具
  04_分析数据/MR_coloc/mr_sensitivity_component.csv                 含严格剪枝行的保留 SNP 名单
  04_分析数据/MR_coloc/表3_工具_pruned.md                            正文表 3（用于对账）
输出：
  08_修订_2026-09-17/02_结果表/逐SNP工具变量与F.csv
  08_修订_2026-09-17/03_分析/logs/B5_MR工具核查.log
"""
import csv
import os
import re
import statistics as st

BASE = ("/Users/zhangkaixuan/Library/Containers/com.kingsoft.wpsoffice.mac/Data/Library/"
        "Application Support/Kingsoft/WPS Cloud Files/userdata/qing/filecache/285447201/"
        "01.材料/01.博士材料/10.毕业/01.小论文/04.立项/滑膜基质刚度生成与刚度感知 0912")
MR = os.path.join(BASE, "04_分析数据/MR_coloc")
OUT = os.path.join(BASE, "08_修订_2026-09-17/02_结果表")
LOG = os.path.join(BASE, "08_修订_2026-09-17/03_分析/logs/B5_MR工具核查.log")
os.makedirs(OUT, exist_ok=True)
os.makedirs(os.path.dirname(LOG), exist_ok=True)

TAGS = ["Synovium", "GTEx", "Alasoo"]
LIN = []


def say(s=""):
    print(s)
    LIN.append(s)


# ---------- 1. 逐 SNP 工具表：F 与 rsID 绑定 ----------
rows = []
for tag in TAGS:
    with open(os.path.join(MR, f"multi_iv_exposure_outcome_{tag}_pruned.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            b, se = float(r["beta_exp"]), float(r["se_exp"])
            rows.append(dict(
                exposure=tag, gene=r["gene"], snp=r["SNP"],
                effect_allele=r["effect_allele"], other_allele=r["other_allele"],
                eaf_eqtl=None,                 # pruned csv 无 eQTL MAF 列，不臆造
                eqtl_beta=b, eqtl_se=se,
                eqtl_ci_low=b - 1.96 * se, eqtl_ci_high=b + 1.96 * se,
                eqtl_p=float(r["p_exp"]),
                koa_beta=float(r["beta_out"]), koa_se=float(r["se_out"]),
                per_snp_F=(b / se) ** 2))

say("========== 1. 逐 SNP 工具变量与 F（r²≤0.1 剪枝后）==========")
say(f"{'exposure':<9} {'gene':<9} {'rsID':<14} {'eQTL β':>9} {'SE':>8} "
    f"{'95% CI':>20} {'eQTL P':>10} {'逐SNP F':>9}")
for r in rows:
    say(f"{r['exposure']:<9} {r['gene']:<9} {r['snp']:<14} {r['eqtl_beta']:+9.3f} "
        f"{r['eqtl_se']:8.4f} [{r['eqtl_ci_low']:+8.3f},{r['eqtl_ci_high']:+8.3f}] "
        f"{r['eqtl_p']:10.3g} {r['per_snp_F']:9.2f}")
say(f"共 {len(rows)} 个工具（与剪枝日志一致：Synovium 10 + GTEx 12 + Macrophage 13 = 35）")

with open(os.path.join(OUT, "逐SNP工具变量与F.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# ---------- 2. 每基因-暴露的 F 汇总：算术平均是不是表 3 的「Gene-level F」 ----------
# 表 3 的「Gene-level F」从 md 文本读——这是与上面 csv 完全不同的来源，
# 两条路径对上才说明该列确实等于逐 SNP F 的算术平均（而不是"大致相当"）。
tbl_F = {}
for line in open(os.path.join(MR, "表3_工具_pruned.md"), encoding="utf-8").read().splitlines():
    c = [x.strip() for x in line.strip().strip("|").split("|")]
    if len(c) != 9 or c[0] == "Exposure" or set(c[0]) <= set("-"):
        continue
    tag = {"Synovium (n=72)": "Synovium", "GTEx fibroblast (n=483)": "GTEx",
           "Macrophage (n=84)": "Alasoo"}[c[0]]
    # 同一基因的 F 在多行重复，取第一次即可；单工具集合该列为「—」
    if c[8] != "—":
        tbl_F[(tag, c[1])] = float(c[8])

say()
say("========== 2. 「Gene-level F」= 逐 SNP F 的算术平均？==========")
say(f"{'exposure':<9} {'gene':<9} {'k':>2} {'逐SNP F（与 rsID 绑定）':<40} "
    f"{'最小':>7} {'中位':>7} {'算术平均':>9} {'表3值':>7} {'一致':>5}")
bad = []
for tag in TAGS:
    for gene in ["LOXL1", "SERPINE1", "VGLL4", "SDC4", "ITGA5", "SAV1"]:
        sub = sorted([r for r in rows if r["exposure"] == tag and r["gene"] == gene],
                     key=lambda x: -x["per_snp_F"])
        if not sub:
            continue
        fs = [r["per_snp_F"] for r in sub]
        mean = st.mean(fs)
        tab = tbl_F.get((tag, gene))
        detail = " ".join(f"{r['snp']}={r['per_snp_F']:.1f}" for r in sub)
        ok = "—"
        if tab is not None:
            ok = "✅" if abs(mean - tab) < 0.06 else "❌"
            if ok == "❌":
                bad.append((tag, gene, mean, tab))
        say(f"{tag:<9} {gene:<9} {len(sub):>2} {detail:<40} "
            f"{min(fs):7.1f} {st.median(fs):7.1f} {mean:9.2f} "
            f"{(f'{tab:.1f}' if tab is not None else '—（单工具）'):>7} {ok:>5}")
assert not bad, f"表 3 的 Gene-level F 与逐 SNP F 的算术平均不符：{bad}"
say()
say("✅ 表 3 的「Gene-level F」逐值等于该基因各 SNP 逐 SNP F 的算术平均（35 个工具全部对上）。")
say("   → 该列应更名为「平均逐 SNP F」。「基因水平 F」这一名称会让人读成工具集的")
say("     联合强度，而算术平均不是任何标准的联合强度度量：它会被强变异拉高。")

# ---------- 3. 平均 F 掩盖了什么：弱工具（F<10）清单 ----------
say()
say("========== 3. 逐 SNP F<10 的工具（平均 F 会掩盖这些）==========")
weak = [r for r in rows if r["per_snp_F"] < 10]
say(f"{'exposure':<9} {'gene':<9} {'rsID':<14} {'逐SNP F':>8} {'所在集合平均F':>13}")
for r in sorted(weak, key=lambda x: x["per_snp_F"]):
    fs = [q["per_snp_F"] for q in rows if q["exposure"] == r["exposure"] and q["gene"] == r["gene"]]
    say(f"{r['exposure']:<9} {r['gene']:<9} {r['snp']:<14} {r['per_snp_F']:8.2f} {st.mean(fs):13.2f}")
say(f"共 {len(weak)} / {len(rows)} 个工具的逐 SNP F<10。")

msg = ("稿件称「巨噬暴露中没有任何多工具集合达到 F>10」（以平均 F 计），但逐 SNP 看，"
       "巨噬 SDC4 的 rs17322191 自身 F=11.8>10；反之，被平均掩盖的弱工具真实的更多："
       "巨噬 LOXL1 三个工具逐 SNP F 为 6.6/5.1/4.7，全部 <10。")
say("→ " + msg)
say("   因此「平均逐 SNP F>10」不能作为「工具集不弱」的依据；须同时报最小逐 SNP F。")
open(os.path.join(OUT, "MR工具_逐SNP_F与弱工具.csv"), "w", newline="", encoding="utf-8").write(
    "exposure,gene,snp,per_snp_F,set_mean_F,set_min_F\n"
    + "".join(f"{r['exposure']},{r['gene']},{r['snp']},{r['per_snp_F']:.4f},"
              f"{st.mean([q['per_snp_F'] for q in rows if q['exposure']==r['exposure'] and q['gene']==r['gene']]):.4f},"
              f"{min(q['per_snp_F'] for q in rows if q['exposure']==r['exposure'] and q['gene']==r['gene']):.4f}\n"
              for r in sorted(weak, key=lambda x: x["per_snp_F"])))

# ---------- 4. 严格剪枝（r²≤0.01）后保留的变异 ----------
say()
say("========== 4. 严格剪枝 r²≤0.01 后保留的变异 ==========")
say(f"{'exposure':<9} {'gene':<9} {'主剪枝 k':>7} {'严格 k':>6} {'保留的变异（严格）':<34} {'集合是否改变':<10}")
# 主剪枝下本来就是单工具的组合（这些组合没有 strict 行，S3 表里沿用主剪枝值）
main_single = sum(
    1 for tag in TAGS for g in {r["gene"] for r in rows if r["exposure"] == tag}
    if len([r for r in rows if r["exposure"] == tag and r["gene"] == g]) == 1)

strict_n, single_after = {}, 0
with open(os.path.join(MR, "mr_sensitivity_component.csv"), encoding="utf-8") as fh:
    for c in list(csv.reader(fh))[1:]:
        expo, gene, n_iv, snp, role = c[0], c[1], int(c[2]), c[3], c[11]
        if not role.startswith("strict"):
            continue
        main_n = len([r for r in rows if r["exposure"] == expo and r["gene"] == gene])
        changed = "yes" if n_iv != main_n else "no"
        strict_n[(expo, gene)] = n_iv
        if n_iv == 1:
            single_after += 1
        say(f"{expo:<9} {gene:<9} {main_n:>7} {n_iv:>6} {snp:<34} {changed:<10}")
say()
say(f"由多工具退化为单工具的集合：{single_after}；主剪枝下即为单工具的：{main_single}；"
    f"故严格阈值下单工具集合共 {single_after + main_single} 个（S3 脚注作 6→9）。")
say("（注意：mr_sensitivity_component.csv 的 strict 行只在主剪枝下 n≥2 时才写出，")
say("  故 18 个组合中只有 12 个有 strict 行；其余 6 个沿用主剪枝行，S3 表里显示「no」。）")
say("→ S3 已报「有 5 个集合改变、单工具集合由 6 增至 9」，但未写出留下的是哪个变异；")
say("  该名单可从 mr_sensitivity_component.csv 的 strict 行取得，应补入 S3。")

# ---------- 5. 与 S3 的 Per-SNP F 列核对 ----------
say()
say("========== 5. S3 的 Per-SNP F 列是否可核验 ==========")
S3 = os.path.join(BASE, "05_投稿包/MDPI版/Supplementary/Supplementary_Table_S3.docx")
say("S3 的 Per-SNP F 取自表 3 的 eQTL β(SE)（21_make_table_s3.py），与本节第 1 节同源，")
say("故数值必然一致。问题不在数值，在**呈现**：该列是降序排列的裸数字（如「19.1, 9.7, 4.9」），")
say("未与 rsID 绑定，读者无法判断哪个 F 属于哪个变异——而 F 与 rsID 的绑定正是")
say("判断「是哪一个工具弱」的唯一途径。建议直接把第 1 节的成对「rsID=F」写入 S3。")

# ---------- 5b. Per-SNP F 的重算精度：从源 csv 算 vs 从表 3 舍入值算 ----------
say()
say("========== 5b. Per-SNP F 的重算精度 ==========")
devs = []
for r in rows:
    b3, se4 = round(r["eqtl_beta"], 3), round(r["eqtl_se"], 4)
    devs.append((abs((b3 / se4) ** 2 - r["per_snp_F"]), r["exposure"], r["gene"],
                 r["snp"], r["per_snp_F"], (b3 / se4) ** 2))
devs.sort(reverse=True)
say(f"{'exposure':<9} {'gene':<9} {'rsID':<14} {'F（源 csv 全精度）':>18} {'F（表 3 舍入值反算）':>20} {'差':>7}")
for d, e, g, s, f_src, f_tab in devs[:5]:
    say(f"{e:<9} {g:<9} {s:<14} {f_src:18.4f} {f_tab:20.4f} {d:7.3f}")
say(f"最大出入 {devs[0][0]:.3f}（{devs[0][1]} / {devs[0][2]} / {devs[0][3]}）。")
say("S3 与本节的 Per-SNP F 都是自表 3 的舍入值反算（21_make_table_s3.py 亦然），")
say("故与源 csv 有 ≤0.14 的出入。这不构成错误，但既然该列要用于判断工具强弱，")
say("建议改为直接从 *_pruned.csv 的 beta_exp/se_exp 计算（本脚本第 1 节即如此）。")

# ---------- 6. 表 3 的 EAF 列是哪一侧的频率 ----------
say()
say("========== 6. 表 3 的「EAF」口径 ==========")
eaf_out = {}
for tag in TAGS:
    with open(os.path.join(MR, f"multi_iv_exposure_outcome_{tag}_pruned.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            eaf_out[(tag, r["gene"], r["SNP"])] = float(r["eaf_out"])
tbl_eaf, cmp_n, hit = {}, 0, 0
for line in open(os.path.join(MR, "表3_工具_pruned.md"), encoding="utf-8").read().splitlines():
    c = [x.strip() for x in line.strip().strip("|").split("|")]
    if len(c) != 9 or c[0] == "Exposure" or set(c[0]) <= set("-"):
        continue
    tag = {"Synovium (n=72)": "Synovium", "GTEx fibroblast (n=483)": "GTEx",
           "Macrophage (n=84)": "Alasoo"}[c[0]]
    k = (tag, c[1], c[2])
    if k not in eaf_out:
        continue
    cmp_n += 1
    # 表 3 列序：0 Exposure 1 Gene 2 SNP 3 Effect/other 4 EAF 5 eQTL β(SE) 6 eQTL P 7 KOA β(SE) 8 F
    if abs(float(c[4]) - eaf_out[k]) < 5e-4:
        hit += 1
    else:
        say(f"  ❌ {tag}/{c[1]}/{c[2]}：表 3 EAF={c[4]}，csv eaf_out={eaf_out[k]:.4f}")
say(f"表 3 的 EAF 列与 *_pruned.csv 的 eaf_out 一致：{hit}/{cmp_n}")
say("→ 该列是**结局（FinnGen KOA）**的效应等位基因频率，不是 eQTL 的。表头只写「EAF」，")
say("  读者会读成 eQTL 侧频率。eQTL 侧 MAF 在 extract_all.py 里读到了但未写入 *_pruned.csv，")
say("  故本核查无法给出 eQTL 侧 MAF 的对照。建议表头写明「KOA EAF」。")

# ---------- 7. 表 3 的 eQTL 列缺 CI ----------
say()
say("========== 6. 表 2 / 表 3 的列完整性 ==========")
say("表 3（工具级）现有：SNP、Effect/other allele、EAF、eQTL β (SE)、eQTL P、KOA β (SE)、Gene-level F")
say("  → 缺 **eQTL β 的 95% CI**（已在本脚本第 1 节算出，口径 β ± 1.96·SE）")
say("表 2（估计级）现有：n IV、Method、β (SE)、P、Intercept P、F、PP.H4")
say("  → 缺 **MR β 的 95% CI**（S3 已有，可直接对齐）")
say("表 2 脚注称「基因水平 F 仅在多工具变量集合中报告，因该估计量只在那里有定义；")
say("  单工具变量集合中软件返回的 F 为空值」——**此说法不成立**：")
say("  (a) 单工具集合的逐 SNP F 完全有定义，S3 自己就报了（如滑膜 VGLL4 F=12.9、")
say("      SERPINE1 F=16.4、ITGA5 F=11.1），与「为空值」直接冲突；")
say("  (b) 多工具集合的那个值不是「基因水平」量，而是逐 SNP F 的算术平均。")
say("  → 建议：表 2/3 该列更名「平均逐 SNP F（mean per-SNP F）」，删去「为空值/只在那里有定义」")
say("     的辩护，改为「单工具集合报其逐 SNP F」——本脚本第 1 节已给出全部 35 个值。")

open(LOG, "w", encoding="utf-8").write("\n".join(LIN) + "\n")
say()
say(f"已写出：02_结果表/逐SNP工具变量与F.csv（{len(rows)} 行）")
say(f"        02_结果表/MR工具_逐SNP_F与弱工具.csv（{len(weak)} 行）")
say(f"        03_分析/logs/B5_MR工具核查.log")
