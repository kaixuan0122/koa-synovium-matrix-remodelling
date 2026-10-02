# -*- coding: utf-8 -*-
"""MR 与 coloc 共用的暴露数据配置与基因坐标表。

为什么单独成文件：基因坐标原先在 extract_multi_iv.py 和 extract_coloc.py 里各抄了一份，
后来只改了前者——SDC4 在 multi_iv 里是 GRCh38 的 45,325,288，在 coloc 里还是 GRCh37 的
43,952,000，于是 SDC4 的共定位窗口几乎没盖住基因本身（实测只取到该基因 cis-eQTL 的约 6%，
剩下 94% 全在窗口之外），正文"SDC4 未共定位"的结论因此没有依据。
一张表，杜绝这类漂移。

关于抓取窗口：实测 eQTL Catalogue 的 cis 区比 ±500 kb 宽（VGLL4 的 cis-eQTL 从 11.0 Mb
一直延伸到 12.7 Mb 以外）。±500 kb 因此是一个有意的、写入 Methods 的方法学选择，
不是"取到了全部"。脚本会打印每个基因命中行相对窗口边界的位置：撞边即说明窗口是限制因素，
解读时需谨慎。没有把窗口放宽到能容下整个 cis 区，是因为 EBI 对 sumstats 限速约
11 KB/s，一个 ±2 Mb 窗口要 1000 秒以上，14 次抓取不可行。
"""

# eQTL Catalogue r7（= GTEx v8）。
# ⚠️ 同一登录号在 r8 pre-release 已被换成 GTEx v10（QTD000216 样本数 483 → 659），
#    本文用的是 r7。引用时必须写明版本，否则日后按登录号复现会拿到另一套数据。
EXPOSURES = {
    "GTEx": dict(
        label="GTEx fibroblast (eQTL Catalogue r7, dataset QTD000216)",
        n=483,
        url="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
            "QTS000015/QTD000216/QTD000216.all.tsv.gz",
    ),
    "Alasoo": dict(
        label="Alasoo 2018 macrophage (eQTL Catalogue r7, dataset QTD000001)",
        n=84,
        url="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
            "QTS000001/QTD000001/QTD000001.all.tsv.gz",
    ),
    # 唯一对口的组织：人 OA 膝关节滑膜。Steinberg 2020, Nat Commun, PMID 33637762。
    # 样本量 72，故功效低；但它测的就是本文所研究的组织，其余两套都不是
    # （GTEx fibroblast 是皮肤成纤维细胞系，Alasoo 是单核细胞来源巨噬细胞）。
    "Synovium": dict(
        label="Steinberg 2020 synovium (eQTL Catalogue r7, dataset QTD000519)",
        n=72,
        url="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
            "QTS000028/QTD000519/QTD000519.all.tsv.gz",
    ),
}
PRIMARY = "Synovium"      # 主分析：与研究对象同组织
SENSITIVITY = ["GTEx", "Alasoo"]   # 敏感性分析：跨细胞类型稳健性

# 基因位置（GRCh38）。CYR61 原写 86,046,000，偏 +465 kb，已按 CCN1 的 GRCh38 坐标更正。
# （不过 CCN1/CYR61 已按正文 §2.4 的说明从本报告中撤出，此处仅保留以备复核。）
GENES = {
    "VGLL4":    ("ENSG00000144560", "3",  11556068),
    "SDC4":     ("ENSG00000124145", "20", 45325288),   # 原 coloc 误用 GRCh37 的 43,952,000
    "LOXL1":    ("ENSG00000129038", "15", 73926461),
    "SERPINE1": ("ENSG00000106366", "7",  101127103),
    "SAV1":     ("ENSG00000151748", "14", 50633579),
    "ITGA5":    ("ENSG00000161638", "12", 54395260),
    "CYR61":    ("ENSG00000142871", "1",  85580761),
}
FETCH = 500_000

MR_GENES = ["VGLL4", "SDC4", "LOXL1", "SERPINE1", "SAV1", "ITGA5"]
# SERPINE1 于 2026-09-17 补入（§5.4）：坐标与 cis-eQTL 缓存本就在，只是从未跑过 coloc。
# 补入不改变下面 FinnGen 的缓存键——GENES_NEEDED = MR_GENES ∪ COLOC_GENES，而 MR_GENES
# 已含 SERPINE1，故 need 集合不变，不会触发 819 MB sumstats 的重扫。
COLOC_GENES = ["VGLL4", "SDC4", "SERPINE1"]     # CYR61 已撤出

FINNGEN = "mr_data/finngen_KOA.gz"     # 结局 GWAS（相对本目录）
# FinnGen R12 M13_ARTHROSIS_KNEE：官方 manifest 载明 61,356 例 / 315,115 对照。
# 与本地 sumstats 反推的病例占比 0.1630（29 万个 SNP，四分位 ±0.0002）吻合到 3 位有效数字。
FINNGEN_CASES = 61356
FINNGEN_CONTROLS = 315115
