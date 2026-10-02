# Analysis code: Criterion-Explicit Reproducibility of Matrix-Remodelling Candidates in Osteoarthritis Synovium

This repository holds the analysis code for the manuscript

> **Criterion-Explicit Reproducibility of Matrix-Remodelling Candidates in Osteoarthritis Synovium**

It is released so that every number and every figure in the paper can be traced to the script that
produced it. The code is archived on Zenodo; the DOI is given in the manuscript's Data Availability
Statement and in `CITATION.cff`.

## What is in scope

Everything here is a Python or an R script. The rule is: **every `.py` and `.R` file** under the
four directories below.

| Directory | Files | What it does |
| --- | --- | --- |
| `code/analysis/` | 37 | Per-cohort differential expression, merged effects and leave-one-cohort-out re-ranking, the cross-linking family under one criterion framework, donor-resolved single-cell localisation, and the MR instrument checks |
| `code/mr_coloc/` | 10 | Instrument extraction, LD pruning, colocalisation and multi-instrument Mendelian randomisation for the *cis*-eQTL assessment |
| `code/figures/` | 6 | The five main figures (English-label set — the set that drew the delivered figures) |
| `code/audit/` | 8 | The scripts that machine-verify the numbers: frozen-input invariant checks, figure-legend structure checks against the figures themselves, and the fixed-point criterion |

## What is deliberately not in scope

- **The manuscript-revision scripts.** These are the scripts whose names begin `F`, `H`, `G` or `R0`
  in the authors' working tree. They edit the manuscript, the supplementary tables and the
  submission package; they produce no analytical number, and including them would suggest that the
  reported results came out of them.
- **The Chinese-label figure code.** The delivered figures were drawn by the English-label set
  included here. A second, Chinese-label variant of the same four figure scripts exists in the
  authors' tree; the two sets are not identical, so only the English-label set is claimed.
- **The data.** All primary data are public and are redistributed by their repositories, not here.
  The cohorts used are GSE89408, GSE206848, GSE1919, GSE55235, GSE55457, GSE82107. The single-cell
  data are GSE152805. The *cis*-eQTL sources are Alasoo and GTEx.
- **The frozen inputs and result tables** that the scripts read. These are the paper's source data
  and supplementary tables.

## Reproducing the analysis

The scripts are numbered by stage. Read them in this order: `A*` (sample and preprocessing audit,
per-cohort effects), `B*` (merged effects, leave-one-out, single-cell localisation, figures),
`C*` (criterion reconciliation and the STROBE-MR checks), `D*` (supplementary table S4 and the
figures), `E*` (figure assembly). Scripts under `code/mr_coloc/` are numbered independently and are
run before the S3 tables are built.

Each script names its own inputs at the top of the file. Paths were written for the authors'
working tree and are not portable as-is; the input a script names is the file it was run against.
The primary data must be downloaded from GEO and the eQTL sources from their respective portals
before the scripts are run.

## Integrity

`MANIFEST.sha256` lists the SHA-256 of every file in this repository, two spaces between the digest
and the path. Verify from the repository root with:

```
shasum -a 256 -c MANIFEST.sha256
```

Every number reported in the manuscript was recomputed from the delivered source data before
submission; `code/audit/` holds the scripts that perform those cross-checks.

## License

Released under the MIT License — see `LICENSE`.

## Citation

See `CITATION.cff`. If you use this code, please cite the manuscript and the archived release.
