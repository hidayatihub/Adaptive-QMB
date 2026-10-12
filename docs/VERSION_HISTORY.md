# Version history

Repository name: **Adaptive-QMB**. Scientific method name: **ACT-QMB** (Adaptive Cut-off Tail Quantile Matching Bootstrap).

Repository-level version numbers describe the publication package and execution workflow. They should not be confused with the locked internal ACT-QMB scientific-method identifier reported by the adaptive core.

## v1.2.4 — Adaptive-QMB public release candidate

Active publication-facing files:

```text
RUN_FULL_ANALYSIS_v1.2.4.R
RUN_INFERENCE_v1.2.1.R
PLOT_CODE_v1.0.1.R
QMB_ADAPTIVE_QDMT_CORE_v1.2.1.R
QMB_KOKIC_CORE.R
COMMON_METRIC_CORE_v1.1.0.R
CHRIPS_GCM_CORE_v1.1.0.R
```

Repository packaging changes for public release:

- aligned README and audit scripts with the actual active filenames;
- retained only the latest active inference runner (`v1.2.1`) in `code/`;
- retained only the latest pilot runner (`v1.0.1`) in `tests/`;
- removed RStudio state, macOS metadata, paired input CSVs, and generated pilot outputs from the public package;
- documented the six-arm major-revision design;
- added citation metadata, license-status notice, SHA-256 code hashes, and a lightweight GitHub Actions audit;
- harmonized the full-runner provenance label with the public filename (`v1.2.4-audited`);
- updated stale textual references in the inference and plotting runners to point to the current canonical runner versions.

The packaging changes above do not alter the statistical formulas or correction algorithms.

### Locked scientific method identity

The adaptive core reports:

```text
QMB_ACT_METHOD_ID = QMB_CV_ADAPTIVE_QDM_TAIL_5Q_STAGE3
QMB_ACT_VERSION   = v1.1.0_STAGE3_20260907
```

This internal method identity is intentionally preserved.

### Inference v1.2.1

`RUN_INFERENCE_v1.2.1.R` supersedes `v1.2.0`. The v1.2.1 runner fixes recursive Stage-3 daily-file discovery by filtering the normalized full path for each `B_####` directory after filename discovery, and it performs a fail-fast 48-cell/TaiESM1 source audit for every stochastic method and bootstrap size.

## Earlier repository drafts

Earlier local/repository drafts used older runner labels and contained superseded runner copies. They are not included in the clean public package to avoid ambiguity about the canonical execution path.
