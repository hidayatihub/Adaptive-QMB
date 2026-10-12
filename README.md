# Adaptive-QMB

**Adaptive Cut-off Tail Quantile Matching Bootstrap (ACT-QMB) for CMIP6 daily rainfall bias correction**

Repository release candidate: **v1.2.4**

**Repository name:** `Adaptive-QMB`  
**Scientific method:** Adaptive Cut-off Tail Quantile Matching Bootstrap (`ACT-QMB`)

This repository contains the publication-facing R implementation of the ACT-QMB major-revision workflow. The package is organized around one canonical analysis path, one common metric engine, an inferential comparison layer, and reproducible figure generation.

> **Important:** the paired CHIRPS–CMIP6 input CSV files and large generated outputs are intentionally not distributed in this code repository. Users must place the required paired inputs under `data/adaptive/` before running the analysis.

## Scientific design

- Calibration/source period: **1990–2014**
- Independent evaluation period: **2015–2025**
- Stage-1 quantiles: **Q10, Q50, Q90, Q95, Q99**
- Wet-day definition: **rainfall > 0.1 mm/day**
- Leap days: removed
- Candidate ACT cutoffs: **0.900–0.950 in steps of 0.005**
- Cutoff selection: blocked pseudo-future cross-validation entirely within **1990–2014**
- Stage-3 bootstrap sizes: **30, 100, 250, 500, 1000, 2000**
- Tail metrics: pooled by calendar month across the 2015–2025 evaluation period
- R quantile estimator: **type 7**

The independent 2015–2025 evaluation observations are not used to select the adaptive cutoff.

## Comparison design

The major-revision workflow evaluates six matched arms:

1. **RAW GCM** — deterministic benchmark
2. **QDM** — deterministic benchmark
3. **QMB-LIN**
4. **Hard P90**
5. **Fixed P94**
6. **ACT-QMB**

QMB-LIN, Hard P90, Fixed P94, and ACT-QMB have paired Stage-3 bootstrap ensembles. RAW GCM and QDM are deterministic and are repeated unchanged at each Stage-3 evaluation layer only to provide matched comparisons; they are not bootstrapped.

The canonical input design is:

```text
3 rainfall regimes × 4 GCMs × 4 SSPs = 48 paired combinations
```

Rainfall regimes:

- Monsoonal
- Equatorial
- Local

CMIP6 models:

- CESM2-WACCM
- CMCC-ESM2
- MRI-ESM2-0
- TaiESM1

SSP pathways:

- SSP126
- SSP245
- SSP370
- SSP585

Project routing used by the analysis:

```text
MRI-ESM2-0  -> location-2
CESM2-WACCM -> location-1
CMCC-ESM2   -> location-1
TaiESM1     -> location-1
```

## Canonical code

```text
code/
├── RUN_FULL_ANALYSIS_v1.2.4.R
├── RUN_INFERENCE_v1.2.1.R
├── PLOT_CODE_v1.0.1.R
├── QMB_ADAPTIVE_QDMT_CORE_v1.2.1.R
├── QMB_KOKIC_CORE.R
├── COMMON_METRIC_CORE_v1.1.0.R
└── CHRIPS_GCM_CORE_v1.1.0.R
```

Roles:

- `RUN_FULL_ANALYSIS_v1.2.4.R` — canonical 48-combination Stage-1/Stage-2/CV/Stage-3 runner.
- `QMB_ADAPTIVE_QDMT_CORE_v1.2.1.R` — ACT-QMB scientific core.
- `QMB_KOKIC_CORE.R` — audited QMB baseline implementation.
- `COMMON_METRIC_CORE_v1.1.0.R` — locked common evaluation definitions for the six-arm comparison.
- `RUN_INFERENCE_v1.2.1.R` — common comparison and inferential analysis, including all-B Stage-3 source auditing.
- `PLOT_CODE_v1.0.1.R` — publication figure generation.
- `CHRIPS_GCM_CORE_v1.1.0.R` — independent RAW GCM diagnostic utility; the matched inferential workflow obtains RAW GCM from the same paired inputs used by the correction methods.

### Scientific-core provenance

The adaptive core reports:

```text
Method ID        : QMB_CV_ADAPTIVE_QDM_TAIL_5Q_STAGE3
Internal version : v1.1.0_STAGE3_20260907
```

The public file `QMB_ADAPTIVE_QDMT_CORE_v1.2.1.R` is the audited repository-facing implementation that preserves this locked internal scientific method identity.

## Repository layout

```text
Adaptive-QMB/
├── code/
├── data/
│   ├── README.md
│   └── adaptive/
├── results/
│   └── README.md
├── tests/
├── docs/
├── archive/
├── .github/workflows/
├── .gitignore
├── .gitattributes
├── README.md
├── CITATION.cff
├── LICENSE_PENDING.md
└── Adaptive-QMB.Rproj
```

## Software requirements

The analysis is written in R. The core ACT-QMB/QMB workflow is primarily based on base R. Additional packages are used by supporting utilities:

- `ggplot2` — required for publication figures and RAW diagnostic plotting;
- `dplyr` — used by selected RAW diagnostic summaries;
- `scales` — optional enhancement used by the plotting script when available.

Package availability can be checked locally before the full run.

## Input data

Place the 48 paired daily CSV files under `data/adaptive/`. Subdirectories are allowed and recommended, for example:

```text
data/adaptive/
├── Monsoonal/
├── Equatorial/
└── Local/
```

Each input file must contain at least:

```text
date
chirps
gcm
```

and filenames should contain `paired` plus the regime, GCM, and SSP identifiers, for example:

```text
paired_Monsoonal_CESM2_WACCM_ssp126.csv
```

The input data are not committed to this repository. See `data/README.md`.

## Preflight audit

Open `Adaptive-QMB.Rproj`, set the repository root as the working directory, and run:

```r
source("tests/AUDIT_REPOSITORY.R")
```

The audit checks required files, parses the R scripts, verifies citation placeholders, and checks the input count when local input data are present. A public clone without input CSVs is valid; a local reproduction run should contain exactly 48 paired CSV files.

For a one-combination execution test after local inputs are available:

```r
source("tests/RUN_ANALYSIS_PILOT_v1.0.1.R")
```

## Canonical run order

After placing the 48 paired input files under `data/adaptive/`:

```r
# 1. Full ACT-QMB workflow
source("code/RUN_FULL_ANALYSIS_v1.2.4.R")

# 2. Six-arm common comparison and inference
source("code/RUN_INFERENCE_v1.2.1.R")

# 3. Publication figures
source("code/PLOT_CODE_v1.0.1.R")
```

## Principal generated outputs

The full runner writes to:

```text
results/full_analysis/
```

The inferential runner writes to:

```text
results/major_revision_inference/
```

The plotting runner writes to:

```text
results/figures/
```

Large outputs are intentionally excluded from Git by default. This repository is the code and reproducibility workflow; it is not a redistribution archive for the underlying CHIRPS/CMIP6 paired data.

## Reproducibility and integrity

- The analysis uses repository-relative paths.
- Input discovery is restricted to `data/adaptive/`.
- The full runner enforces the expected 48-combination design.
- The inferential runner v1.2.1 includes a fail-fast audit for Stage-3 daily files at every bootstrap size.
- `docs/CODE_SHA256SUMS.txt` records SHA-256 hashes for the publication-facing code files in this release candidate.
- GitHub Actions performs an R syntax/repository audit on pushes and pull requests.

See `docs/VERSION_HISTORY.md` for version provenance.

## Citation

Citation metadata are provided in `CITATION.cff`.

Repository URL:

<https://github.com/hidayatihub/Adaptive-QMB>

A repository/software DOI is not claimed until a permanent archive such as Zenodo is created.

## License status

A software license has **not yet been granted** for this release candidate. The repository may be viewed for scientific review and reproducibility assessment, but reuse, modification, or redistribution should not be assumed until a license is explicitly added.

See `LICENSE_PENDING.md`.
