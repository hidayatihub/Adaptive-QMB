# Tests and preflight checks

## Static repository audit

From the repository root:

```r
source("tests/AUDIT_REPOSITORY.R")
```

This check does not require the private/local input CSVs. If no CSV inputs are present, the input-count check is reported as `SKIP`. If inputs are present, exactly 48 paired CSV files are expected.

## One-combination pilot

After the local paired inputs are available:

```r
source("tests/RUN_ANALYSIS_PILOT_v1.0.1.R")
```

The pilot executes one matched combination at `B = 30` and writes temporary/generated output under `results/pilot_analysis/`, which is ignored by Git.
