# Input data

The public code repository does **not** distribute the paired CHIRPS–CMIP6 input CSV files.

For local reproduction, place exactly 48 paired CSV files under:

```text
data/adaptive/
```

Subdirectories are allowed. The canonical design is:

```text
3 rainfall regimes × 4 GCMs × 4 SSPs = 48 files
```

Each file should contain at least `date`, `chirps`, and `gcm`, cover the required historical/evaluation periods, and have a filename containing `paired` plus regime/GCM/SSP identifiers.

The repository `.gitignore` excludes these CSV files by default. Redistribution of source-derived data remains subject to the terms of the original data providers.
