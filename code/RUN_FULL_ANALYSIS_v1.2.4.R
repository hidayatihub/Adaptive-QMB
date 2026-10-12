# =============================================================================
# RUN_FULL_ANALYSIS_v1.2.4.R
# QMB Adaptive QDM-Tail with Cross-Validated Cutoff + Stage-3 Bootstrap
# Publication-facing reproducibility runner for the ACT-QMB study
# Repository release: v1.2.4-audited
# =============================================================================
#
# CANONICAL WORKFLOW
#   DATA
#     -> QMB Stage 1 (1990-2014)
#     -> QMB-LIN Stage 2
#     -> blocked pseudo-future CV for adaptive QDM-tail cutoff
#     -> Adaptive Stage 2
#     -> paired QMB-LIN / Hard-P90 / Fixed-P94 / ACT-QMB Stage 3 bootstrap
#     -> deterministic QDM reference
#     -> final accuracy, P95/P99, coverage, and interval-width outputs
#
# EXPECTED REPOSITORY STRUCTURE
#
# qmb-adaptive-interpolation/
# |-- code/
# |   |-- RUN_FULL_ANALYSIS_v1.2.4.R
# |   |-- QMB_ADAPTIVE_QDMT_CORE_v1.2.1.R
# |   `-- QMB_KOKIC_CORE.R
# |-- data/
# |   `-- paired input CSV files (subfolders allowed)
# `-- results/
#
# RUN
#   From the repository root:
#       source("code/RUN_FULL_ANALYSIS_v1.2.4.R")
#
# Or:
#       Rscript code/RUN_FULL_ANALYSIS_v1.2.4.R
#
# IMPORTANT
# - The adaptive core is the publication-facing filename for the scientifically
#   unchanged internal implementation formerly named
#   QMB_ADAPTIVE_QDMT_CV_v1.1.0_STAGE3.R.
# - QDM is implemented inside QMB_ADAPTIVE_QDMT_CORE_v1.2.1.R and is the QDM reference
#   used by the adaptive algorithm and final paper comparison.
# - The standalone QDM code is therefore not required by this canonical run.
# =============================================================================

options(stringsAsFactors = FALSE)

RUNNER_RELEASE <- "v1.2.4-audited"
EXPECTED_N_COMBINATIONS <- 48L
STRICT_EXPECTED_DESIGN <- TRUE

# -----------------------------------------------------------------------------
# 0. Resolve repository root without hard-coded user paths
# -----------------------------------------------------------------------------
.get_this_script <- function() {
  # Rscript --file=
  ca <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", ca, value = TRUE)
  if (length(hit)) {
    f <- sub("^--file=", "", hit[1])
    if (nzchar(f)) return(normalizePath(f, winslash = "/", mustWork = FALSE))
  }

  # source(...)
  nf <- sys.nframe()
  if (nf > 0L) {
    for (i in rev(seq_len(nf))) {
      f <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
      if (!is.null(f) && length(f) == 1L && nzchar(f)) {
        return(normalizePath(f, winslash = "/", mustWork = FALSE))
      }
    }
  }

  NA_character_
}

THIS_SCRIPT <- .get_this_script()
SCRIPT_DIR <- if (is.na(THIS_SCRIPT)) {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
} else {
  dirname(THIS_SCRIPT)
}

# If runner lives in /code, repository root is its parent.
REPO_ROOT <- if (tolower(basename(SCRIPT_DIR)) == "code") {
  dirname(SCRIPT_DIR)
} else {
  SCRIPT_DIR
}
REPO_ROOT <- normalizePath(REPO_ROOT, winslash = "/", mustWork = FALSE)

CODE_DIR <- file.path(REPO_ROOT, "code")
DATA_ROOT <- file.path(REPO_ROOT, "data", "adaptive")
RESULTS_ROOT <- file.path(REPO_ROOT, "results")
OUTPUT_ROOT <- file.path(RESULTS_ROOT, "full_analysis")

ADAPTIVE_CORE_FILE <- file.path(CODE_DIR, "QMB_ADAPTIVE_QDMT_CORE_v1.2.1.R")
KOKIC_CORE_FILE <- file.path(CODE_DIR, "QMB_KOKIC_CORE.R")

# Robust provenance path for the runner itself.
# When source() no longer exposes sys.frame()$ofile at the provenance step,
# THIS_SCRIPT can be NA. Resolve the versioned repository filename explicitly
# instead of falling back to the obsolete RUN_FULL_ANALYSIS.R name.
RUNNER_FILE <- if (!is.na(THIS_SCRIPT) && file.exists(THIS_SCRIPT)) {
  normalizePath(THIS_SCRIPT, winslash = "/", mustWork = TRUE)
} else {
  candidate <- file.path(CODE_DIR, "RUN_FULL_ANALYSIS_v1.2.4.R")
  if (file.exists(candidate)) {
    normalizePath(candidate, winslash = "/", mustWork = TRUE)
  } else {
    NA_character_
  }
}

# -----------------------------------------------------------------------------
# 1. Locked scientific configuration
# -----------------------------------------------------------------------------
DATE_COL <- "date"
OBS_COL <- "chirps"
GCM_COL <- "gcm"

BASE_YEARS <- c(1990L, 2014L)
TARGET_YEARS <- c(2015L, 2025L)
WET_THRESHOLD <- 0.1
NO_LEAP <- TRUE

RANDOM_STATE <- 42L

# Adaptive QDM-tail cutoff selection
CV_TOLERANCE_MM <- 1.0
RUN_STRICT_05MM_SENSITIVITY <- TRUE

# Stage 3
RUN_STAGE3_BOOTSTRAP <- TRUE
B_VALUES <- c(30L, 100L, 250L, 500L, 1000L, 2000L)
STAGE3_INTERVAL <- 0.95
SAVE_STAGE3_DAILY_ALL_B <- TRUE
REQUIRE_STAGE3_DAILY_ALL_B <- TRUE
KEEP_STAGE3_ENSEMBLE_RDS <- FALSE

# IMPORTANT:
# Stage 2 is deterministic. It has point predictions/diagnostic metrics but no
# bootstrap prediction interval. Prediction-interval coverage, width and
# probabilistic calibration are Stage-3 quantities only.

# -----------------------------------------------------------------------------
# 2. Preflight
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("QMB ADAPTIVE INTERPOLATION — FULL REPRODUCIBILITY RUN\n")
cat("Repository release :", RUNNER_RELEASE, "\n")
cat("Repository root    :", REPO_ROOT, "\n")
cat("Code directory     :", CODE_DIR, "\n")
cat("Data directory     :", DATA_ROOT, "\n")
cat("Results directory  :", OUTPUT_ROOT, "\n")
cat("============================================================\n\n")

required_files <- c(ADAPTIVE_CORE_FILE, KOKIC_CORE_FILE)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    paste0(
      "Required code file(s) not found:\n - ",
      paste(missing_files, collapse = "\n - "),
      "\n\nExpected files:\n",
      "  code/QMB_ADAPTIVE_QDMT_CORE_v1.2.1.R\n",
      "  code/QMB_KOKIC_CORE.R"
    ),
    call. = FALSE
  )
}

if (!dir.exists(DATA_ROOT)) {
  stop(
    paste0(
      "Data directory not found:\n  ", DATA_ROOT,
      "\n\nCreate 'data/adaptive' under the repository root and place ",
      "the paired CSV inputs there."
    ),
    call. = FALSE
  )
}

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 3. Load publication core
# -----------------------------------------------------------------------------
source(ADAPTIVE_CORE_FILE)

needed_objects <- c(
  "QMB_ACT_VERSION", "QMB_ACT_METHOD_ID", "QMB_ACT_TAUS", "QMB_ACT_CUTOFFS"
)
missing_objects <- needed_objects[
  !vapply(needed_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_objects)) {
  stop(
    "Adaptive core loaded but expected object(s) are missing: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

needed_functions <- c(
  "default_location_id",
  "validate_manifest",
  "run_qmb_adaptive_qdm_tail"
)
missing_functions <- needed_functions[
  !vapply(needed_functions, exists, logical(1), mode = "function", inherits = TRUE)
]
if (length(missing_functions)) {
  stop(
    "Adaptive core loaded but expected function(s) are missing: ",
    paste(missing_functions, collapse = ", "),
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 4. Discover and construct the 3-regime x 4-GCM x 4-SSP manifest
# -----------------------------------------------------------------------------
files <- list.files(
  DATA_ROOT,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
files <- files[grepl("paired", basename(files), ignore.case = TRUE)]

if (!length(files)) {
  stop(
    paste0(
      "No paired CSV files were found under:\n  ", DATA_ROOT,
      "\n\nInput filenames must contain the word 'paired'."
    ),
    call. = FALSE
  )
}

parse_one <- function(path) {
  txt <- toupper(gsub("[- ]", "_", path))

  # Scientific regime labels plus legacy geographical filenames/folders.
  regime <- if (grepl("MONSOONAL|JAWA_TIMUR|EAST_JAVA", txt)) {
    "Monsoonal"
  } else if (grepl("EQUATORIAL|SUMATERA_BARAT|WEST_SUMATRA", txt)) {
    "Equatorial"
  } else if (grepl("LOCAL|MALUKU", txt)) {
    "Local"
  } else {
    NA_character_
  }

  gcm <- if (grepl("CESM2_WACCM", txt)) {
    "CESM2_WACCM"
  } else if (grepl("MRI_ESM2_0", txt)) {
    "MRI_ESM2_0"
  } else if (grepl("TAIESM1", txt)) {
    "TaiESM1"
  } else if (grepl("CMCC_ESM2", txt)) {
    "CMCC_ESM2"
  } else {
    NA_character_
  }

  ssp <- if (grepl("SSP126", txt)) {
    "ssp126"
  } else if (grepl("SSP245", txt)) {
    "ssp245"
  } else if (grepl("SSP370", txt)) {
    "ssp370"
  } else if (grepl("SSP585", txt)) {
    "ssp585"
  } else {
    NA_character_
  }

  if (is.na(regime) || is.na(gcm) || is.na(ssp)) return(NULL)

  data.frame(
    REGIME = regime,
    LOCATION_ID = default_location_id(regime, gcm),
    GCM = gcm,
    SSP = ssp,
    FILE = normalizePath(path, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

manifest_list <- lapply(files, parse_one)
manifest_list <- manifest_list[
  !vapply(manifest_list, is.null, logical(1))
]

if (!length(manifest_list)) {
  stop(
    paste0(
      "Paired CSV files were found, but none could be parsed into ",
      "REGIME/GCM/SSP combinations."
    ),
    call. = FALSE
  )
}

MANIFEST <- do.call(rbind, manifest_list)
MANIFEST <- validate_manifest(MANIFEST)

cat("Detected analysis manifest:\n")
print(
  MANIFEST[, c("REGIME", "LOCATION_ID", "GCM", "SSP", "FILE")],
  row.names = FALSE
)
cat("\nDetected combinations:", nrow(MANIFEST), "\n")

if (isTRUE(STRICT_EXPECTED_DESIGN) &&
    nrow(MANIFEST) != EXPECTED_N_COMBINATIONS) {
  stop(
    sprintf(
      paste0(
        "Final paper design requires %d combinations ",
        "(3 regimes x 4 GCM x 4 SSP), but %d were detected."
      ),
      EXPECTED_N_COMBINATIONS,
      nrow(MANIFEST)
    ),
    call. = FALSE
  )
}

# Preserve the manifest detected by the publication runner.
utils::write.csv(
  MANIFEST,
  file.path(OUTPUT_ROOT, "RUNNER_INPUT_MANIFEST.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5. Record run configuration before computation
# -----------------------------------------------------------------------------
config_lines <- c(
  paste0("RUNNER_RELEASE=", RUNNER_RELEASE),
  paste0("CORE_METHOD_ID=", QMB_ACT_METHOD_ID),
  paste0("CORE_INTERNAL_VERSION=", QMB_ACT_VERSION),
  paste0("BASE_YEARS=", paste(BASE_YEARS, collapse = "-")),
  paste0("TARGET_YEARS=", paste(TARGET_YEARS, collapse = "-")),
  paste0("TAUS=", paste(QMB_ACT_TAUS, collapse = ",")),
  paste0("WET_THRESHOLD=", WET_THRESHOLD),
  paste0("RANDOM_STATE=", RANDOM_STATE),
  paste0("CV_CUTOFFS=", paste(sprintf("%.3f", QMB_ACT_CUTOFFS), collapse = ",")),
  paste0("CV_TOLERANCE_MM=", CV_TOLERANCE_MM),
  paste0("STRICT_05MM_SENSITIVITY=", RUN_STRICT_05MM_SENSITIVITY),
  paste0("STAGE3_BOOTSTRAP=", RUN_STAGE3_BOOTSTRAP),
  paste0("STAGE3_B_VALUES=", paste(B_VALUES, collapse = ",")),
  paste0("STAGE3_INTERVAL=", STAGE3_INTERVAL),
  paste0("SAVE_STAGE3_DAILY_ALL_B=", SAVE_STAGE3_DAILY_ALL_B),
  paste0("REQUIRE_STAGE3_DAILY_ALL_B=", REQUIRE_STAGE3_DAILY_ALL_B),
  paste0("QDM_REFERENCE=embedded_empirical_calendar_month_QDM_in_adaptive_core"),
  paste0("N_COMBINATIONS=", nrow(MANIFEST))
)

writeLines(
  config_lines,
  file.path(OUTPUT_ROOT, "RUN_CONFIGURATION.txt")
)

# -----------------------------------------------------------------------------
# 6. Run the complete paper analysis
# -----------------------------------------------------------------------------
run_started <- Sys.time()

RESULT <- run_qmb_adaptive_qdm_tail(
  manifest = MANIFEST,
  output_root = OUTPUT_ROOT,
  core_file = KOKIC_CORE_FILE,
  date_col = DATE_COL,
  obs_col = OBS_COL,
  gcm_col = GCM_COL,
  base_years = BASE_YEARS,
  target_years = TARGET_YEARS,
  no_leap = NO_LEAP,
  wet_threshold = WET_THRESHOLD,
  cutoffs = QMB_ACT_CUTOFFS,
  cv_tolerance_mm = CV_TOLERANCE_MM,
  also_run_strict_05mm = RUN_STRICT_05MM_SENSITIVITY,
  random_state = RANDOM_STATE,
  run_stage3_bootstrap = RUN_STAGE3_BOOTSTRAP,
  b_values = B_VALUES,
  stage3_interval = STAGE3_INTERVAL,
  save_stage3_daily_all_b = SAVE_STAGE3_DAILY_ALL_B,
  keep_stage3_ensemble_rds = KEEP_STAGE3_ENSEMBLE_RDS,
  save_rds = TRUE
)

run_finished <- Sys.time()
elapsed_seconds <- as.numeric(
  difftime(run_finished, run_started, units = "secs")
)

# -----------------------------------------------------------------------------
# 7. Reproducibility metadata
# -----------------------------------------------------------------------------
repro_dir <- file.path(OUTPUT_ROOT, "reproducibility")
dir.create(repro_dir, recursive = TRUE, showWarnings = FALSE)

capture.output(
  sessionInfo(),
  file = file.path(repro_dir, "sessionInfo.txt")
)

writeLines(
  c(
    paste0("RUN_STARTED=", format(run_started, "%Y-%m-%d %H:%M:%S %Z")),
    paste0("RUN_FINISHED=", format(run_finished, "%Y-%m-%d %H:%M:%S %Z")),
    paste0("ELAPSED_SECONDS=", format(elapsed_seconds, scientific = FALSE)),
    paste0("RUNNER_RELEASE=", RUNNER_RELEASE),
    paste0("CORE_METHOD_ID=", QMB_ACT_METHOD_ID),
    paste0("CORE_INTERNAL_VERSION=", QMB_ACT_VERSION)
  ),
  file.path(repro_dir, "RUN_PROVENANCE.txt")
)

# SHA-256 for the three canonical code files, using system utilities when
# available. This does not affect the scientific computation.
.write_sha256 <- function(files, outfile) {
  files <- unique(files[!is.na(files) & nzchar(files)])
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    warning(
      "Skipping missing provenance file(s): ",
      paste(basename(missing), collapse = ", "),
      call. = FALSE
    )
  }
  files <- files[file.exists(files)]
  if (!length(files)) {
    warning("No existing files available for SHA-256 provenance.", call. = FALSE)
    return(invisible(FALSE))
  }
  files <- normalizePath(files, winslash = "/", mustWork = TRUE)

  sha256sum <- Sys.which("sha256sum")
  shasum <- Sys.which("shasum")

  if (nzchar(sha256sum)) {
    rows <- unlist(lapply(files, function(f) {
      system2(sha256sum, shQuote(f), stdout = TRUE)
    }))
    writeLines(rows, outfile)
    return(invisible(TRUE))
  }

  if (nzchar(shasum)) {
    rows <- unlist(lapply(files, function(f) {
      system2(shasum, c("-a", "256", shQuote(f)), stdout = TRUE)
    }))
    writeLines(rows, outfile)
    return(invisible(TRUE))
  }

  warning(
    "No SHA-256 command found; SHA256SUMS.txt was not generated.",
    call. = FALSE
  )
  invisible(FALSE)
}

.write_sha256(
  c(
    ADAPTIVE_CORE_FILE,
    KOKIC_CORE_FILE,
    RUNNER_FILE
  ),
  file.path(repro_dir, "SHA256SUMS.txt")
)



# -----------------------------------------------------------------------------
# 7a. Paper-ready six-arm comparison including RAW GCM
# -----------------------------------------------------------------------------
# RAW is a deterministic benchmark, exactly like QDM. It is evaluated from the
# same matched 2015-2025 daily records already written by the adaptive core.
# It is NOT assigned prediction-interval metrics because it has no Stage-3
# bootstrap ensemble.

.paper_pct_improvement <- function(method_value, raw_value) {
  ifelse(
    is.finite(method_value) & is.finite(raw_value) & raw_value != 0,
    100 * (raw_value - method_value) / raw_value,
    NA_real_
  )
}

stage2_daily_files <- list.files(
  file.path(OUTPUT_ROOT, "FINAL_STAGE2"),
  pattern = "stage2_daily_qmb_adaptive_qdm_tail\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(stage2_daily_files) != EXPECTED_N_COMBINATIONS) {
  stop(
    "RAW paper-table build failed: expected ", EXPECTED_N_COMBINATIONS,
    " Stage-2 daily files, found ", length(stage2_daily_files), ".",
    call. = FALSE
  )
}

raw_daily_parts <- vector("list", length(stage2_daily_files))
raw_month_parts <- vector("list", length(stage2_daily_files))

for (i in seq_along(stage2_daily_files)) {
  f <- stage2_daily_files[i]
  d <- utils::read.csv(f, check.names = FALSE)

  required_raw_cols <- c("date", "obs", "gcm", "REGIME", "LOCATION_ID", "GCM", "SSP")
  miss <- setdiff(required_raw_cols, names(d))
  if (length(miss)) {
    stop(
      "RAW paper-table build failed in ", basename(f),
      ": missing column(s): ", paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  d$date <- as.Date(d$date)
  if (nrow(d) != 4015L) {
    stop(
      "RAW paper-table build failed: ", f,
      " contains ", nrow(d), " rows; expected 4015.",
      call. = FALSE
    )
  }

  raw_daily_parts[[i]] <- d[, required_raw_cols, drop = FALSE]

  mm <- calendar_month_metrics(
    daily = d,
    pred_col = "gcm",
    obs_col = "obs"
  )
  mm$REGIME <- d$REGIME[1]
  mm$LOCATION_ID <- d$LOCATION_ID[1]
  mm$GCM <- d$GCM[1]
  mm$SSP <- d$SSP[1]
  mm$METHOD <- "RAW_GCM"
  raw_month_parts[[i]] <- mm
}

raw_daily_all <- do.call(rbind, raw_daily_parts)
raw_month_all <- do.call(rbind, raw_month_parts)

raw_overall <- data.frame(
  METHOD = "RAW_GCM",
  DAILY_MAE = .act_mae(raw_daily_all$obs, raw_daily_all$gcm),
  DAILY_RMSE = .act_rmse(raw_daily_all$obs, raw_daily_all$gcm),
  DAILY_BIAS = .act_bias(raw_daily_all$obs, raw_daily_all$gcm),
  MEAN_MONTHLY_P95_ABS_ERROR = mean(raw_month_all$P95_ABS_ERROR, na.rm = TRUE),
  MEAN_MONTHLY_P99_ABS_ERROR = mean(raw_month_all$P99_ABS_ERROR, na.rm = TRUE),
  stringsAsFactors = FALSE
)

# ---- Stage 2: six-arm paper table --------------------------------------------
paper_stage2 <- rbind(
  raw_overall,
  RESULT$final_overall_comparison
)

paper_stage2$METHOD <- as.character(paper_stage2$METHOD)
paper_stage2$METHOD[paper_stage2$METHOD == "QMB"] <- "QMB_LIN"
paper_stage2$METHOD[paper_stage2$METHOD == "ADAPTIVE_CV"] <- "ACT_QMB"

paper_order <- c("RAW_GCM", "QMB_LIN", "HARD90", "FIXED94", "ACT_QMB", "QDM")
paper_stage2 <- paper_stage2[
  match(paper_order, paper_stage2$METHOD),
  , drop = FALSE
]

raw_mae <- raw_overall$DAILY_MAE[1]
raw_rmse <- raw_overall$DAILY_RMSE[1]
raw_abs_bias <- abs(raw_overall$DAILY_BIAS[1])
raw_p95 <- raw_overall$MEAN_MONTHLY_P95_ABS_ERROR[1]
raw_p99 <- raw_overall$MEAN_MONTHLY_P99_ABS_ERROR[1]

paper_stage2$ABS_DAILY_BIAS <- abs(paper_stage2$DAILY_BIAS)
paper_stage2$MAE_IMPROVEMENT_VS_RAW_PCT <-
  .paper_pct_improvement(paper_stage2$DAILY_MAE, raw_mae)
paper_stage2$RMSE_IMPROVEMENT_VS_RAW_PCT <-
  .paper_pct_improvement(paper_stage2$DAILY_RMSE, raw_rmse)
paper_stage2$ABS_BIAS_IMPROVEMENT_VS_RAW_PCT <-
  .paper_pct_improvement(paper_stage2$ABS_DAILY_BIAS, raw_abs_bias)
paper_stage2$P95_IMPROVEMENT_VS_RAW_PCT <-
  .paper_pct_improvement(paper_stage2$MEAN_MONTHLY_P95_ABS_ERROR, raw_p95)
paper_stage2$P99_IMPROVEMENT_VS_RAW_PCT <-
  .paper_pct_improvement(paper_stage2$MEAN_MONTHLY_P99_ABS_ERROR, raw_p99)

# RAW is the reference: report 0 rather than tiny floating-point noise.
paper_stage2[
  paper_stage2$METHOD == "RAW_GCM",
  grep("_IMPROVEMENT_VS_RAW_PCT$", names(paper_stage2))
] <- 0

utils::write.csv(
  paper_stage2,
  file.path(OUTPUT_ROOT, "PAPER_STAGE2_6ARM_COMPARISON.csv"),
  row.names = FALSE
)

utils::write.csv(
  raw_month_all,
  file.path(OUTPUT_ROOT, "PAPER_RAW_MONTHLY_METRICS_BY_COMBINATION.csv"),
  row.names = FALSE
)

# ---- Stage 3: six-arm all-B paper table --------------------------------------
# Four stochastic methods come from Stage 3. QDM and RAW are deterministic
# benchmarks and therefore repeat unchanged across B; their interval metrics
# are intentionally NA.
if (isTRUE(RUN_STAGE3_BOOTSTRAP)) {
  s3 <- RESULT$stage3_overall_comparison

  raw_s3 <- do.call(
    rbind,
    lapply(B_VALUES, function(B) {
      data.frame(
        METHOD = "RAW_GCM_DETERMINISTIC",
        N_BOOT = B,
        DAILY_MAE = raw_overall$DAILY_MAE,
        DAILY_RMSE = raw_overall$DAILY_RMSE,
        DAILY_BIAS = raw_overall$DAILY_BIAS,
        MEAN_MONTHLY_P95_ABS_ERROR = raw_overall$MEAN_MONTHLY_P95_ABS_ERROR,
        MEAN_MONTHLY_P99_ABS_ERROR = raw_overall$MEAN_MONTHLY_P99_ABS_ERROR,
        COVERAGE = NA_real_,
        MEAN_INTERVAL_WIDTH = NA_real_,
        stringsAsFactors = FALSE
      )
    })
  )

  paper_stage3 <- rbind(s3, raw_s3)

  # Stable publication labels.
  paper_stage3$METHOD <- as.character(paper_stage3$METHOD)
  paper_stage3$METHOD[paper_stage3$METHOD == "QMB_LIN_STAGE3"] <- "QMB_LIN_STAGE3"
  paper_stage3$METHOD[paper_stage3$METHOD == "HARD90_STAGE3"] <- "HARD90_STAGE3"
  paper_stage3$METHOD[paper_stage3$METHOD == "FIXED94_STAGE3"] <- "FIXED94_STAGE3"
  paper_stage3$METHOD[paper_stage3$METHOD == "ADAPTIVE_CV_STAGE3"] <- "ACT_QMB_STAGE3"

  paper_stage3$ABS_DAILY_BIAS <- abs(paper_stage3$DAILY_BIAS)
  paper_stage3$MAE_IMPROVEMENT_VS_RAW_PCT <-
    .paper_pct_improvement(paper_stage3$DAILY_MAE, raw_mae)
  paper_stage3$RMSE_IMPROVEMENT_VS_RAW_PCT <-
    .paper_pct_improvement(paper_stage3$DAILY_RMSE, raw_rmse)
  paper_stage3$ABS_BIAS_IMPROVEMENT_VS_RAW_PCT <-
    .paper_pct_improvement(paper_stage3$ABS_DAILY_BIAS, raw_abs_bias)
  paper_stage3$P95_IMPROVEMENT_VS_RAW_PCT <-
    .paper_pct_improvement(paper_stage3$MEAN_MONTHLY_P95_ABS_ERROR, raw_p95)
  paper_stage3$P99_IMPROVEMENT_VS_RAW_PCT <-
    .paper_pct_improvement(paper_stage3$MEAN_MONTHLY_P99_ABS_ERROR, raw_p99)

  paper_stage3[
    paper_stage3$METHOD == "RAW_GCM_DETERMINISTIC",
    grep("_IMPROVEMENT_VS_RAW_PCT$", names(paper_stage3))
  ] <- 0

  method_rank <- c(
    RAW_GCM_DETERMINISTIC = 1L,
    QDM_DETERMINISTIC = 2L,
    QMB_LIN_STAGE3 = 3L,
    HARD90_STAGE3 = 4L,
    FIXED94_STAGE3 = 5L,
    ACT_QMB_STAGE3 = 6L
  )
  paper_stage3$.METHOD_ORDER <- unname(method_rank[paper_stage3$METHOD])
  paper_stage3 <- paper_stage3[
    order(paper_stage3$N_BOOT, paper_stage3$.METHOD_ORDER),
    , drop = FALSE
  ]
  paper_stage3$.METHOD_ORDER <- NULL

  utils::write.csv(
    paper_stage3,
    file.path(OUTPUT_ROOT, "PAPER_STAGE3_6ARM_ALL_B_COMPARISON.csv"),
    row.names = FALSE
  )
}

# Keep paper-ready tables in RESULT for the final console report.
RESULT$paper_stage2_6arm_comparison <- paper_stage2
if (isTRUE(RUN_STAGE3_BOOTSTRAP)) {
  RESULT$paper_stage3_6arm_all_b_comparison <- paper_stage3
}

# -----------------------------------------------------------------------------
# 7b. Major-revision Stage-3 daily-output audit
# -----------------------------------------------------------------------------
if (isTRUE(RUN_STAGE3_BOOTSTRAP) && isTRUE(REQUIRE_STAGE3_DAILY_ALL_B)) {

  expected_methods_daily <- c(
    "qmb_lin_stage3_daily_interval.csv",
    "hard90_stage3_daily_interval.csv",
    "fixed94_stage3_daily_interval.csv",
    "adaptive_cv_stage3_daily_interval.csv"
  )

  s3_root <- file.path(OUTPUT_ROOT, "FINAL_STAGE3")
  s3_daily <- list.files(
    s3_root,
    pattern = "_stage3_daily_interval\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  expected_n <- 48L * length(B_VALUES) * length(expected_methods_daily)

  cat("\n============================================================\n")
  cat("STAGE-3 DAILY ALL-B AUDIT\n")
  cat("============================================================\n")
  cat("Expected daily files :", expected_n, "\n")
  cat("Found daily files    :", length(s3_daily), "\n")

  if (length(s3_daily) != expected_n) {
    stop(
      "Stage-3 daily all-B audit failed: expected ", expected_n,
      " files (=48 combinations x ", length(B_VALUES),
      " B values x 4 stochastic methods), found ", length(s3_daily),
      ". Do not run inferential analysis until this is corrected.",
      call. = FALSE
    )
  }

  # Audit every B independently.
  b_count <- vapply(B_VALUES, function(B) {
    token <- sprintf("/B_%04d/", B)
    sum(grepl(token, gsub("\\\\", "/", s3_daily), fixed = TRUE))
  }, integer(1))

  cat("\nFiles by B:\n")
  for (i in seq_along(B_VALUES)) {
    cat(sprintf("  B=%-4d : %3d / 192 %s\n",
                B_VALUES[i], b_count[i],
                if (b_count[i] == 192L) "PASS" else "FAIL"))
  }

  if (any(b_count != 192L)) {
    stop(
      "Stage-3 daily all-B audit failed for B values: ",
      paste(B_VALUES[b_count != 192L], collapse = ", "),
      call. = FALSE
    )
  }

  # Audit the four stochastic Stage-3 methods separately.
  method_count <- vapply(expected_methods_daily, function(fn) {
    sum(basename(s3_daily) == fn)
  }, integer(1))

  cat("\nFiles by stochastic Stage-3 method:\n")
  for (i in seq_along(expected_methods_daily)) {
    cat(sprintf("  %-39s : %3d / 288 %s\n",
                expected_methods_daily[i], method_count[i],
                if (method_count[i] == 288L) "PASS" else "FAIL"))
  }

  if (any(method_count != 288L)) {
    stop(
      "Stage-3 daily method audit failed. Expected 288 files per stochastic ",
      "method (=48 combinations x 6 B values).",
      call. = FALSE
    )
  }

  # Check that every daily file contains the complete 2015-2025 no-leap target.
  nrow_daily <- vapply(s3_daily, function(f) {
    nrow(utils::read.csv(f, check.names = FALSE))
  }, integer(1))

  cat("\nDaily row-count range :", min(nrow_daily), "-", max(nrow_daily), "\n")
  if (any(nrow_daily != 4015L)) {
    bad <- s3_daily[nrow_daily != 4015L]
    stop(
      "Stage-3 daily row-count audit failed: ",
      length(bad), " file(s) do not contain 4015 evaluation days.",
      call. = FALSE
    )
  }

  # Summary outputs must also span all six B values.
  summary_checks <- list(
    accuracy = file.path(OUTPUT_ROOT, "FINAL_STAGE3_ACCURACY_BY_COMBINATION.csv"),
    interval = file.path(OUTPUT_ROOT, "FINAL_STAGE3_INTERVAL_METRICS_BY_COMBINATION.csv"),
    monthly = file.path(OUTPUT_ROOT, "FINAL_STAGE3_MONTHLY_METRICS.csv"),
    calibration = file.path(OUTPUT_ROOT, "MASTER_STAGE3_PROBABILISTIC_CALIBRATION.csv")
  )

  for (nm in names(summary_checks)) {
    f <- summary_checks[[nm]]
    if (!file.exists(f)) {
      stop("Missing Stage-3 summary output: ", basename(f), call. = FALSE)
    }
    d <- utils::read.csv(f, check.names = FALSE)
    if (!"N_BOOT" %in% names(d)) {
      stop("Missing N_BOOT in ", basename(f), call. = FALSE)
    }
    got_b <- sort(unique(as.integer(d$N_BOOT[is.finite(as.numeric(d$N_BOOT))])))
    if (!identical(got_b, sort(as.integer(B_VALUES)))) {
      stop(
        "Stage-3 summary B audit failed in ", basename(f),
        ". Found B={", paste(got_b, collapse = ","), "}; expected B={",
        paste(B_VALUES, collapse = ","), "}.",
        call. = FALSE
      )
    }
  }

  cat("\nSTAGE-3 DAILY ALL-B AUDIT: PASS\n")
  cat("  48 combinations x 6 B x 4 stochastic methods = 1152 daily files\n")
  cat("  Every daily file = 4015 evaluation days\n")
  cat("  Accuracy/interval/monthly/calibration summaries contain all B values\n")
  cat("============================================================\n\n")
}

# -----------------------------------------------------------------------------
# 8. Final console report
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("FULL ANALYSIS COMPLETED\n")
cat("============================================================\n")
cat("Output :", RESULT$output_root, "\n")
cat("Elapsed:", round(elapsed_seconds / 60, 2), "minutes\n")
cat("Core   :", QMB_ACT_METHOD_ID, "|", QMB_ACT_VERSION, "\n\n")

cat("Selected cutoff distribution:\n")
print(table(RESULT$selected_cutoff$SELECTED_CUTOFF))

cat("\nStage-2 six-arm comparison — PAPER-READY:\n")
print(RESULT$paper_stage2_6arm_comparison, row.names = FALSE)

if (isTRUE(RUN_STAGE3_BOOTSTRAP)) {
  cat("\nStage-3 six-arm all-B comparison — PRIMARY PAPER ANALYSIS:\n")
  print(RESULT$paper_stage3_6arm_all_b_comparison, row.names = FALSE)
}

cat("\nKey paper-facing outputs are under:\n")
cat("  ", OUTPUT_ROOT, "\n", sep = "")
cat("\nPaper-ready comparison files:\n")
cat("  PAPER_STAGE2_6ARM_COMPARISON.csv\n")
if (isTRUE(RUN_STAGE3_BOOTSTRAP)) {
  cat("  PAPER_STAGE3_6ARM_ALL_B_COMPARISON.csv\n")
}
cat("  PAPER_RAW_MONTHLY_METRICS_BY_COMBINATION.csv\n")
cat("============================================================\n")
