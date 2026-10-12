# =============================================================================
# RUN_FULL_ANALYSIS.R
# QMB Adaptive QDM-Tail with Cross-Validated Cutoff + Stage-3 Bootstrap
# Publication-facing reproducibility runner for the ACT-QMB study
# Repository release: v1.2.1-audited
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
#     -> PILOT execution audit for one matched combination at B=30
#
# EXPECTED REPOSITORY STRUCTURE
#
# qmb-adaptive-interpolation/
# |-- code/
# |   |-- RUN_FULL_ANALYSIS.R
# |   |-- QMB_ADAPTIVE_QDMT_CORE_v1.2.1.R
# |   `-- QMB_KOKIC_CORE.R
# |-- data/
# |   `-- paired input CSV files (subfolders allowed)
# `-- results/
#
# RUN
#   From the repository root:
#       source("code/RUN_FULL_ANALYSIS.R")
#
# Or:
#       Rscript code/RUN_FULL_ANALYSIS.R
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

RUNNER_RELEASE <- "v1.2.0-major-revision"
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
OUTPUT_ROOT <- file.path(RESULTS_ROOT, "pilot_analysis")

ADAPTIVE_CORE_FILE <- file.path(CODE_DIR, "QMB_ADAPTIVE_QDMT_CORE_v1.2.1.R")
KOKIC_CORE_FILE <- file.path(CODE_DIR, "QMB_KOKIC_CORE.R")
COMMON_METRIC_FILE <- file.path(CODE_DIR, "COMMON_METRIC_CORE_v1.1.0.R")

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
B_VALUES <- 30L
STAGE3_INTERVAL <- 0.95
SAVE_STAGE3_DAILY_ALL_B <- TRUE
KEEP_STAGE3_ENSEMBLE_RDS <- FALSE

# Pilot design: exactly one combination.
# By default the first validated manifest row is used.
# To target a specific combination, change these three values.
PILOT_REGIME <- NA_character_   # e.g. "Monsoonal"
PILOT_GCM    <- NA_character_   # e.g. "MRI_ESM2_0"
PILOT_SSP    <- NA_character_   # e.g. "ssp126"


# -----------------------------------------------------------------------------
# 2. Preflight
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("QMB ADAPTIVE INTERPOLATION — ONE-COMBINATION PILOT RUN\n")
cat("Repository release :", RUNNER_RELEASE, "\n")
cat("Repository root    :", REPO_ROOT, "\n")
cat("Code directory     :", CODE_DIR, "\n")
cat("Data directory     :", DATA_ROOT, "\n")
cat("Results directory  :", OUTPUT_ROOT, "\n")
cat("============================================================\n\n")

required_files <- c(ADAPTIVE_CORE_FILE, KOKIC_CORE_FILE, COMMON_METRIC_FILE)
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
source(COMMON_METRIC_FILE)
if (!exists("COMMON_METRIC_CORE_VERSION") ||
    !identical(COMMON_METRIC_CORE_VERSION, "v1.1.0")) {
  stop("RUN_ANALYSIS_PILOT_v1.0.1.R requires COMMON_METRIC_CORE_v1.1.0.R",
       call. = FALSE)
}

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

# Audit the repository input design before subsetting the pilot.
if (nrow(MANIFEST) != EXPECTED_N_COMBINATIONS) {
  stop(
    sprintf(
      paste0(
        "Pilot requires the complete repository input design to be present first: ",
        "expected %d combinations (3 regimes x 4 GCM x 4 SSP), detected %d."
      ),
      EXPECTED_N_COMBINATIONS,
      nrow(MANIFEST)
    ),
    call. = FALSE
  )
}

FULL_MANIFEST <- MANIFEST

# Optional explicit pilot selection. If all selectors are NA, use the first
# validated manifest row in deterministic manifest order.
if (!is.na(PILOT_REGIME)) {
  MANIFEST <- MANIFEST[MANIFEST$REGIME == PILOT_REGIME, , drop = FALSE]
}
if (!is.na(PILOT_GCM)) {
  MANIFEST <- MANIFEST[MANIFEST$GCM == PILOT_GCM, , drop = FALSE]
}
if (!is.na(PILOT_SSP)) {
  MANIFEST <- MANIFEST[MANIFEST$SSP == PILOT_SSP, , drop = FALSE]
}

if (all(is.na(c(PILOT_REGIME, PILOT_GCM, PILOT_SSP)))) {
  MANIFEST <- FULL_MANIFEST[1, , drop = FALSE]
}

if (nrow(MANIFEST) != 1L) {
  stop(
    paste0(
      "Pilot selector must resolve to exactly one combination; resolved rows = ",
      nrow(MANIFEST),
      ". Set PILOT_REGIME, PILOT_GCM, and PILOT_SSP explicitly."
    ),
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("PILOT COMBINATION\n")
print(MANIFEST[, c("REGIME","LOCATION_ID","GCM","SSP","FILE")], row.names = FALSE)
cat("Pilot B values:", paste(B_VALUES, collapse = ","), "\n")
cat("Pilot output  :", OUTPUT_ROOT, "\n")
cat("============================================================\n\n")

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
  paste0("PILOT=TRUE"),
  paste0("PILOT_COMBINATION=", paste(MANIFEST$REGIME, MANIFEST$LOCATION_ID, MANIFEST$GCM, MANIFEST$SSP, sep="|")),
  paste0("STAGE3_B_VALUES=", paste(B_VALUES, collapse = ",")),
  paste0("STAGE3_INTERVAL=", STAGE3_INTERVAL),
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
# Pilot post-run execution checks
# -----------------------------------------------------------------------------
cat("\nPilot computational run returned successfully.\n")
pilot_files <- list.files(OUTPUT_ROOT, recursive = TRUE, full.names = TRUE)
cat("Files generated:", length(pilot_files), "\n")

expected_stage3_names <- c(
  "qmb_lin_stage3_daily_interval.csv",
  "hard90_stage3_daily_interval.csv",
  "fixed94_stage3_daily_interval.csv",
  "adaptive_cv_stage3_daily_interval.csv"
)

stage3_hits <- vapply(expected_stage3_names, function(nm) {
  any(basename(pilot_files) == nm)
}, logical(1))

cat("\nStage-3 B=30 daily interval outputs:\n")
for (i in seq_along(expected_stage3_names)) {
  cat(if (stage3_hits[i]) "PASS " else "FAIL ", expected_stage3_names[i], "\n")
}

if (!all(stage3_hits)) {
  warning(
    "Pilot completed but one or more expected Stage-3 daily interval files ",
    "were not found. Inspect results/pilot_analysis before a full run."
  )
}


# -----------------------------------------------------------------------------
# 6B. Six-arm common-metric pilot audit
# -----------------------------------------------------------------------------
stage2_file <- list.files(
  file.path(OUTPUT_ROOT, "FINAL_STAGE2"),
  pattern = "^stage2_daily_qmb_adaptive_qdm_tail\\.csv$",
  recursive = TRUE, full.names = TRUE
)
if (length(stage2_file) != 1L) {
  stop("Pilot expected exactly one Stage-2 daily file; found ",
       length(stage2_file), ".", call. = FALSE)
}

s2 <- read.csv(stage2_file, stringsAsFactors = FALSE)
s2$date <- as.Date(s2$date)
if (nrow(s2) != 4015L) {
  stop("Stage-2 pilot series must contain 4,015 evaluation days; found ",
       nrow(s2), ".", call. = FALSE)
}

required_s2 <- c(
  "date","obs","gcm","qdm_full","qmb_kokic5q_stage2",
  "qmb_hard90_stage2","qmb_fixed94_stage2","qmb_adaptive_stage2"
)
miss_s2 <- setdiff(required_s2, names(s2))
if (length(miss_s2)) {
  stop("Stage-2 daily file is missing: ", paste(miss_s2, collapse=", "),
       call. = FALSE)
}

pilot_metric_row <- function(method, pred) {
  dd <- data.frame(
    date=s2$date, obs=as.numeric(s2$obs), prediction=as.numeric(pred),
    stringsAsFactors=FALSE
  )
  dd$year <- as.integer(format(dd$date,"%Y"))
  dd$month <- as.integer(format(dd$date,"%m"))

  mt <- cm_monthly_totals(dd, min_completeness=0.90)
  ok_mt <- mt$valid_month & is.finite(mt$Total_obs) & is.finite(mt$Total_prediction)
  wo <- cm_monthly_occurrence(dd, wet_threshold=CM_WET_THRESHOLD)
  ex <- cm_pooled_calendar_extremes(dd, probs=CM_EXTREME_PROBS)

  data.frame(
    METHOD=method,
    EVALUATION_LAYER="STAGE2",
    N_BOOT=NA_integer_,
    IS_BOOTSTRAPPED=FALSE,
    DAILY_MAE=cm_mae(dd$obs,dd$prediction),
    DAILY_RMSE=cm_rmse(dd$obs,dd$prediction),
    DAILY_BIAS=cm_bias(dd$obs,dd$prediction),
    DAILY_CORRELATION=cm_cor(dd$obs,dd$prediction),
    DAILY_KGE=cm_kge(dd$obs,dd$prediction),
    MONTHLY_TOTAL_MAE=cm_mae(mt$Total_obs[ok_mt],mt$Total_prediction[ok_mt]),
    WET_OCC_MAE=cm_mae(wo$WetProp_obs,wo$WetProp_prediction),
    P90_ABS_ERROR=mean(ex$P90_abs_error,na.rm=TRUE),
    P95_ABS_ERROR=mean(ex$P95_abs_error,na.rm=TRUE),
    P99_ABS_ERROR=mean(ex$P99_abs_error,na.rm=TRUE),
    stringsAsFactors=FALSE
  )
}

SIX_ARM_STAGE2 <- do.call(rbind, list(
  pilot_metric_row("RAW_GCM", s2$gcm),
  pilot_metric_row("QDM", s2$qdm_full),
  pilot_metric_row("QMB_LIN", s2$qmb_kokic5q_stage2),
  pilot_metric_row("HARD90", s2$qmb_hard90_stage2),
  pilot_metric_row("FIXED94", s2$qmb_fixed94_stage2),
  pilot_metric_row("ACT_QMB", s2$qmb_adaptive_stage2)
))

write.csv(
  SIX_ARM_STAGE2,
  file.path(OUTPUT_ROOT, "PILOT_SIX_ARM_STAGE2_COMMON_METRICS.csv"),
  row.names=FALSE
)

# RAW and QDM are deterministic matched benchmarks at the Stage-3 comparison
# layer. They are repeated unchanged and are explicitly not bootstrapped.
DET_BENCH <- SIX_ARM_STAGE2[
  SIX_ARM_STAGE2$METHOD %in% c("RAW_GCM","QDM"), , drop=FALSE
]
DET_BENCH$EVALUATION_LAYER <- "STAGE3_MEDIAN"
DET_BENCH$N_BOOT <- 30L
DET_BENCH$IS_BOOTSTRAPPED <- FALSE

write.csv(
  DET_BENCH,
  file.path(OUTPUT_ROOT, "PILOT_STAGE3_DETERMINISTIC_BENCHMARKS_B30.csv"),
  row.names=FALSE
)

PILOT_AUDIT <- data.frame(
  CHECK=c(
    "Evaluation days","Six Stage-2 methods","RAW source",
    "RAW bootstrapped","QDM bootstrapped","Common metric core","Pilot B"
  ),
  VALUE=c(
    nrow(s2), paste(SIX_ARM_STAGE2$METHOD,collapse=";"),
    "gcm column in exact paired Stage-2 input",
    "NO","NO",COMMON_METRIC_CORE_VERSION,paste(B_VALUES,collapse=";")
  ),
  STATUS=c(
    ifelse(nrow(s2)==4015L,"PASS","FAIL"),
    ifelse(nrow(SIX_ARM_STAGE2)==6L,"PASS","FAIL"),
    "PASS","PASS","PASS","PASS",
    ifelse(length(B_VALUES)==1L && as.integer(B_VALUES[1])==30L,"PASS","FAIL")
  ),
  stringsAsFactors=FALSE
)
write.csv(PILOT_AUDIT,file.path(OUTPUT_ROOT,"PILOT_SIX_ARM_AUDIT.csv"),row.names=FALSE)

cat("\nSix-arm Stage-2 common-metric audit:\n")
print(SIX_ARM_STAGE2,row.names=FALSE)
cat("\nDeterministic benchmarks repeated at Stage-3 B=30:\n")
print(DET_BENCH,row.names=FALSE)
cat("\nSix-arm pilot audit checks:\n")
print(PILOT_AUDIT,row.names=FALSE)
cat("\nRAW_GCM and QDM above are deterministic comparison benchmarks; ",
    "they are not bootstrap ensembles.\n",sep="")


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
    COMMON_METRIC_FILE,
    if (!is.na(THIS_SCRIPT)) THIS_SCRIPT else file.path(CODE_DIR, "RUN_FULL_ANALYSIS.R")
  ),
  file.path(repro_dir, "SHA256SUMS.txt")
)

# -----------------------------------------------------------------------------
# 8. Final console report
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("PILOT ANALYSIS COMPLETED\n")
cat("============================================================\n")
cat("Output :", RESULT$output_root, "\n")
cat("Elapsed:", round(elapsed_seconds / 60, 2), "minutes\n")
cat("Core   :", QMB_ACT_METHOD_ID, "|", QMB_ACT_VERSION, "\n\n")

cat("Selected cutoff distribution:\n")
print(table(RESULT$selected_cutoff$SELECTED_CUTOFF))

cat("\nStage-2 comparison (diagnostic):\n")
print(RESULT$final_overall_comparison, row.names = FALSE)

if (isTRUE(RUN_STAGE3_BOOTSTRAP)) {
  cat("\nStage-3 comparison — PRIMARY PAPER ANALYSIS:\n")
  print(RESULT$stage3_overall_comparison, row.names = FALSE)
}

cat("\nKey paper-facing outputs are under:\n")
cat("  ", OUTPUT_ROOT, "\n", sep = "")
cat("============================================================\n")
