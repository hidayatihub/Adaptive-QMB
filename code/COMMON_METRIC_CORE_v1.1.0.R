# =============================================================================
# COMMON_METRIC_CORE_v1.1.0.R
# =============================================================================
# Locked common evaluation layer for the six-arm major-revision design:
#   RAW GCM      (deterministic benchmark)
#   QDM          (deterministic benchmark)
#   QMB-LIN      (Stage 2 + Stage 3 median)
#   HARD90       (Stage 2 + Stage 3 median)
#   FIXED94      (Stage 2 + Stage 3 median)
#   ACT-QMB      (Stage 2 + Stage 3 median)
#
# RAW GCM and QDM are repeated unchanged at each Stage-3 evaluation layer only
# to provide matched benchmarks. They are NEVER treated as bootstrapped methods.
#
# Major-revision locked definitions:
#   evaluation years : 2015-2025
#   leap day         : removed
#   wet day          : rainfall > 0.1 mm/day
#   quantile type    : R type = 7
#   tail pooling     : pooled calendar month across evaluation years
#   tail days        : ALL days, including zero/dry days
#   nominal interval : 0.95 unless explicitly overridden
#
# The purpose is NOT to change any correction algorithm.
# The purpose is to re-evaluate all daily series using one identical metric
# implementation and one identical aggregation protocol.
#
# Locked evaluation protocol
# --------------------------
# Target period       : 2015-2025
# Leap day            : removed
# Wet day             : rainfall > 0.1 mm/day
#
# Daily metrics:
#   MAE, RMSE, Bias, Pearson r, KGE
#
# Occurrence metrics:
#   monthly wet-day proportion -> MAE, RMSE, Bias across year-months
#
# Monthly-total metrics:
#   year-month rainfall totals -> MAE, RMSE, Bias, Pearson r, KGE
#
# Extreme metrics:
#   P90, P95, P99 are calculated with the locked definition:
#   "pooled_calendar_month_across_target_years"
#
#   Example:
#     all January daily values from 2015-2025 are pooled -> January P95
#     ...
#     all December daily values from 2015-2025 are pooled -> December P95
#
#   MAE/RMSE/Bias are then calculated across the 12 calendar-month percentiles.
#
# Interval metrics (only if lower/upper exist):
#   Coverage, coverage error, mean/median interval width,
#   mean/median interval score.
#
# Canonical input accepted by common_metric_evaluate():
#   REGIME, LOCATION_ID, GCM, SSP, METHOD,
#   EVALUATION_LAYER (preferred) or STAGE (legacy alias),
#   N_BOOT, IS_BOOTSTRAPPED,
#   date, obs, prediction, lower, upper
#
# EVALUATION_LAYER:
#   STAGE2
#   STAGE3_MEDIAN
#
# N_BOOT is NA for Stage 2. For matched Stage-3 comparison it is one of:
#   30, 100, 250, 500, 1000, 2000
#
# IS_BOOTSTRAPPED:
#   FALSE for RAW_GCM and QDM
#   TRUE  for QMB_LIN, HARD90, FIXED94, ACT_QMB at Stage 3
#
# =============================================================================


# -----------------------------------------------------------------------------
# v1.1.0 locked design constants
# -----------------------------------------------------------------------------
COMMON_METRIC_CORE_VERSION <- "v1.1.0"
CM_TARGET_YEARS <- c(2015L, 2025L)
CM_WET_THRESHOLD <- 0.1
CM_EXTREME_PROBS <- c(0.90, 0.95, 0.99)
CM_QUANTILE_TYPE <- 7L
CM_B_VALUES <- c(30L, 100L, 250L, 500L, 1000L, 2000L)
CM_METHODS <- c("RAW_GCM", "QDM", "QMB_LIN", "HARD90", "FIXED94", "ACT_QMB")
CM_BOOTSTRAPPED_METHODS <- c("QMB_LIN", "HARD90", "FIXED94", "ACT_QMB")
CM_DETERMINISTIC_BENCHMARKS <- c("RAW_GCM", "QDM")

cm_prepare_major_revision_schema <- function(df) {
  x <- df
  if (!"EVALUATION_LAYER" %in% names(x) && "STAGE" %in% names(x)) {
    x$EVALUATION_LAYER <- as.character(x$STAGE)
  }
  if (!"STAGE" %in% names(x) && "EVALUATION_LAYER" %in% names(x)) {
    x$STAGE <- as.character(x$EVALUATION_LAYER)
  }
  if (!"N_BOOT" %in% names(x)) x$N_BOOT <- NA_integer_
  if (!"IS_BOOTSTRAPPED" %in% names(x)) {
    x$IS_BOOTSTRAPPED <- x$METHOD %in% CM_BOOTSTRAPPED_METHODS &
      toupper(as.character(x$EVALUATION_LAYER)) == "STAGE3_MEDIAN"
  }
  x$METHOD <- toupper(gsub("[^A-Z0-9]+", "_", as.character(x$METHOD)))
  x$METHOD[x$METHOD %in% c("QMB", "QMB_KOKIC", "QMB_KOKIC5Q")] <- "QMB_LIN"
  x$METHOD[x$METHOD %in% c("ADAPTIVE_CV", "ACT", "ACTQMB")] <- "ACT_QMB"
  x$METHOD[x$METHOD %in% c("RAW", "GCM", "RAWGCM")] <- "RAW_GCM"
  x
}

cm_assert_benchmark_invariance <- function(canonical_daily, tolerance = 1e-12) {
  d <- cm_prepare_major_revision_schema(canonical_daily)
  d$date <- cm_parse_date(d$date)
  for (meth in CM_DETERMINISTIC_BENCHMARKS) {
    z <- d[d$METHOD == meth, , drop = FALSE]
    if (!nrow(z)) next
    key <- interaction(z$REGIME, z$LOCATION_ID, z$GCM, z$SSP, z$date, drop = TRUE)
    sp <- split(seq_len(nrow(z)), key)
    bad <- vapply(sp, function(ix) {
      v <- z$prediction[ix]
      v <- v[is.finite(v)]
      length(v) > 1L && (max(v) - min(v) > tolerance)
    }, logical(1))
    if (any(bad)) {
      stop(meth, " changed across evaluation layers/B. Deterministic benchmark invariance failed.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# General utilities
# -----------------------------------------------------------------------------

cm_norm <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z <- gsub("[^A-Z0-9]+", "_", z)
  z <- gsub("_+", "_", z)
  gsub("^_|_$", "", z)
}

cm_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

cm_parse_date <- function(x) {
  if (inherits(x, "Date")) return(x)

  if (inherits(x, "POSIXt")) {
    return(as.Date(x))
  }

  z <- trimws(as.character(x))
  z[z %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_

  out <- as.Date(rep(NA_character_, length(z)))

  # Numeric dates need special treatment because several Stage-3 outputs
  # store R Date values as their underlying integer representation.
  num <- suppressWarnings(as.numeric(z))

  # R Date serial: days since 1970-01-01.
  # Values around 16,000 correspond to 2013-2016.
  idx_rdate <- is.na(out) & !is.na(num) & num >= 10000 & num < 30000
  if (any(idx_rdate)) {
    out[idx_rdate] <- as.Date(num[idx_rdate], origin = "1970-01-01")
  }

  # Excel serial date.
  idx_excel <- is.na(out) & !is.na(num) & num >= 30000 & num < 80000
  if (any(idx_excel)) {
    out[idx_excel] <- as.Date(num[idx_excel], origin = "1899-12-30")
  }

  # ISO date or ISO datetime.
  idx <- is.na(out) & !is.na(z) &
    grepl("^\\d{4}-\\d{1,2}-\\d{1,2}", z)

  if (any(idx)) {
    first10 <- substr(z[idx], 1, 10)
    out[idx] <- suppressWarnings(
      as.Date(first10, format = "%Y-%m-%d")
    )
  }

  fmts <- c(
    "%Y/%m/%d",
    "%Y.%m.%d",
    "%d/%m/%Y",
    "%m/%d/%Y",
    "%d-%m-%Y",
    "%m-%d-%Y",
    "%d.%m.%Y",
    "%m.%d.%Y",
    "%Y%m%d"
  )

  for (fmt in fmts) {
    idx <- is.na(out) & !is.na(z)
    if (!any(idx)) break

    parsed <- suppressWarnings(
      as.Date(z[idx], format = fmt)
    )
    out[idx] <- parsed
  }

  out
}

cm_mean <- function(x) {
  x <- cm_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

cm_sd <- function(x) {
  x <- cm_num(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

cm_mae <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (!any(ok)) return(NA_real_)
  mean(abs(pred[ok] - obs[ok]))
}

cm_rmse <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((pred[ok] - obs[ok])^2))
}

cm_bias <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (!any(ok)) return(NA_real_)
  mean(pred[ok] - obs[ok])
}

cm_cor <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 3L) return(NA_real_)

  o <- obs[ok]
  p <- pred[ok]

  if (stats::sd(o) <= 0 || stats::sd(p) <= 0) return(NA_real_)

  suppressWarnings(stats::cor(o, p, method = "pearson"))
}

cm_kge <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 3L) return(NA_real_)

  o <- obs[ok]
  p <- pred[ok]

  r <- cm_cor(o, p)
  if (!is.finite(r)) return(NA_real_)

  mo <- mean(o)
  mp <- mean(p)
  sdo <- stats::sd(o)
  sdp <- stats::sd(p)

  if (!is.finite(mo) || !is.finite(mp) ||
      !is.finite(sdo) || !is.finite(sdp) ||
      mo == 0 || sdo == 0) {
    return(NA_real_)
  }

  alpha <- sdp / sdo
  beta <- mp / mo

  1 - sqrt(
    (r - 1)^2 +
      (alpha - 1)^2 +
      (beta - 1)^2
  )
}

cm_quantile <- function(x, probs) {
  x <- cm_num(x)
  x <- x[is.finite(x)]

  if (!length(x)) {
    return(rep(NA_real_, length(probs)))
  }

  as.numeric(stats::quantile(
    x,
    probs = probs,
    type = CM_QUANTILE_TYPE,
    names = FALSE,
    na.rm = TRUE
  ))
}

cm_wet <- function(x, threshold = 0.1) {
  is.finite(x) & x > threshold
}

cm_days_in_month_noleap <- function(month) {
  c(
    31L, 28L, 31L, 30L,
    31L, 30L, 31L, 31L,
    30L, 31L, 30L, 31L
  )[as.integer(month)]
}

cm_rbind <- function(lst) {
  if (!length(lst)) return(data.frame())

  keep <- vapply(
    lst,
    function(x) is.data.frame(x) && nrow(x) > 0L,
    logical(1)
  )

  lst <- lst[keep]
  if (!length(lst)) return(data.frame())

  cols <- unique(unlist(lapply(lst, names), use.names = FALSE))

  lst <- lapply(lst, function(x) {
    miss <- setdiff(cols, names(x))
    for (nm in miss) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  })

  out <- do.call(rbind, lst)
  rownames(out) <- NULL
  out
}

cm_key <- function(location, gcm, ssp) {
  paste(location, gcm, ssp, sep = "||")
}

# -----------------------------------------------------------------------------
# Metadata parsing
# -----------------------------------------------------------------------------

cm_parse_regime <- function(x) {
  z <- cm_norm(x)

  if (grepl("MONSOONAL|MONSUN|JAWA_TIMUR|EAST_JAVA", z)) {
    return("Monsoonal")
  }

  if (grepl("EQUATORIAL|EKUATORIAL|SUMATERA_BARAT|WEST_SUMATRA", z)) {
    return("Equatorial")
  }

  if (grepl("LOCAL|LOKAL|MALUKU", z)) {
    return("Local")
  }

  NA_character_
}

cm_parse_gcm <- function(x) {
  z <- cm_norm(x)

  if (grepl("CESM2_WACCM|CESM2WACCM", z)) return("CESM2_WACCM")
  if (grepl("CMCC_ESM2|CMCCESM2", z)) return("CMCC_ESM2")
  if (grepl("MRI_ESM2_0|MRIESM2_0|MRIESM20", z)) return("MRI_ESM2_0")
  if (grepl("TAIESM1|TAI_ESM1", z)) return("TaiESM1")

  NA_character_
}

cm_parse_ssp <- function(x) {
  z <- cm_norm(x)

  if (grepl("SSP126|SSP1_2_6", z)) return("ssp126")
  if (grepl("SSP245|SSP2_4_5", z)) return("ssp245")
  if (grepl("SSP370|SSP3_7_0", z)) return("ssp370")
  if (grepl("SSP585|SSP5_8_5", z)) return("ssp585")

  NA_character_
}

cm_parse_location <- function(x, regime = NA_character_) {
  z <- cm_norm(x)

  patterns <- list(
    Monsoonal_1 = c("MONSOONAL_1", "MONSUN_1", "JAWA_TIMUR_1", "EAST_JAVA_1"),
    Monsoonal_2 = c("MONSOONAL_2", "MONSUN_2", "JAWA_TIMUR_2", "EAST_JAVA_2"),
    Equatorial_1 = c("EQUATORIAL_1", "EKUATORIAL_1", "SUMATERA_BARAT_1", "WEST_SUMATRA_1"),
    Equatorial_2 = c("EQUATORIAL_2", "EKUATORIAL_2", "SUMATERA_BARAT_2", "WEST_SUMATRA_2"),
    Local_1 = c("LOCAL_1", "LOKAL_1", "MALUKU_1"),
    Local_2 = c("LOCAL_2", "LOKAL_2", "MALUKU_2")
  )

  for (nm in names(patterns)) {
    hit <- vapply(
      patterns[[nm]],
      function(p) grepl(p, z, fixed = TRUE),
      logical(1)
    )

    if (any(hit)) return(nm)
  }

  if (!is.na(regime)) {
    if (grepl("LOCATION_1|LOC_1|LOKASI_1|GRID_1|SITE_1|POINT_1", z)) {
      return(paste0(regime, "_1"))
    }

    if (grepl("LOCATION_2|LOC_2|LOKASI_2|GRID_2|SITE_2|POINT_2", z)) {
      return(paste0(regime, "_2"))
    }
  }

  NA_character_
}

cm_location_from_routing <- function(regime, gcm) {
  if (is.na(regime) || is.na(gcm)) return(NA_character_)

  loc_no <- switch(
    gcm,
    "CESM2_WACCM" = 1L,
    "CMCC_ESM2" = 1L,
    "TaiESM1" = 1L,
    "MRI_ESM2_0" = 2L,
    NA_integer_
  )

  if (is.na(loc_no)) return(NA_character_)
  paste0(regime, "_", loc_no)
}

cm_metadata_from_path <- function(path, allow_routing = TRUE) {
  txt <- normalizePath(
    path,
    winslash = "/",
    mustWork = FALSE
  )

  regime <- cm_parse_regime(txt)
  gcm <- cm_parse_gcm(txt)
  ssp <- cm_parse_ssp(txt)
  location <- cm_parse_location(txt, regime)

  if (is.na(regime) && !is.na(location)) {
    regime <- sub("_[12]$", "", location)
  }

  if (is.na(location) &&
      isTRUE(allow_routing) &&
      !is.na(regime) &&
      !is.na(gcm)) {
    location <- cm_location_from_routing(regime, gcm)
  }

  data.frame(
    REGIME = regime,
    LOCATION_ID = location,
    GCM = gcm,
    SSP = ssp,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Canonicalization helpers
# -----------------------------------------------------------------------------

cm_first_col <- function(df, candidates) {
  n <- cm_norm(names(df))
  c <- cm_norm(candidates)

  hit <- which(n %in% c)
  if (!length(hit)) return(NA_character_)
  names(df)[hit[1]]
}

cm_canonicalize <- function(df,
                            metadata,
                            method,
                            stage,
                            n_boot = NA_integer_,
                            date_col,
                            obs_col,
                            pred_col,
                            lower_col = NA_character_,
                            upper_col = NA_character_,
                            source_file = NA_character_,
                            source_priority = 100L) {
  if (is.na(date_col) || is.na(obs_col) || is.na(pred_col)) {
    stop("date_col, obs_col, and pred_col must exist.", call. = FALSE)
  }

  out <- data.frame(
    REGIME = metadata$REGIME[1],
    LOCATION_ID = metadata$LOCATION_ID[1],
    GCM = metadata$GCM[1],
    SSP = metadata$SSP[1],
    METHOD = method,
    STAGE = stage,
    N_BOOT = n_boot,
    date = cm_parse_date(df[[date_col]]),
    obs = cm_num(df[[obs_col]]),
    prediction = cm_num(df[[pred_col]]),
    lower = if (!is.na(lower_col) && lower_col %in% names(df)) {
      cm_num(df[[lower_col]])
    } else {
      rep(NA_real_, nrow(df))
    },
    upper = if (!is.na(upper_col) && upper_col %in% names(df)) {
      cm_num(df[[upper_col]])
    } else {
      rep(NA_real_, nrow(df))
    },
    SOURCE_FILE = source_file,
    SOURCE_PRIORITY = as.integer(source_priority),
    stringsAsFactors = FALSE
  )

  out
}

cm_read_csv_safe <- function(path) {
  tryCatch(
    utils::read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) NULL
  )
}

# Dedicated RAW output created by RUN_CHRIPS_GCM_ALL_IN_ONE_v1.1.0.R
cm_adapter_raw_master <- function(path) {
  df <- cm_read_csv_safe(path)
  if (is.null(df) || !nrow(df)) return(data.frame())

  required <- c(
    "REGIME", "LOCATION_ID", "GCM", "SSP",
    "date", "obs", "prediction"
  )

  if (!all(required %in% names(df))) return(data.frame())

  if (!"METHOD" %in% names(df)) df$METHOD <- "RAW_GCM"
  if (!"STAGE" %in% names(df)) df$STAGE <- "RAW"
  if (!"N_BOOT" %in% names(df)) df$N_BOOT <- NA_integer_
  if (!"lower" %in% names(df)) df$lower <- NA_real_
  if (!"upper" %in% names(df)) df$upper <- NA_real_

  df$SOURCE_FILE <- path
  df$SOURCE_PRIORITY <- 200L

  df[, c(
    "REGIME", "LOCATION_ID", "GCM", "SSP",
    "METHOD", "STAGE", "N_BOOT",
    "date", "obs", "prediction", "lower", "upper",
    "SOURCE_FILE", "SOURCE_PRIORITY"
  )]
}

cm_adapter_qdm <- function(path) {
  df <- cm_read_csv_safe(path)
  if (is.null(df) || !nrow(df)) return(data.frame())

  metadata <- cm_metadata_from_path(path, allow_routing = TRUE)

  date_col <- cm_first_col(df, c("date", "tanggal"))
  obs_col <- cm_first_col(df, c("obs", "chirps", "chrips", "observed"))
  pred_col <- cm_first_col(df, c("qdm", "qdm_full", "corrected", "prediction"))

  if (any(is.na(c(date_col, obs_col, pred_col)))) return(data.frame())

  cm_canonicalize(
    df,
    metadata = metadata,
    method = "QDM",
    stage = "FINAL_DETERMINISTIC",
    date_col = date_col,
    obs_col = obs_col,
    pred_col = pred_col,
    source_file = path,
    source_priority = 200L
  )
}

cm_adapter_qmb_stage2 <- function(path) {
  df <- cm_read_csv_safe(path)
  if (is.null(df) || !nrow(df)) return(data.frame())

  metadata <- cm_metadata_from_path(path, allow_routing = TRUE)

  date_col <- cm_first_col(df, c("date", "tanggal"))
  obs_col <- cm_first_col(df, c("obs", "chirps", "chrips", "observed"))
  pred_col <- cm_first_col(df, c("qmb_stage2", "qmb_kokic5q_stage2"))

  if (any(is.na(c(date_col, obs_col, pred_col)))) return(data.frame())

  cm_canonicalize(
    df,
    metadata = metadata,
    method = "QMB_KOKIC",
    stage = "STAGE2",
    date_col = date_col,
    obs_col = obs_col,
    pred_col = pred_col,
    source_file = path,
    source_priority = 200L
  )
}

cm_adapter_qmb_stage3 <- function(path, n_boot = 2000L) {
  df <- cm_read_csv_safe(path)
  if (is.null(df) || !nrow(df)) return(data.frame())

  metadata <- cm_metadata_from_path(path, allow_routing = TRUE)

  date_col <- cm_first_col(df, c("date", "tanggal"))
  obs_col <- cm_first_col(df, c("obs", "chirps", "chrips", "observed"))
  pred_col <- cm_first_col(df, c("median", "prediction", "pred", "stage3_median"))
  lower_col <- cm_first_col(df, c("lower", "lower95", "lwr"))
  upper_col <- cm_first_col(df, c("upper", "upper95", "upr"))

  if (any(is.na(c(date_col, obs_col, pred_col)))) return(data.frame())

  cm_canonicalize(
    df,
    metadata = metadata,
    method = "QMB_KOKIC",
    stage = paste0("STAGE3_MEDIAN_B", n_boot),
    n_boot = n_boot,
    date_col = date_col,
    obs_col = obs_col,
    pred_col = pred_col,
    lower_col = lower_col,
    upper_col = upper_col,
    source_file = path,
    source_priority = 300L
  )
}

# ACT stage-2 file contains RAW, QDM, QMB-5Q Stage-2, and ACT Stage-2.
# It is used as a fallback source only; dedicated RAW/QDM/QMB files receive
# higher priority.
cm_adapter_act_stage2_multi <- function(path) {
  df <- cm_read_csv_safe(path)
  if (is.null(df) || !nrow(df)) return(list())

  # Prefer metadata columns inside ACT file because they are explicit.
  metadata_path <- cm_metadata_from_path(path, allow_routing = TRUE)

  regime_col <- cm_first_col(df, c("REGIME"))
  location_col <- cm_first_col(df, c("LOCATION_ID"))
  gcm_meta_col <- cm_first_col(df, c("GCM"))
  ssp_col <- cm_first_col(df, c("SSP"))

  metadata <- metadata_path

  if (!is.na(regime_col)) {
    v <- unique(as.character(df[[regime_col]]))
    v <- v[!is.na(v)]
    if (length(v)) metadata$REGIME <- cm_parse_regime(v[1])
  }

  if (!is.na(location_col)) {
    v <- unique(as.character(df[[location_col]]))
    v <- v[!is.na(v)]
    if (length(v)) {
      loc_txt <- gsub("\\s+", "_", trimws(v[1]))
      loc <- cm_parse_location(loc_txt, metadata$REGIME)
      if (!is.na(loc)) metadata$LOCATION_ID <- loc
    }
  }

  if (!is.na(gcm_meta_col)) {
    v <- unique(as.character(df[[gcm_meta_col]]))
    v <- v[!is.na(v)]
    if (length(v)) metadata$GCM <- cm_parse_gcm(v[1])
  }

  if (!is.na(ssp_col)) {
    v <- unique(as.character(df[[ssp_col]]))
    v <- v[!is.na(v)]
    if (length(v)) metadata$SSP <- cm_parse_ssp(v[1])
  }

  if (is.na(metadata$LOCATION_ID) &&
      !is.na(metadata$REGIME) &&
      !is.na(metadata$GCM)) {
    metadata$LOCATION_ID <- cm_location_from_routing(
      metadata$REGIME,
      metadata$GCM
    )
  }

  date_col <- cm_first_col(df, c("date"))
  obs_col <- cm_first_col(df, c("obs", "chirps", "chrips"))

  if (any(is.na(c(date_col, obs_col)))) return(list())

  out <- list()

  raw_col <- cm_first_col(df, c("gcm"))
  if (!is.na(raw_col)) {
    out$RAW_GCM <- cm_canonicalize(
      df,
      metadata,
      method = "RAW_GCM",
      stage = "RAW",
      date_col = date_col,
      obs_col = obs_col,
      pred_col = raw_col,
      source_file = path,
      source_priority = 50L
    )
  }

  qdm_col <- cm_first_col(df, c("qdm_full", "qdm"))
  if (!is.na(qdm_col)) {
    out$QDM <- cm_canonicalize(
      df,
      metadata,
      method = "QDM",
      stage = "FINAL_DETERMINISTIC",
      date_col = date_col,
      obs_col = obs_col,
      pred_col = qdm_col,
      source_file = path,
      source_priority = 50L
    )
  }

  qmb_col <- cm_first_col(df, c("qmb_kokic5q_stage2", "qmb_stage2"))
  if (!is.na(qmb_col)) {
    out$QMB_STAGE2 <- cm_canonicalize(
      df,
      metadata,
      method = "QMB_KOKIC",
      stage = "STAGE2",
      date_col = date_col,
      obs_col = obs_col,
      pred_col = qmb_col,
      source_file = path,
      source_priority = 50L
    )
  }

  act_col <- cm_first_col(df, c("qmb_adaptive_stage2", "act_qmb_stage2"))
  if (!is.na(act_col)) {
    out$ACT_STAGE2 <- cm_canonicalize(
      df,
      metadata,
      method = "ACT_QMB",
      stage = "STAGE2",
      date_col = date_col,
      obs_col = obs_col,
      pred_col = act_col,
      source_file = path,
      source_priority = 250L
    )
  }

  out
}

# Generic ACT Stage-3 adapter.
# This supports files containing date/obs/median/lower/upper and an ACT/adaptive
# path token. The comparison runner is responsible for filtering ACT paths.
cm_adapter_act_stage3 <- function(path, n_boot = 2000L) {
  df <- cm_read_csv_safe(path)
  if (is.null(df) || !nrow(df)) return(data.frame())

  metadata <- cm_metadata_from_path(path, allow_routing = TRUE)

  # Allow metadata columns to override path.
  for (pair in list(
    c("REGIME", "REGIME"),
    c("LOCATION_ID", "LOCATION_ID"),
    c("GCM", "GCM"),
    c("SSP", "SSP")
  )) {
    col <- cm_first_col(df, pair[1])
    if (!is.na(col)) {
      v <- unique(as.character(df[[col]]))
      v <- v[!is.na(v)]
      if (length(v)) {
        if (pair[2] == "REGIME") metadata$REGIME <- cm_parse_regime(v[1])
        if (pair[2] == "LOCATION_ID") {
          loc <- cm_parse_location(
            gsub("\\s+", "_", v[1]),
            metadata$REGIME
          )
          if (!is.na(loc)) metadata$LOCATION_ID <- loc
        }
        if (pair[2] == "GCM") metadata$GCM <- cm_parse_gcm(v[1])
        if (pair[2] == "SSP") metadata$SSP <- cm_parse_ssp(v[1])
      }
    }
  }

  if (is.na(metadata$LOCATION_ID) &&
      !is.na(metadata$REGIME) &&
      !is.na(metadata$GCM)) {
    metadata$LOCATION_ID <- cm_location_from_routing(
      metadata$REGIME,
      metadata$GCM
    )
  }

  date_col <- cm_first_col(df, c("date", "tanggal"))
  obs_col <- cm_first_col(df, c("obs", "chirps", "chrips", "observed"))
  pred_col <- cm_first_col(
    df,
    c(
      "median",
      "act_median",
      "adaptive_median",
      "prediction",
      "stage3_median",
      "qmb_adaptive_stage3_median"
    )
  )
  lower_col <- cm_first_col(df, c("lower", "lower95", "lwr"))
  upper_col <- cm_first_col(df, c("upper", "upper95", "upr"))

  if (any(is.na(c(date_col, obs_col, pred_col)))) return(data.frame())

  cm_canonicalize(
    df,
    metadata,
    method = "ACT_QMB",
    stage = paste0("STAGE3_MEDIAN_B", n_boot),
    n_boot = n_boot,
    date_col = date_col,
    obs_col = obs_col,
    pred_col = pred_col,
    lower_col = lower_col,
    upper_col = upper_col,
    source_file = path,
    source_priority = 350L
  )
}

# -----------------------------------------------------------------------------
# Daily-series standardization
# -----------------------------------------------------------------------------

cm_standardize_daily <- function(df,
                                 target_years = CM_TARGET_YEARS,
                                 no_leap = TRUE) {
  df <- cm_prepare_major_revision_schema(df)
  required <- c(
    "REGIME", "LOCATION_ID", "GCM", "SSP",
    "METHOD", "STAGE", "EVALUATION_LAYER", "N_BOOT", "IS_BOOTSTRAPPED",
    "date", "obs", "prediction", "lower", "upper"
  )

  miss <- setdiff(required, names(df))
  if (length(miss)) {
    stop(
      paste0(
        "Canonical daily table is missing: ",
        paste(miss, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  d <- df

  original_date <- d$date
  d$date <- cm_parse_date(d$date)

  n_bad_date <- sum(is.na(d$date) & !is.na(original_date))
  if (n_bad_date > 0L) {
    bad_examples <- unique(as.character(original_date[is.na(d$date) & !is.na(original_date)]))
    bad_examples <- utils::head(bad_examples, 10L)
    warning(
      paste0(
        "Date parser could not parse ", n_bad_date,
        " row(s). Examples: ",
        paste(bad_examples, collapse = " | ")
      ),
      call. = FALSE
    )
  }
  d$obs <- cm_num(d$obs)
  d$prediction <- cm_num(d$prediction)
  d$lower <- cm_num(d$lower)
  d$upper <- cm_num(d$upper)

  d <- d[!is.na(d$date), , drop = FALSE]
  d$year <- as.integer(format(d$date, "%Y"))
  d$month <- as.integer(format(d$date, "%m"))
  d$day <- as.integer(format(d$date, "%d"))

  d <- d[
    d$year >= target_years[1] &
      d$year <= target_years[2],
    ,
    drop = FALSE
  ]

  if (isTRUE(no_leap)) {
    d <- d[!(d$month == 2L & d$day == 29L), , drop = FALSE]
  }

  d <- d[order(d$date), , drop = FALSE]
  d <- d[!duplicated(d$date), , drop = FALSE]
  rownames(d) <- NULL
  d
}

# -----------------------------------------------------------------------------
# Common metric components
# -----------------------------------------------------------------------------

cm_monthly_occurrence <- function(d, wet_threshold = 0.1) {
  key <- interaction(
    d$year,
    d$month,
    drop = TRUE,
    lex.order = TRUE
  )

  spl <- split(seq_len(nrow(d)), key)

  cm_rbind(lapply(spl, function(idx) {
    g <- d[idx, , drop = FALSE]
    ok <- is.finite(g$obs) & is.finite(g$prediction)

    o <- g$obs[ok]
    p <- g$prediction[ok]

    data.frame(
      year = g$year[1],
      month = g$month[1],
      n_days = sum(ok),
      WetProp_obs = if (sum(ok)) {
        mean(cm_wet(o, wet_threshold))
      } else {
        NA_real_
      },
      WetProp_prediction = if (sum(ok)) {
        mean(cm_wet(p, wet_threshold))
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
}

cm_monthly_totals <- function(d, min_completeness = 0.90) {
  key <- interaction(
    d$year,
    d$month,
    drop = TRUE,
    lex.order = TRUE
  )

  spl <- split(seq_len(nrow(d)), key)

  out <- cm_rbind(lapply(spl, function(idx) {
    g <- d[idx, , drop = FALSE]
    ok <- is.finite(g$obs) & is.finite(g$prediction)

    expected <- cm_days_in_month_noleap(g$month[1])
    n_pair <- sum(ok)
    completeness <- n_pair / expected

    data.frame(
      year = g$year[1],
      month = g$month[1],
      expected_days = expected,
      paired_days = n_pair,
      completeness = completeness,
      valid_month = completeness >= min_completeness,
      Total_obs = if (n_pair) sum(g$obs[ok]) else NA_real_,
      Total_prediction = if (n_pair) sum(g$prediction[ok]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  out
}

cm_pooled_calendar_extremes <- function(d,
                                        probs = c(0.90, 0.95, 0.99)) {
  rows <- vector("list", 12L)

  for (m in 1:12) {
    g <- d[
      d$month == m &
        is.finite(d$obs) &
        is.finite(d$prediction),
      ,
      drop = FALSE
    ]

    qo <- cm_quantile(g$obs, probs)
    qp <- cm_quantile(g$prediction, probs)

    row <- data.frame(
      month = m,
      month_name = month.abb[m],
      n_days = nrow(g),
      stringsAsFactors = FALSE
    )

    for (j in seq_along(probs)) {
      qn <- paste0("P", round(100 * probs[j]))
      row[[paste0(qn, "_obs")]] <- qo[j]
      row[[paste0(qn, "_prediction")]] <- qp[j]
      row[[paste0(qn, "_error")]] <- qp[j] - qo[j]
      row[[paste0(qn, "_abs_error")]] <- abs(qp[j] - qo[j])
    }

    rows[[m]] <- row
  }

  cm_rbind(rows)
}

cm_interval_score <- function(obs, lower, upper, nominal = 0.95) {
  alpha <- 1 - nominal
  ok <- is.finite(obs) & is.finite(lower) & is.finite(upper)

  if (!any(ok)) return(rep(NA_real_, length(obs)))

  y <- obs
  l <- lower
  u <- upper

  score <- rep(NA_real_, length(y))
  score[ok] <- (
    (u[ok] - l[ok]) +
      (2 / alpha) * (l[ok] - y[ok]) * (y[ok] < l[ok]) +
      (2 / alpha) * (y[ok] - u[ok]) * (y[ok] > u[ok])
  )

  score
}

# -----------------------------------------------------------------------------
# Main common evaluator for one canonical method-combination
# -----------------------------------------------------------------------------

common_metric_evaluate <- function(df,
                                   target_years = CM_TARGET_YEARS,
                                   wet_threshold = CM_WET_THRESHOLD,
                                   min_month_completeness = 0.90,
                                   extreme_probs = CM_EXTREME_PROBS,
                                   nominal_interval = 0.95) {
  df <- cm_prepare_major_revision_schema(df)
  d <- cm_standardize_daily(
    df,
    target_years = target_years,
    no_leap = TRUE
  )

  if (!nrow(d)) {
    stop(
      paste0(
        "No rows remained after date standardization for ",
        paste(unique(df$METHOD), collapse = ";"),
        " | ",
        paste(unique(df$LOCATION_ID), collapse = ";"),
        " | ",
        paste(unique(df$GCM), collapse = ";"),
        " | ",
        paste(unique(df$SSP), collapse = ";"),
        ". Check date encoding."
      ),
      call. = FALSE
    )
  }

  meta_cols <- c(
    "REGIME", "LOCATION_ID", "GCM", "SSP",
    "METHOD", "STAGE", "EVALUATION_LAYER", "N_BOOT", "IS_BOOTSTRAPPED"
  )

  meta <- d[1, meta_cols, drop = FALSE]

  ok <- is.finite(d$obs) & is.finite(d$prediction)
  o <- d$obs[ok]
  p <- d$prediction[ok]

  occurrence <- cm_monthly_occurrence(
    d,
    wet_threshold = wet_threshold
  )

  totals <- cm_monthly_totals(
    d,
    min_completeness = min_month_completeness
  )

  valid_tot <- totals$valid_month &
    is.finite(totals$Total_obs) &
    is.finite(totals$Total_prediction)

  extremes <- cm_pooled_calendar_extremes(
    d,
    probs = extreme_probs
  )

  result <- cbind(
    meta,
    data.frame(
      N_DAYS = sum(ok),
      EvaluationStartYear = target_years[1],
      EvaluationEndYear = target_years[2],
      WetThreshold = wet_threshold,

      Daily_MAE = cm_mae(o, p),
      Daily_RMSE = cm_rmse(o, p),
      Daily_Bias = cm_bias(o, p),
      Daily_r = cm_cor(o, p),
      Daily_KGE = cm_kge(o, p),

      Occurrence_NMonths = sum(
        is.finite(occurrence$WetProp_obs) &
          is.finite(occurrence$WetProp_prediction)
      ),
      WetProp_MAE = cm_mae(
        occurrence$WetProp_obs,
        occurrence$WetProp_prediction
      ),
      WetProp_RMSE = cm_rmse(
        occurrence$WetProp_obs,
        occurrence$WetProp_prediction
      ),
      WetProp_Bias = cm_bias(
        occurrence$WetProp_obs,
        occurrence$WetProp_prediction
      ),

      MonthlyTotal_NMonths = sum(valid_tot),
      MonthlyTotal_MAE = cm_mae(
        totals$Total_obs[valid_tot],
        totals$Total_prediction[valid_tot]
      ),
      MonthlyTotal_RMSE = cm_rmse(
        totals$Total_obs[valid_tot],
        totals$Total_prediction[valid_tot]
      ),
      MonthlyTotal_Bias = cm_bias(
        totals$Total_obs[valid_tot],
        totals$Total_prediction[valid_tot]
      ),
      MonthlyTotal_r = cm_cor(
        totals$Total_obs[valid_tot],
        totals$Total_prediction[valid_tot]
      ),
      MonthlyTotal_KGE = cm_kge(
        totals$Total_obs[valid_tot],
        totals$Total_prediction[valid_tot]
      ),

      ExtremeMetricDefinition = "pooled_calendar_month_across_target_years",
      ExtremeMetricNMonths = nrow(extremes),

      stringsAsFactors = FALSE
    )
  )

  for (prob in extreme_probs) {
    qn <- paste0("P", round(100 * prob))

    obs_col <- paste0(qn, "_obs")
    pred_col <- paste0(qn, "_prediction")

    result[[paste0(qn, "_MAE")]] <- cm_mae(
      extremes[[obs_col]],
      extremes[[pred_col]]
    )

    result[[paste0(qn, "_RMSE")]] <- cm_rmse(
      extremes[[obs_col]],
      extremes[[pred_col]]
    )

    result[[paste0(qn, "_Bias")]] <- cm_bias(
      extremes[[obs_col]],
      extremes[[pred_col]]
    )
  }

  # Interval metrics only if lower/upper are genuinely available.
  int_ok <- is.finite(d$obs) &
    is.finite(d$lower) &
    is.finite(d$upper)

  if (any(int_ok)) {
    coverage <- mean(
      d$obs[int_ok] >= d$lower[int_ok] &
        d$obs[int_ok] <= d$upper[int_ok]
    )

    width <- d$upper[int_ok] - d$lower[int_ok]

    score <- cm_interval_score(
      d$obs,
      d$lower,
      d$upper,
      nominal = nominal_interval
    )

    result$NominalInterval <- nominal_interval
    result$Coverage <- coverage
    result$CoverageErrorSigned <- coverage - nominal_interval
    result$CoverageErrorAbs <- abs(coverage - nominal_interval)
    result$MeanIntervalWidth <- mean(width)
    result$MedianIntervalWidth <- stats::median(width)
    result$MeanIntervalScore <- mean(score, na.rm = TRUE)
    result$MedianIntervalScore <- stats::median(score, na.rm = TRUE)
  } else {
    result$NominalInterval <- NA_real_
    result$Coverage <- NA_real_
    result$CoverageErrorSigned <- NA_real_
    result$CoverageErrorAbs <- NA_real_
    result$MeanIntervalWidth <- NA_real_
    result$MedianIntervalWidth <- NA_real_
    result$MeanIntervalScore <- NA_real_
    result$MedianIntervalScore <- NA_real_
  }

  list(
    metrics = result,
    occurrence = occurrence,
    monthly_totals = totals,
    pooled_extremes = extremes,
    standardized_daily = d
  )
}

# -----------------------------------------------------------------------------
# Group evaluator
# -----------------------------------------------------------------------------

common_metric_evaluate_all <- function(canonical_daily,
                                       target_years = c(2015L, 2025L),
                                       wet_threshold = 0.1,
                                       min_month_completeness = 0.90,
                                       extreme_probs = c(0.90, 0.95, 0.99),
                                       nominal_interval = 0.95) {
  if (!nrow(canonical_daily)) {
    return(list(
      metrics = data.frame(),
      occurrence = data.frame(),
      monthly_totals = data.frame(),
      pooled_extremes = data.frame()
    ))
  }

  group_cols <- c(
    "REGIME", "LOCATION_ID", "GCM", "SSP",
    "METHOD", "STAGE", "N_BOOT"
  )

  # Do not use interaction() directly here because N_BOOT is NA for
  # RAW/QDM/Stage-2 and NA interaction levels would otherwise be dropped.
  key_parts <- canonical_daily[, group_cols, drop = FALSE]
  key_parts$N_BOOT <- ifelse(
    is.na(key_parts$N_BOOT),
    "NA",
    as.character(key_parts$N_BOOT)
  )

  key <- do.call(
    paste,
    c(key_parts, sep = "||")
  )

  spl <- split(seq_len(nrow(canonical_daily)), key)

  met_list <- list()
  occ_list <- list()
  tot_list <- list()
  ext_list <- list()

  for (idx in spl) {
    d <- canonical_daily[idx, , drop = FALSE]

    res <- common_metric_evaluate(
      d,
      target_years = target_years,
      wet_threshold = wet_threshold,
      min_month_completeness = min_month_completeness,
      extreme_probs = extreme_probs,
      nominal_interval = nominal_interval
    )

    met <- res$metrics
    meta <- met[1, group_cols, drop = FALSE]

    occ <- cbind(
      meta[rep(1, nrow(res$occurrence)), , drop = FALSE],
      res$occurrence
    )

    tot <- cbind(
      meta[rep(1, nrow(res$monthly_totals)), , drop = FALSE],
      res$monthly_totals
    )

    ext <- cbind(
      meta[rep(1, nrow(res$pooled_extremes)), , drop = FALSE],
      res$pooled_extremes
    )

    met_list[[length(met_list) + 1L]] <- met
    occ_list[[length(occ_list) + 1L]] <- occ
    tot_list[[length(tot_list) + 1L]] <- tot
    ext_list[[length(ext_list) + 1L]] <- ext
  }

  list(
    metrics = cm_rbind(met_list),
    occurrence = cm_rbind(occ_list),
    monthly_totals = cm_rbind(tot_list),
    pooled_extremes = cm_rbind(ext_list)
  )
}

# -----------------------------------------------------------------------------
# Deduplicate source overlap by priority
# -----------------------------------------------------------------------------

cm_deduplicate_sources <- function(canonical_daily) {
  if (!nrow(canonical_daily)) return(canonical_daily)

  canonical_daily$COMBO_METHOD_STAGE <- paste(
    canonical_daily$LOCATION_ID,
    canonical_daily$GCM,
    canonical_daily$SSP,
    canonical_daily$METHOD,
    canonical_daily$STAGE,
    sep = "||"
  )

  # Select one source file per combination-method-stage using highest priority.
  source_map <- unique(
    canonical_daily[, c(
      "COMBO_METHOD_STAGE",
      "SOURCE_FILE",
      "SOURCE_PRIORITY"
    )]
  )

  key <- split(
    seq_len(nrow(source_map)),
    source_map$COMBO_METHOD_STAGE
  )

  selected <- vapply(key, function(idx) {
    x <- source_map[idx, , drop = FALSE]
    x$SOURCE_PRIORITY[is.na(x$SOURCE_PRIORITY)] <- 0L
    idx[which.max(x$SOURCE_PRIORITY)]
  }, integer(1))

  keep_source <- source_map[selected, , drop = FALSE]

  join_key <- paste(
    canonical_daily$COMBO_METHOD_STAGE,
    canonical_daily$SOURCE_FILE,
    sep = "##"
  )

  keep_key <- paste(
    keep_source$COMBO_METHOD_STAGE,
    keep_source$SOURCE_FILE,
    sep = "##"
  )

  out <- canonical_daily[join_key %in% keep_key, , drop = FALSE]
  out$COMBO_METHOD_STAGE <- NULL

  rownames(out) <- NULL
  out
}

# -----------------------------------------------------------------------------
# Improvement table
# -----------------------------------------------------------------------------

cm_improvement_vs_raw <- function(final_metrics) {
  if (!nrow(final_metrics)) return(data.frame())

  error_metrics <- intersect(
    c(
      "Daily_MAE",
      "Daily_RMSE",
      "P90_MAE",
      "P95_MAE",
      "P99_MAE",
      "WetProp_MAE",
      "WetProp_RMSE",
      "MonthlyTotal_MAE",
      "MonthlyTotal_RMSE"
    ),
    names(final_metrics)
  )

  abs_error_metrics <- intersect(
    c(
      "Daily_Bias",
      "P90_Bias",
      "P95_Bias",
      "P99_Bias",
      "WetProp_Bias",
      "MonthlyTotal_Bias"
    ),
    names(final_metrics)
  )

  key <- cm_key(
    final_metrics$LOCATION_ID,
    final_metrics$GCM,
    final_metrics$SSP
  )

  raw <- final_metrics[final_metrics$METHOD == "RAW_GCM", , drop = FALSE]
  if (!nrow(raw)) return(data.frame())

  raw_key <- cm_key(raw$LOCATION_ID, raw$GCM, raw$SSP)

  corrected <- final_metrics[
    final_metrics$METHOD != "RAW_GCM",
    ,
    drop = FALSE
  ]

  rows <- vector("list", nrow(corrected))

  for (i in seq_len(nrow(corrected))) {
    r <- corrected[i, , drop = FALSE]
    k <- cm_key(r$LOCATION_ID, r$GCM, r$SSP)
    j <- match(k, raw_key)

    if (is.na(j)) next

    out <- r[, c(
      "REGIME", "LOCATION_ID", "GCM", "SSP",
      "METHOD", "STAGE", "N_BOOT"
    ), drop = FALSE]

    for (m in error_metrics) {
      base <- raw[[m]][j]
      val <- r[[m]][1]

      out[[paste0(m, "_RAW")]] <- base
      out[[paste0(m, "_METHOD")]] <- val

      out[[paste0(m, "_IMPROVEMENT_PCT")]] <- if (
        is.finite(base) && base != 0 && is.finite(val)
      ) {
        100 * (base - val) / base
      } else {
        NA_real_
      }
    }

    for (m in abs_error_metrics) {
      base <- abs(raw[[m]][j])
      val <- abs(r[[m]][1])

      out[[paste0("ABS_", m, "_RAW")]] <- base
      out[[paste0("ABS_", m, "_METHOD")]] <- val

      out[[paste0("ABS_", m, "_IMPROVEMENT_PCT")]] <- if (
        is.finite(base) && base != 0 && is.finite(val)
      ) {
        100 * (base - val) / base
      } else {
        NA_real_
      }
    }

    rows[[i]] <- out
  }

  cm_rbind(rows)
}


# -----------------------------------------------------------------------------
# Major-revision convenience API
# -----------------------------------------------------------------------------
common_metric_evaluate_major_revision <- function(canonical_daily,
                                                  target_years = CM_TARGET_YEARS,
                                                  wet_threshold = CM_WET_THRESHOLD,
                                                  min_month_completeness = 0.90,
                                                  nominal_interval = 0.95,
                                                  assert_benchmarks = TRUE) {
  d <- cm_prepare_major_revision_schema(canonical_daily)
  if (isTRUE(assert_benchmarks)) cm_assert_benchmark_invariance(d)

  common_metric_evaluate_all(
    canonical_daily = d,
    target_years = target_years,
    wet_threshold = wet_threshold,
    min_month_completeness = min_month_completeness,
    extreme_probs = CM_EXTREME_PROBS,
    nominal_interval = nominal_interval
  )
}

cat("Loaded COMMON_METRIC_CORE_v1.1.0.R — six-arm major-revision metric engine.\n")
