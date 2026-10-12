# =============================================================================
# CHRIPS_GCM_CORE_v1.1.0.R
# Core functions for RAW GCM vs CHIRPS metrics + visualization
# =============================================================================
#
# Repository utility for RAW GCM vs CHIRPS diagnostics.
# Active project input is repository-relative under data/adaptive/.
# Generated outputs should be written under results/.
#
# Scientific location labels:
#   Monsoonal_1, Monsoonal_2
#   Equatorial_1, Equatorial_2
#   Local_1, Local_2
#
# Metadata resolution is intentionally robust:
#   1) parse filename/path;
#   2) parse metadata columns inside CSV;
#   3) infer regime from explicit location;
#   4) if location is absent, use the fixed project routing:
#        CESM2-WACCM / CMCC-ESM2 / TaiESM1 -> location 1
#        MRI-ESM2-0                         -> location 2
#
# Therefore files named, for example,
#   Jawa_Timur_CESM2_WACCM_ssp126_....csv
# can still be resolved as:
#   REGIME = Monsoonal
#   LOCATION_ID = Monsoonal_1
#   GCM = CESM2_WACCM
#   SSP = ssp126
#
# =============================================================================

# -----------------------------------------------------------------------------
# Utility helpers
# -----------------------------------------------------------------------------

cg_norm <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

cg_safe_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

cg_mean <- function(x) {
  x <- cg_safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

cg_sd <- function(x) {
  x <- cg_safe_num(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

cg_mae <- function(obs, sim) {
  ok <- is.finite(obs) & is.finite(sim)
  if (!any(ok)) return(NA_real_)
  mean(abs(sim[ok] - obs[ok]))
}

cg_rmse <- function(obs, sim) {
  ok <- is.finite(obs) & is.finite(sim)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((sim[ok] - obs[ok])^2))
}

cg_bias <- function(obs, sim) {
  ok <- is.finite(obs) & is.finite(sim)
  if (!any(ok)) return(NA_real_)
  mean(sim[ok] - obs[ok])
}

cg_cor <- function(obs, sim) {
  ok <- is.finite(obs) & is.finite(sim)
  if (sum(ok) < 3L) return(NA_real_)
  if (stats::sd(obs[ok]) <= 0 || stats::sd(sim[ok]) <= 0) return(NA_real_)
  suppressWarnings(stats::cor(obs[ok], sim[ok], method = "pearson"))
}

cg_kge <- function(obs, sim) {
  ok <- is.finite(obs) & is.finite(sim)
  if (sum(ok) < 3L) return(NA_real_)

  o <- obs[ok]
  s <- sim[ok]

  r <- cg_cor(o, s)
  if (!is.finite(r)) return(NA_real_)

  mo <- mean(o)
  ms <- mean(s)
  sdo <- stats::sd(o)
  sds <- stats::sd(s)

  if (!is.finite(mo) || !is.finite(ms) ||
      !is.finite(sdo) || !is.finite(sds) ||
      mo == 0 || sdo == 0) {
    return(NA_real_)
  }

  alpha <- sds / sdo
  beta  <- ms / mo

  1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
}

cg_pinball <- function(y, qhat, tau) {
  ok <- is.finite(y) & is.finite(qhat)
  if (!any(ok)) return(NA_real_)
  e <- y[ok] - qhat[ok]
  mean(pmax(tau * e, (tau - 1) * e))
}

cg_quantile <- function(x, probs) {
  x <- cg_safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(
    x,
    probs = probs,
    type = 7,
    na.rm = TRUE,
    names = FALSE
  ))
}

cg_qname <- function(tau) {
  paste0("Q", sprintf("%02d", round(100 * tau)))
}

cg_parse_date <- function(x) {
  if (inherits(x, "Date")) return(x)

  z <- as.character(x)
  out <- suppressWarnings(as.Date(z))

  fmts <- c(
    "%Y/%m/%d",
    "%d/%m/%Y",
    "%m/%d/%Y",
    "%d-%m-%Y",
    "%m-%d-%Y",
    "%Y%m%d"
  )

  for (fmt in fmts) {
    miss <- is.na(out) & !is.na(z)
    if (!any(miss)) break
    out[miss] <- suppressWarnings(as.Date(z[miss], format = fmt))
  }

  out
}

cg_wet <- function(x, threshold = 0.1, rule = "GT") {
  x <- cg_safe_num(x)

  if (toupper(rule) == "GT") {
    return(is.finite(x) & x > threshold)
  }

  if (toupper(rule) == "GE") {
    return(is.finite(x) & x >= threshold)
  }

  stop("WET_DAY_RULE must be 'GT' or 'GE'.", call. = FALSE)
}

cg_days_in_month_noleap <- function(month) {
  c(31L, 28L, 31L, 30L, 31L, 30L,
    31L, 31L, 30L, 31L, 30L, 31L)[as.integer(month)]
}

cg_rbind <- function(lst) {
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

# -----------------------------------------------------------------------------
# Metadata parsing
# -----------------------------------------------------------------------------

cg_parse_regime <- function(x) {
  z <- cg_norm(x)

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

cg_parse_gcm <- function(x) {
  z <- cg_norm(x)

  if (grepl("CESM2_WACCM|CESM2WACCM", z)) return("CESM2_WACCM")
  if (grepl("CMCC_ESM2|CMCCESM2", z))     return("CMCC_ESM2")
  if (grepl("MRI_ESM2_0|MRIESM2_0|MRIESM20", z)) return("MRI_ESM2_0")
  if (grepl("TAIESM1|TAI_ESM1", z))       return("TaiESM1")

  NA_character_
}

cg_parse_ssp <- function(x) {
  z <- cg_norm(x)

  if (grepl("SSP126|SSP1_2_6", z)) return("ssp126")
  if (grepl("SSP245|SSP2_4_5", z)) return("ssp245")
  if (grepl("SSP370|SSP3_7_0", z)) return("ssp370")
  if (grepl("SSP585|SSP5_8_5", z)) return("ssp585")

  NA_character_
}

cg_parse_location <- function(x, regime = NA_character_) {
  z <- cg_norm(x)

  explicit <- list(
    Monsoonal_1 = c("MONSOONAL_1", "MONSUN_1", "JAWA_TIMUR_1", "EAST_JAVA_1"),
    Monsoonal_2 = c("MONSOONAL_2", "MONSUN_2", "JAWA_TIMUR_2", "EAST_JAVA_2"),
    Equatorial_1 = c("EQUATORIAL_1", "EKUATORIAL_1", "SUMATERA_BARAT_1", "WEST_SUMATRA_1"),
    Equatorial_2 = c("EQUATORIAL_2", "EKUATORIAL_2", "SUMATERA_BARAT_2", "WEST_SUMATRA_2"),
    Local_1 = c("LOCAL_1", "LOKAL_1", "MALUKU_1"),
    Local_2 = c("LOCAL_2", "LOKAL_2", "MALUKU_2")
  )

  for (nm in names(explicit)) {
    if (any(vapply(explicit[[nm]], function(p) grepl(p, z, fixed = TRUE), logical(1)))) {
      return(nm)
    }
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

cg_regime_from_location <- function(location_id) {
  if (is.na(location_id)) return(NA_character_)
  sub("_[12]$", "", location_id)
}

cg_fixed_routing_location <- function(regime, gcm) {
  if (is.na(regime) || is.na(gcm)) return(NA_character_)

  loc_no <- switch(
    gcm,
    "CESM2_WACCM" = 1L,
    "CMCC_ESM2"   = 1L,
    "MRI_ESM2_0"  = 2L,
    "TaiESM1"     = 1L,
    NA_integer_
  )

  if (is.na(loc_no)) return(NA_character_)
  paste0(regime, "_", loc_no)
}

cg_find_col <- function(nms, candidates) {
  nn <- cg_norm(nms)
  cc <- cg_norm(candidates)

  j <- which(nn %in% cc)
  if (!length(j)) return(NA_integer_)
  j[1]
}

cg_first_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (!length(x)) NA_character_ else x[1]
}

cg_read_header_sample <- function(path, n = 200L) {
  tryCatch(
    utils::read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      nrows = n
    ),
    error = function(e) NULL
  )
}

cg_metadata_from_columns <- function(path) {
  h <- cg_read_header_sample(path)

  out <- list(
    regime = NA_character_,
    location = NA_character_,
    gcm = NA_character_,
    ssp = NA_character_,
    detail = character()
  )

  if (is.null(h) || !nrow(h)) return(out)

  nms <- names(h)

  regime_j <- cg_find_col(
    nms,
    c("region", "regime", "rainfall_regime", "wilayah", "province", "provinsi", "area")
  )

  location_j <- cg_find_col(
    nms,
    c("location_id", "location", "lokasi", "loc", "site", "grid_id", "grid", "point")
  )

  gcm_j <- cg_find_col(
    nms,
    c("gcm", "model", "gcm_model", "model_name", "climate_model")
  )

  ssp_j <- cg_find_col(
    nms,
    c("ssp", "scenario", "experiment", "experiment_id", "scenario_id")
  )

  if (!is.na(regime_j)) {
    v <- cg_first_nonempty(h[[regime_j]])
    out$regime <- cg_parse_regime(v)
    if (!is.na(out$regime)) {
      out$detail <- c(out$detail, paste0("regime_column=", nms[regime_j]))
    }
  }

  if (!is.na(gcm_j)) {
    v <- cg_first_nonempty(h[[gcm_j]])
    out$gcm <- cg_parse_gcm(v)
    if (!is.na(out$gcm)) {
      out$detail <- c(out$detail, paste0("gcm_column=", nms[gcm_j]))
    }
  }

  if (!is.na(ssp_j)) {
    v <- cg_first_nonempty(h[[ssp_j]])
    out$ssp <- cg_parse_ssp(v)
    if (!is.na(out$ssp)) {
      out$detail <- c(out$detail, paste0("ssp_column=", nms[ssp_j]))
    }
  }

  if (!is.na(location_j)) {
    vals <- unique(h[[location_j]])
    vals <- vals[!is.na(vals)]

    for (v in vals) {
      loc <- cg_parse_location(v, out$regime)

      if (is.na(loc)) {
        vv <- trimws(as.character(v))
        if (!is.na(out$regime) && vv %in% c("1", "01")) {
          loc <- paste0(out$regime, "_1")
        }
        if (!is.na(out$regime) && vv %in% c("2", "02")) {
          loc <- paste0(out$regime, "_2")
        }
      }

      if (!is.na(loc)) {
        out$location <- loc
        out$detail <- c(out$detail, paste0("location_column=", nms[location_j]))
        break
      }
    }
  }

  if (is.na(out$regime) && !is.na(out$location)) {
    out$regime <- cg_regime_from_location(out$location)
  }

  out
}

cg_build_manifest_row <- function(path, allow_routing = TRUE) {
  txt <- normalizePath(path, winslash = "/", mustWork = FALSE)

  regime <- cg_parse_regime(txt)
  gcm <- cg_parse_gcm(txt)
  ssp <- cg_parse_ssp(txt)
  location <- cg_parse_location(txt, regime)

  src_regime <- if (!is.na(regime)) "path" else NA_character_
  src_gcm <- if (!is.na(gcm)) "path" else NA_character_
  src_ssp <- if (!is.na(ssp)) "path" else NA_character_
  src_location <- if (!is.na(location)) "path" else NA_character_

  meta_col <- cg_metadata_from_columns(path)

  if (is.na(regime) && !is.na(meta_col$regime)) {
    regime <- meta_col$regime
    src_regime <- "column"
  }

  if (is.na(gcm) && !is.na(meta_col$gcm)) {
    gcm <- meta_col$gcm
    src_gcm <- "column"
  }

  if (is.na(ssp) && !is.na(meta_col$ssp)) {
    ssp <- meta_col$ssp
    src_ssp <- "column"
  }

  if (is.na(location) && !is.na(meta_col$location)) {
    location <- meta_col$location
    src_location <- "column"
  }

  if (is.na(regime) && !is.na(location)) {
    regime <- cg_regime_from_location(location)
    src_regime <- "location"
  }

  if (is.na(location) && isTRUE(allow_routing) &&
      !is.na(regime) && !is.na(gcm)) {
    location <- cg_fixed_routing_location(regime, gcm)
    if (!is.na(location)) src_location <- "fixed_project_routing"
  }

  missing_meta <- c(
    if (is.na(regime)) "REGIME",
    if (is.na(location)) "LOCATION_ID",
    if (is.na(gcm)) "GCM",
    if (is.na(ssp)) "SSP"
  )

  status <- if (!length(missing_meta)) "OK" else {
    paste0("UNRESOLVED_", paste(missing_meta, collapse = "_"))
  }

  data.frame(
    REGIME = regime,
    LOCATION_ID = location,
    GCM = gcm,
    SSP = ssp,
    REGIME_SOURCE = src_regime,
    LOCATION_SOURCE = src_location,
    GCM_SOURCE = src_gcm,
    SSP_SOURCE = src_ssp,
    FILE = txt,
    FILE_NAME = basename(path),
    PARSE_STATUS = status,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Paired rainfall column detection and loading
# -----------------------------------------------------------------------------

cg_detect_data_columns <- function(df) {
  nms <- names(df)

  date_j <- cg_find_col(
    nms,
    c("date", "tanggal", "time", "datetime")
  )

  obs_j <- cg_find_col(
    nms,
    c(
      "chirps", "chrips", "obs", "observed", "observation",
      "precip_obs", "rain_obs", "rainfall_obs",
      "pr_obs", "chirps_rainfall", "chirps_precip"
    )
  )

  sim_j <- cg_find_col(
    nms,
    c(
      "gcm", "raw_gcm", "sim", "simulation", "model_value",
      "precip_gcm", "rain_gcm", "rainfall_gcm",
      "pr_gcm", "pred", "prediction"
    )
  )

  # Fallback: if CHIRPS is identifiable and there is exactly one other
  # plausible numeric rainfall column, use it as GCM.
  if (!is.na(obs_j) && is.na(sim_j)) {
    numeric_like <- which(vapply(df, function(x) {
      y <- suppressWarnings(as.numeric(as.character(x)))
      sum(is.finite(y)) >= max(5L, floor(0.5 * length(y)))
    }, logical(1)))

    excluded <- unique(c(
      date_j,
      obs_j,
      cg_find_col(nms, c("year")),
      cg_find_col(nms, c("month")),
      cg_find_col(nms, c("day")),
      cg_find_col(nms, c("lat", "latitude")),
      cg_find_col(nms, c("lon", "long", "longitude"))
    ))

    excluded <- excluded[is.finite(excluded)]
    candidates <- setdiff(numeric_like, excluded)

    if (length(candidates) == 1L) {
      sim_j <- candidates[1]
    }
  }

  list(
    date = date_j,
    obs = obs_j,
    sim = sim_j
  )
}

cg_load_pair <- function(path,
                         target_years = c(2015L, 2025L),
                         no_leap = TRUE) {
  raw <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  cols <- cg_detect_data_columns(raw)

  if (any(is.na(unlist(cols)))) {
    stop(
      paste0(
        "Could not detect date / CHIRPS / GCM rainfall columns. ",
        "Available columns: ",
        paste(names(raw), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  d <- data.frame(
    date = cg_parse_date(raw[[cols$date]]),
    chirps = cg_safe_num(raw[[cols$obs]]),
    gcm = cg_safe_num(raw[[cols$sim]]),
    stringsAsFactors = FALSE
  )

  d <- d[!is.na(d$date), , drop = FALSE]
  d <- d[order(d$date), , drop = FALSE]
  d <- d[!duplicated(d$date), , drop = FALSE]

  d$year <- as.integer(format(d$date, "%Y"))
  d$month <- as.integer(format(d$date, "%m"))
  d$day <- as.integer(format(d$date, "%d"))

  if (isTRUE(no_leap)) {
    d <- d[!(d$month == 2L & d$day == 29L), , drop = FALSE]
  }

  d <- d[
    d$year >= target_years[1] &
      d$year <= target_years[2],
    ,
    drop = FALSE
  ]

  d$pair_finite <- is.finite(d$chirps) & is.finite(d$gcm)

  if (!nrow(d)) {
    stop(
      sprintf(
        "No usable data inside target period %d-%d.",
        target_years[1],
        target_years[2]
      ),
      call. = FALSE
    )
  }

  d
}

# -----------------------------------------------------------------------------
# Monthly, annual, and QQ structures
# -----------------------------------------------------------------------------

cg_make_monthly <- function(daily,
                            wet_threshold = 0.1,
                            wet_rule = "GT",
                            taus = c(0.10, 0.50, 0.90, 0.95, 0.99),
                            min_completeness = 0.90) {
  key <- interaction(
    daily$year,
    daily$month,
    drop = TRUE,
    lex.order = TRUE
  )

  spl <- split(seq_len(nrow(daily)), key)
  rows <- vector("list", length(spl))
  ii <- 0L

  for (idx in spl) {
    ii <- ii + 1L
    g <- daily[idx, , drop = FALSE]

    paired <- is.finite(g$chirps) & is.finite(g$gcm)

    obs_all <- g$chirps[paired]
    sim_all <- g$gcm[paired]

    obs_wet <- obs_all[cg_wet(obs_all, wet_threshold, wet_rule)]
    sim_wet <- sim_all[cg_wet(sim_all, wet_threshold, wet_rule)]

    expected_days <- cg_days_in_month_noleap(g$month[1])
    paired_days <- sum(paired)
    completeness <- paired_days / expected_days

    q_obs <- cg_quantile(obs_wet, taus)
    q_sim <- cg_quantile(sim_wet, taus)

    upper_obs <- cg_quantile(obs_all, c(0.95, 0.99))
    upper_sim <- cg_quantile(sim_all, c(0.95, 0.99))

    row <- data.frame(
      year = g$year[1],
      month = g$month[1],
      expected_days = expected_days,
      paired_days = paired_days,
      completeness = completeness,
      chirps_total = if (paired_days) sum(obs_all) else NA_real_,
      gcm_total = if (paired_days) sum(sim_all) else NA_real_,
      wet_prop_chirps = if (paired_days) mean(cg_wet(obs_all, wet_threshold, wet_rule)) else NA_real_,
      wet_prop_gcm = if (paired_days) mean(cg_wet(sim_all, wet_threshold, wet_rule)) else NA_real_,
      P95_ALL_CHIRPS = upper_obs[1],
      P95_ALL_GCM = upper_sim[1],
      P99_ALL_CHIRPS = upper_obs[2],
      P99_ALL_GCM = upper_sim[2],
      stringsAsFactors = FALSE
    )

    for (j in seq_along(taus)) {
      qn <- cg_qname(taus[j])
      row[[paste0(qn, "_CHIRPS_WET")]] <- q_obs[j]
      row[[paste0(qn, "_GCM_WET")]] <- q_sim[j]
    }

    rows[[ii]] <- row
  }

  out <- cg_rbind(rows)
  out <- out[order(out$year, out$month), , drop = FALSE]
  out$valid_total_month <- out$completeness >= min_completeness
  out$monthly_total_bias <- out$gcm_total - out$chirps_total
  rownames(out) <- NULL
  out
}

cg_make_annual <- function(monthly) {
  x <- monthly[
    monthly$valid_total_month &
      is.finite(monthly$chirps_total) &
      is.finite(monthly$gcm_total),
    ,
    drop = FALSE
  ]

  if (!nrow(x)) return(data.frame())

  spl <- split(seq_len(nrow(x)), x$year)

  cg_rbind(lapply(spl, function(idx) {
    g <- x[idx, , drop = FALSE]
    data.frame(
      year = g$year[1],
      n_valid_months = nrow(g),
      chirps_annual_total = sum(g$chirps_total),
      gcm_annual_total = sum(g$gcm_total),
      annual_bias = sum(g$gcm_total) - sum(g$chirps_total),
      stringsAsFactors = FALSE
    )
  }))
}

cg_make_qq <- function(daily,
                       wet_threshold = 0.1,
                       wet_rule = "GT",
                       probs = unique(c(seq(0.01, 0.90, 0.01), seq(0.91, 0.99, 0.01)))) {
  p <- daily[daily$pair_finite, c("chirps", "gcm"), drop = FALSE]

  obs <- p$chirps[cg_wet(p$chirps, wet_threshold, wet_rule)]
  sim <- p$gcm[cg_wet(p$gcm, wet_threshold, wet_rule)]

  data.frame(
    prob = probs,
    CHIRPS = cg_quantile(obs, probs),
    GCM = cg_quantile(sim, probs),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Pooled-calendar-month extreme metrics
# -----------------------------------------------------------------------------
#
# FINAL cross-method P90/P95/P99 definition:
#   For each calendar month m = 1,...,12, pool ALL daily values from all
#   evaluation years (2015-2025 by default), compute the percentile once for
#   CHIRPS and once for the prediction, then calculate MAE/RMSE/Bias across the
#   12 calendar-month percentile pairs.
#
# This is intentionally different from a year-month diagnostic (132 values).
# The pooled definition is the locked definition for comparison with
# QMB-Kokic, QDM, and ACT-QMB.

cg_pooled_calendar_extremes <- function(daily,
                                        probs = c(0.90, 0.95, 0.99)) {
  rows <- vector("list", 12L)

  for (m in 1:12) {
    d <- daily[
      daily$month == m &
        is.finite(daily$chirps) &
        is.finite(daily$gcm),
      ,
      drop = FALSE
    ]

    qo <- cg_quantile(d$chirps, probs)
    qs <- cg_quantile(d$gcm, probs)

    row <- data.frame(
      month = m,
      month_name = month.abb[m],
      n_days = nrow(d),
      stringsAsFactors = FALSE
    )

    for (j in seq_along(probs)) {
      qn <- paste0("P", round(100 * probs[j]))
      row[[paste0(qn, "_CHIRPS")]] <- qo[j]
      row[[paste0(qn, "_GCM")]] <- qs[j]
      row[[paste0(qn, "_ERROR")]] <- qs[j] - qo[j]
      row[[paste0(qn, "_ABS_ERROR")]] <- abs(qs[j] - qo[j])
    }

    rows[[m]] <- row
  }

  cg_rbind(rows)
}

cg_extreme_summary <- function(pooled_table,
                               probs = c(0.90, 0.95, 0.99)) {
  out <- data.frame(
    ExtremeMetricDefinition = "pooled_calendar_month_across_target_years",
    ExtremeMetricNMonths = nrow(pooled_table),
    stringsAsFactors = FALSE
  )

  for (p in probs) {
    qn <- paste0("P", round(100 * p))
    obs <- pooled_table[[paste0(qn, "_CHIRPS")]]
    sim <- pooled_table[[paste0(qn, "_GCM")]]

    out[[paste0(qn, "_MAE")]] <- cg_mae(obs, sim)
    out[[paste0(qn, "_RMSE")]] <- cg_rmse(obs, sim)
    out[[paste0(qn, "_Bias")]] <- cg_bias(obs, sim)
  }

  out
}

# -----------------------------------------------------------------------------
# Metrics
# -----------------------------------------------------------------------------

cg_quantile_metrics <- function(daily,
                                monthly,
                                taus = c(0.10, 0.50, 0.90, 0.95, 0.99),
                                wet_threshold = 0.1,
                                wet_rule = "GT") {
  rows <- vector("list", length(taus))

  for (j in seq_along(taus)) {
    tau <- taus[j]
    qn <- cg_qname(tau)

    obs_col <- paste0(qn, "_CHIRPS_WET")
    sim_col <- paste0(qn, "_GCM_WET")

    obs_q <- monthly[[obs_col]]
    sim_q <- monthly[[sim_col]]

    look <- monthly[, c("year", "month", sim_col), drop = FALSE]

    d <- merge(
      daily[, c("year", "month", "chirps"), drop = FALSE],
      look,
      by = c("year", "month"),
      all.x = TRUE,
      sort = FALSE
    )

    d <- d[
      cg_wet(d$chirps, wet_threshold, wet_rule) &
        is.finite(d[[sim_col]]),
      ,
      drop = FALSE
    ]

    rows[[j]] <- data.frame(
      tau = tau,
      quantile = qn,
      n_months = sum(is.finite(obs_q) & is.finite(sim_q)),
      n_wet_days_pinball = nrow(d),
      pinball_loss = cg_pinball(d$chirps, d[[sim_col]], tau),
      MAE = cg_mae(obs_q, sim_q),
      RMSE = cg_rmse(obs_q, sim_q),
      Bias = cg_bias(obs_q, sim_q),
      stringsAsFactors = FALSE
    )
  }

  cg_rbind(rows)
}

cg_occurrence_metrics <- function(monthly) {
  wet_obs <- monthly$wet_prop_chirps
  wet_sim <- monthly$wet_prop_gcm

  dry_obs <- 1 - wet_obs
  dry_sim <- 1 - wet_sim

  data.frame(
    n_months = sum(is.finite(wet_obs) & is.finite(wet_sim)),
    Wet_MAE = cg_mae(wet_obs, wet_sim),
    Wet_RMSE = cg_rmse(wet_obs, wet_sim),
    Wet_Bias = cg_bias(wet_obs, wet_sim),
    Dry_MAE = cg_mae(dry_obs, dry_sim),
    Dry_RMSE = cg_rmse(dry_obs, dry_sim),
    Dry_Bias = cg_bias(dry_obs, dry_sim),
    stringsAsFactors = FALSE
  )
}

cg_main_metrics <- function(daily, monthly) {
  d <- daily[daily$pair_finite, , drop = FALSE]

  vm <- monthly$valid_total_month &
    is.finite(monthly$chirps_total) &
    is.finite(monthly$gcm_total)

  pooled <- cg_pooled_calendar_extremes(
    daily,
    probs = c(0.90, 0.95, 0.99)
  )

  extreme <- cg_extreme_summary(
    pooled,
    probs = c(0.90, 0.95, 0.99)
  )

  daily_part <- data.frame(
    n_days = nrow(d),

    Daily_MAE = cg_mae(d$chirps, d$gcm),
    Daily_RMSE = cg_rmse(d$chirps, d$gcm),
    Daily_Bias = cg_bias(d$chirps, d$gcm),
    Daily_r = cg_cor(d$chirps, d$gcm),
    Daily_KGE = cg_kge(d$chirps, d$gcm),

    MonthlyTotal_n = sum(vm),
    MonthlyTotal_MAE = cg_mae(monthly$chirps_total[vm], monthly$gcm_total[vm]),
    MonthlyTotal_RMSE = cg_rmse(monthly$chirps_total[vm], monthly$gcm_total[vm]),
    MonthlyTotal_Bias = cg_bias(monthly$chirps_total[vm], monthly$gcm_total[vm]),
    MonthlyTotal_r = cg_cor(monthly$chirps_total[vm], monthly$gcm_total[vm]),
    MonthlyTotal_KGE = cg_kge(monthly$chirps_total[vm], monthly$gcm_total[vm]),

    stringsAsFactors = FALSE
  )

  cbind(daily_part, extreme)
}


cg_calendar_month_metrics <- function(daily, monthly,
                                      wet_threshold = 0.1,
                                      wet_rule = "GT") {
  rows <- vector("list", 12L)

  for (m in 1:12) {
    d <- daily[daily$month == m & daily$pair_finite, , drop = FALSE]
    mm <- monthly[monthly$month == m, , drop = FALSE]

    valid <- mm$valid_total_month &
      is.finite(mm$chirps_total) &
      is.finite(mm$gcm_total)

    rows[[m]] <- data.frame(
      month = m,
      month_name = month.abb[m],
      n_days = nrow(d),

      Daily_MAE = cg_mae(d$chirps, d$gcm),
      Daily_RMSE = cg_rmse(d$chirps, d$gcm),
      Daily_Bias = cg_bias(d$chirps, d$gcm),
      Daily_r = cg_cor(d$chirps, d$gcm),
      Daily_KGE = cg_kge(d$chirps, d$gcm),

      WetProp_MAE = cg_mae(mm$wet_prop_chirps, mm$wet_prop_gcm),
      WetProp_RMSE = cg_rmse(mm$wet_prop_chirps, mm$wet_prop_gcm),
      WetProp_Bias = cg_bias(mm$wet_prop_chirps, mm$wet_prop_gcm),

      MonthlyTotal_MAE = cg_mae(mm$chirps_total[valid], mm$gcm_total[valid]),
      MonthlyTotal_RMSE = cg_rmse(mm$chirps_total[valid], mm$gcm_total[valid]),
      MonthlyTotal_Bias = cg_bias(mm$chirps_total[valid], mm$gcm_total[valid]),
      MonthlyTotal_r = cg_cor(mm$chirps_total[valid], mm$gcm_total[valid]),

      P95_MAE = cg_mae(mm$P95_ALL_CHIRPS, mm$P95_ALL_GCM),
      P95_RMSE = cg_rmse(mm$P95_ALL_CHIRPS, mm$P95_ALL_GCM),
      P95_Bias = cg_bias(mm$P95_ALL_CHIRPS, mm$P95_ALL_GCM),

      P99_MAE = cg_mae(mm$P99_ALL_CHIRPS, mm$P99_ALL_GCM),
      P99_RMSE = cg_rmse(mm$P99_ALL_CHIRPS, mm$P99_ALL_GCM),
      P99_Bias = cg_bias(mm$P99_ALL_CHIRPS, mm$P99_ALL_GCM),

      stringsAsFactors = FALSE
    )
  }

  cg_rbind(rows)
}

# -----------------------------------------------------------------------------
# Grouped summary
# -----------------------------------------------------------------------------

cg_group_summary <- function(df, groups, metrics) {
  if (!nrow(df)) return(data.frame())

  key <- interaction(
    df[, groups, drop = FALSE],
    drop = TRUE,
    lex.order = TRUE
  )

  spl <- split(seq_len(nrow(df)), key)

  rows <- lapply(spl, function(idx) {
    d <- df[idx, , drop = FALSE]

    row <- as.list(d[1, groups, drop = FALSE])
    row$n_combinations <- nrow(d)

    for (m in metrics) {
      row[[paste0(m, "_mean")]] <- cg_mean(d[[m]])
      row[[paste0(m, "_sd")]] <- cg_sd(d[[m]])
    }

    as.data.frame(
      row,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  cg_rbind(rows)
}

# -----------------------------------------------------------------------------
# Plot helpers
# -----------------------------------------------------------------------------

cg_theme_pub <- function(base_size = 10.5, base_family = "Helvetica") {
  ggplot2::theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_line(
        linewidth = 0.25,
        colour = "grey88"
      ),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    )
}

cg_save_plot <- function(plot, stem,
                         width = 8,
                         height = 5,
                         save_pdf = TRUE,
                         save_tiff = TRUE,
                         tiff_dpi = 600) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(save_pdf)) {
    ggplot2::ggsave(
      paste0(stem, ".pdf"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      bg = "white",
      limitsize = FALSE
    )
  }

  if (isTRUE(save_tiff)) {
    ggplot2::ggsave(
      paste0(stem, ".tiff"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = tiff_dpi,
      compression = "lzw",
      bg = "white",
      limitsize = FALSE
    )
  }

  invisible(stem)
}

cg_build_combo_plots <- function(daily, monthly, annual, qq,
                                 location_id, gcm, ssp) {
  tag <- paste(location_id, gcm, toupper(ssp), sep = " | ")

  clim <- monthly[monthly$valid_total_month, , drop = FALSE]

  clim_long <- rbind(
    data.frame(
      year = clim$year,
      month = clim$month,
      Source = "CHIRPS",
      Total = clim$chirps_total
    ),
    data.frame(
      year = clim$year,
      month = clim$month,
      Source = "RAW GCM",
      Total = clim$gcm_total
    )
  )

  clim_sum <- clim_long |>
    dplyr::group_by(Source, month) |>
    dplyr::summarise(
      Mean = mean(Total, na.rm = TRUE),
      Q25 = stats::quantile(Total, 0.25, na.rm = TRUE),
      Q75 = stats::quantile(Total, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  p1 <- ggplot2::ggplot(
    clim_sum,
    ggplot2::aes(x = month, y = Mean, linetype = Source)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = Q25, ymax = Q75, group = Source),
      alpha = 0.12
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_x_continuous(breaks = 1:12, labels = month.abb) +
    ggplot2::labs(
      title = "Monthly rainfall climatology",
      subtitle = tag,
      x = NULL,
      y = "Monthly total (mm)"
    ) +
    cg_theme_pub()

  if (nrow(annual)) {
    ann_long <- rbind(
      data.frame(
        year = annual$year,
        Source = "CHIRPS",
        Total = annual$chirps_annual_total
      ),
      data.frame(
        year = annual$year,
        Source = "RAW GCM",
        Total = annual$gcm_annual_total
      )
    )
  } else {
    ann_long <- data.frame(year = numeric(), Source = character(), Total = numeric())
  }

  p2 <- ggplot2::ggplot(
    ann_long,
    ggplot2::aes(x = year, y = Total, linetype = Source)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::labs(
      title = "Annual rainfall totals",
      subtitle = tag,
      x = "Year",
      y = "Annual total (mm)"
    ) +
    cg_theme_pub()

  scat <- monthly[
    monthly$valid_total_month &
      is.finite(monthly$chirps_total) &
      is.finite(monthly$gcm_total),
    ,
    drop = FALSE
  ]

  ann_txt <- sprintf(
    "r = %.3f\nRMSE = %.2f mm\nBias = %.2f mm",
    cg_cor(scat$chirps_total, scat$gcm_total),
    cg_rmse(scat$chirps_total, scat$gcm_total),
    cg_bias(scat$chirps_total, scat$gcm_total)
  )

  p3 <- ggplot2::ggplot(
    scat,
    ggplot2::aes(x = chirps_total, y = gcm_total)
  ) +
    ggplot2::geom_point(alpha = 0.55, size = 1.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = ann_txt,
      hjust = 1.05,
      vjust = -0.25,
      size = 3
    ) +
    ggplot2::labs(
      title = "Monthly-total agreement",
      subtitle = tag,
      x = "CHIRPS monthly total (mm)",
      y = "RAW GCM monthly total (mm)"
    ) +
    cg_theme_pub()

  p4 <- ggplot2::ggplot(
    qq,
    ggplot2::aes(x = CHIRPS, y = GCM)
  ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = "Wet-day Q-Q",
      subtitle = tag,
      x = "CHIRPS wet-day quantile (mm/day)",
      y = "RAW GCM wet-day quantile (mm/day)"
    ) +
    cg_theme_pub()

  hm <- monthly
  hm$month_f <- factor(hm$month, levels = 1:12, labels = month.abb)

  p5 <- ggplot2::ggplot(
    hm,
    ggplot2::aes(
      x = month_f,
      y = factor(year),
      fill = monthly_total_bias
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(midpoint = 0, name = "Bias\n(mm)") +
    ggplot2::labs(
      title = "Year-month total bias",
      subtitle = "RAW GCM minus CHIRPS",
      x = NULL,
      y = "Year"
    ) +
    cg_theme_pub() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "right")

  wet_long <- rbind(
    data.frame(
      month = monthly$month,
      Source = "CHIRPS",
      WetProp = monthly$wet_prop_chirps
    ),
    data.frame(
      month = monthly$month,
      Source = "RAW GCM",
      WetProp = monthly$wet_prop_gcm
    )
  )

  wet_sum <- wet_long |>
    dplyr::group_by(Source, month) |>
    dplyr::summarise(
      WetProp = mean(WetProp, na.rm = TRUE),
      .groups = "drop"
    )

  p6 <- ggplot2::ggplot(
    wet_sum,
    ggplot2::aes(x = month, y = WetProp, linetype = Source)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_x_continuous(breaks = 1:12, labels = month.abb) +
    ggplot2::labs(
      title = "Monthly wet-day frequency",
      subtitle = tag,
      x = NULL,
      y = "Wet-day proportion"
    ) +
    cg_theme_pub()

  list(
    climatology = p1,
    annual = p2,
    monthly_scatter = p3,
    wet_qq = p4,
    bias_heatmap = p5,
    wet_frequency = p6
  )
}
