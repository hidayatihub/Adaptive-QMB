# =============================================================================
# QMB KOKIC - SINGLE REGION + SINGLE GCM (4 SSP) - R VERSION
# REVISION: target sparse-GCM fallback + robust handling/diagnostics for degenerate knots
# =============================================================================
# Baseline Quantile Matching Bootstrap (QMB) for daily rainfall.
# R translation of qmb_kokic_single_region_gcm.py with the SAME input structure:
#     date, chirps, gcm
#
# Design:
#     Base / calibration : 1990-2014
#     Target / evaluation: 2015-2025
#     Quantiles           : 0.10, 0.50, 0.90, 0.95, 0.99
#     Stage 1             : monthly wet-day quantile projection (linear regression)
#     Stage 2             : Kokic quantile matching + linear interpolation
#                           + residual/error correction
#     Stage 3             : seasonally matched month-block bootstrap
#
# IMPORTANT:
# - Stage-1 is fitted ONCE per region-GCM using 1990-2014.
# - The same fit is applied to all four SSP target series, 2015-2025.
# - The four SSP files may each contain the full 1990-2025 period.
# - The script verifies that the 1990-2014 historical portions are identical.
# - Uses BASE R only (no external packages required).
#
# Reference:
# Kokic, P., Jin, H., & Crimp, S. (2013). Climate Dynamics, 41, 853-866.
# DOI: 10.1007/s00382-013-1791-z
# =============================================================================

# =============================================================================
# BASIC UTILITIES
# =============================================================================

.qname <- function(tau) {
  s <- sprintf("%.3f", as.numeric(tau))
  s <- sub("0+$", "", s)
  s <- sub("\\.$", "", s)
  paste0("q", s)
}

.slug <- function(text) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(text))
  gsub("^_+|_+$", "", x)
}

.elapsed <- function(start_time) {
  as.numeric(proc.time()[["elapsed"]] - start_time)
}

.logit <- function(p, eps = 1e-6) {
  p <- as.numeric(p)
  p <- pmin(pmax(p, eps), 1 - eps)
  log(p / (1 - p))
}

.inv_logit <- function(x) {
  x <- as.numeric(x)
  out <- numeric(length(x))
  pos <- x >= 0
  out[pos] <- 1 / (1 + exp(-x[pos]))
  ex <- exp(x[!pos])
  out[!pos] <- ex / (1 + ex)
  out
}

.fit_ols <- function(x, y, min_n = 8L) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < min_n) {
    stop(sprintf("Training sample too small: n=%d < %d", n, min_n), call. = FALSE)
  }
  X <- cbind(1, x)
  # SVD least-squares solution, matching numpy.linalg.lstsq behavior more closely
  # than qr.solve when the design is close to rank-deficient.
  z <- svd(X)
  tol <- max(dim(X)) * max(z$d) * .Machine$double.eps
  keep <- z$d > tol
  if (!any(keep)) stop("OLS design matrix has zero numerical rank", call. = FALSE)
  beta <- z$v[, keep, drop = FALSE] %*% ((crossprod(z$u[, keep, drop = FALSE], y)) / z$d[keep])
  beta <- as.numeric(beta)
  list(intercept = beta[1], slope = beta[2], n = n)
}

.mae <- function(y, yhat) {
  y <- as.numeric(y); yhat <- as.numeric(yhat)
  ok <- is.finite(y) & is.finite(yhat)
  if (!any(ok)) return(NA_real_)
  mean(abs(y[ok] - yhat[ok]))
}

.rmse <- function(y, yhat) {
  y <- as.numeric(y); yhat <- as.numeric(yhat)
  ok <- is.finite(y) & is.finite(yhat)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((y[ok] - yhat[ok])^2))
}

.bias <- function(y, yhat) {
  y <- as.numeric(y); yhat <- as.numeric(yhat)
  ok <- is.finite(y) & is.finite(yhat)
  if (!any(ok)) return(NA_real_)
  mean(yhat[ok] - y[ok])
}

.pinball <- function(y, qhat, tau) {
  y <- as.numeric(y); qhat <- as.numeric(qhat)
  ok <- is.finite(y) & is.finite(qhat)
  if (!any(ok)) return(NA_real_)
  e <- y[ok] - qhat[ok]
  mean(pmax(tau * e, (tau - 1) * e))
}

.quantile_safe <- function(values, probs) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values) == 0L) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(values, probs = probs, names = FALSE, type = 7, na.rm = TRUE))
}

.parse_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  out <- suppressWarnings(as.Date(as.character(x)))
  if (anyNA(out)) {
    formats <- c("%Y/%m/%d", "%d/%m/%Y", "%m/%d/%Y", "%d-%m-%Y", "%m-%d-%Y")
    raw <- as.character(x)
    for (fmt in formats) {
      miss <- is.na(out) & !is.na(raw)
      if (!any(miss)) break
      out[miss] <- suppressWarnings(as.Date(raw[miss], format = fmt))
    }
  }
  out
}

.safe_rbind <- function(lst) {
  if (length(lst) == 0L) return(data.frame())
  do.call(rbind, lst)
}

# =============================================================================
# DATA PREPARATION
# =============================================================================

prepare_daily_pair <- function(df,
                               date_col = "date",
                               obs_col = "chirps",
                               gcm_col = "gcm",
                               no_leap = TRUE) {
  required <- c(date_col, obs_col, gcm_col)
  missing <- setdiff(required, names(df))
  if (length(missing) > 0L) {
    stop(sprintf("Missing required columns: %s. Available: %s",
                 paste(missing, collapse = ", "), paste(names(df), collapse = ", ")),
         call. = FALSE)
  }

  out <- df[, required, drop = FALSE]
  out[[date_col]] <- .parse_date(out[[date_col]])
  out[[obs_col]] <- suppressWarnings(as.numeric(out[[obs_col]]))
  out[[gcm_col]] <- suppressWarnings(as.numeric(out[[gcm_col]]))
  out <- out[!is.na(out[[date_col]]), , drop = FALSE]
  out <- out[order(out[[date_col]]), , drop = FALSE]
  out <- out[!duplicated(out[[date_col]]), , drop = FALSE]

  if (isTRUE(no_leap)) {
    mm <- as.integer(format(out[[date_col]], "%m"))
    dd <- as.integer(format(out[[date_col]], "%d"))
    out <- out[!(mm == 2L & dd == 29L), , drop = FALSE]
  }

  out$year <- as.integer(format(out[[date_col]], "%Y"))
  out$month <- as.integer(format(out[[date_col]], "%m"))
  out$day <- as.integer(format(out[[date_col]], "%d"))
  rownames(out) <- NULL
  out
}

load_pair_csv <- function(path,
                          date_col = "date",
                          obs_col = "chirps",
                          gcm_col = "gcm",
                          no_leap = TRUE) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  prepare_daily_pair(df, date_col, obs_col, gcm_col, no_leap)
}

split_period <- function(daily,
                         base_years = c(1990L, 2014L),
                         target_years = c(2015L, 2025L)) {
  base <- daily[daily$year >= base_years[1] & daily$year <= base_years[2], , drop = FALSE]
  target <- daily[daily$year >= target_years[1] & daily$year <= target_years[2], , drop = FALSE]
  rownames(base) <- NULL; rownames(target) <- NULL
  if (nrow(base) == 0L) stop(sprintf("No base-period data found for %d-%d", base_years[1], base_years[2]), call. = FALSE)
  if (nrow(target) == 0L) stop(sprintf("No target-period data found for %d-%d", target_years[1], target_years[2]), call. = FALSE)
  list(base = base, target = target)
}

.expected_noleap_days <- function(years) {
  (as.integer(years[2]) - as.integer(years[1]) + 1L) * 365L
}

validate_period_coverage <- function(df, years, label,
                                     no_leap = TRUE,
                                     strict = FALSE,
                                     obs_col = "chirps") {
  sub <- df[df$year >= years[1] & df$year <= years[2], , drop = FALSE]
  expected <- if (isTRUE(no_leap)) .expected_noleap_days(years) else NA_integer_
  info <- data.frame(
    label = label,
    start = if (nrow(sub)) min(sub$date) else as.Date(NA),
    end = if (nrow(sub)) max(sub$date) else as.Date(NA),
    n_rows = nrow(sub),
    expected_rows_noleap = expected,
    n_missing_obs = if (obs_col %in% names(sub)) sum(is.na(sub[[obs_col]])) else NA_integer_,
    stringsAsFactors = FALSE
  )
  if (isTRUE(strict) && isTRUE(no_leap) && nrow(sub) != expected) {
    stop(sprintf("%s: expected %d no-leap rows, found %d", label, expected, nrow(sub)), call. = FALSE)
  }
  info
}

verify_shared_base_across_ssps <- function(base_by_ssp,
                                           date_col = "date",
                                           obs_col = "chirps",
                                           gcm_col = "gcm",
                                           atol = 1e-10) {
  ssps <- names(base_by_ssp)
  if (length(ssps) == 0L) stop("base_by_ssp is empty", call. = FALSE)
  ref_ssp <- ssps[1]
  ref <- base_by_ssp[[ref_ssp]][, c(date_col, obs_col, gcm_col), drop = FALSE]
  rownames(ref) <- NULL

  rows <- vector("list", length(ssps))
  for (i in seq_along(ssps)) {
    ssp <- ssps[i]
    cur <- base_by_ssp[[ssp]][, c(date_col, obs_col, gcm_col), drop = FALSE]
    rownames(cur) <- NULL
    same_n <- nrow(cur) == nrow(ref)
    same_date <- same_n && identical(cur[[date_col]], ref[[date_col]])
    if (same_n && same_date) {
      obs_a <- cur[[obs_col]]; obs_b <- ref[[obs_col]]
      gcm_a <- cur[[gcm_col]]; gcm_b <- ref[[gcm_col]]
      obs_missing_same <- identical(is.na(obs_a), is.na(obs_b))
      gcm_missing_same <- identical(is.na(gcm_a), is.na(gcm_b))
      d_obs <- abs(obs_a - obs_b)
      d_gcm <- abs(gcm_a - gcm_b)
      obs_ok <- obs_missing_same && all(d_obs[!is.na(d_obs)] <= atol)
      gcm_ok <- gcm_missing_same && all(d_gcm[!is.na(d_gcm)] <= atol)
      max_obs <- if (all(is.na(d_obs))) 0 else max(d_obs, na.rm = TRUE)
      max_gcm <- if (all(is.na(d_gcm))) 0 else max(d_gcm, na.rm = TRUE)
    } else {
      obs_ok <- FALSE; gcm_ok <- FALSE; max_obs <- NA_real_; max_gcm <- NA_real_
    }
    rows[[i]] <- data.frame(
      reference_ssp = ref_ssp,
      ssp = ssp,
      same_n_rows = same_n,
      same_dates = same_date,
      same_chirps_base = obs_ok,
      same_gcm_historical_base = gcm_ok,
      max_abs_chirps_diff = max_obs,
      max_abs_gcm_diff = max_gcm,
      stringsAsFactors = FALSE
    )
  }
  check <- do.call(rbind, rows)
  if (!all(check$same_dates) || !all(check$same_chirps_base) || !all(check$same_gcm_historical_base)) {
    stop(paste0(
      "The 1990-2014 base period is not identical across SSP files. ",
      "Do not share one Stage-1 fit until the historical data are reconciled.\n",
      paste(capture.output(print(check, row.names = FALSE)), collapse = "\n")
    ), call. = FALSE)
  }
  check
}

# =============================================================================
# MONTHLY RAINFALL FEATURES
# =============================================================================

monthly_rainfall_stats <- function(daily, value_col, taus, wet_threshold) {
  keys <- interaction(daily$year, daily$month, drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(daily)), keys)
  rows <- vector("list", length(groups))
  j <- 0L

  for (idx in groups) {
    g <- daily[idx, , drop = FALSE]
    v <- suppressWarnings(as.numeric(g[[value_col]]))
    v <- v[is.finite(v)]
    if (length(v) == 0L) next
    wet_values <- v[v > wet_threshold]
    j <- j + 1L
    row <- list(
      year = as.integer(g$year[1]),
      month = as.integer(g$month[1]),
      n_days = length(v),
      n_wet = length(wet_values),
      wet_prop = length(wet_values) / length(v),
      dry_prop = 1 - length(wet_values) / length(v)
    )
    for (tau in taus) {
      row[[.qname(tau)]] <- if (length(wet_values)) .quantile_safe(wet_values, tau) else NA_real_
    }
    rows[[j]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }
  rows <- rows[seq_len(j)]
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out)) {
    out <- out[order(out$year, out$month), , drop = FALSE]
    rownames(out) <- NULL
  }
  out
}

paired_monthly_stats <- function(daily, obs_col, gcm_col, taus, wet_threshold) {
  obs <- monthly_rainfall_stats(daily, obs_col, taus, wet_threshold)
  gcm <- monthly_rainfall_stats(daily, gcm_col, taus, wet_threshold)
  out <- merge(obs, gcm, by = c("year", "month"), all = FALSE, suffixes = c("_obs", "_gcm"), sort = TRUE)
  out <- out[order(out$year, out$month), , drop = FALSE]
  rownames(out) <- NULL
  out
}


# -----------------------------------------------------------------------------
# TARGET-ONLY FALLBACK FOR ZERO/SPARSE GCM WET-DAY MONTHS
# -----------------------------------------------------------------------------
# This safeguards Stage 1 when a target GCM month has no wet-day value at or
# above wet_threshold, which otherwise makes q*_gcm = NA and can later collapse
# the Stage-2 quantile-matching knots.
#
# The fallback uses GCM information only:
#   1) same calendar month from nearby TARGET years,
#   2) same calendar month from all TARGET years,
#   3) historical GCM same-month quantile climatology from fit$base_monthly.
# Target CHIRPS values are never used to fill GCM predictor quantiles.
.repair_sparse_target_gcm_quantiles <- function(target_daily,
                                                monthly,
                                                gcm_col,
                                                taus,
                                                wet_threshold,
                                                base_monthly = NULL,
                                                enabled = TRUE,
                                                min_wet_target_quantile = 1L,
                                                fallback_min_pool_wet = 10L,
                                                fallback_max_year_radius = 5L) {
  out <- monthly

  qcols <- vapply(taus, function(tt) paste0(.qname(tt), "_gcm"), character(1))
  raw_qcols <- paste0(qcols, "_raw")

  for (j in seq_along(qcols)) out[[raw_qcols[j]]] <- out[[qcols[j]]]

  out$q_fallback_gcm <- FALSE
  out$q_source_gcm <- "original_month"
  out$q_pool_radius_gcm <- 0L
  out$n_wet_pool_gcm <- as.integer(out$n_wet_gcm)

  if (!isTRUE(enabled)) return(out)

  min_wet_target_quantile <- max(1L, as.integer(min_wet_target_quantile))
  fallback_min_pool_wet <- max(1L, as.integer(fallback_min_pool_wet))
  fallback_max_year_radius <- max(1L, as.integer(fallback_max_year_radius))

  for (i in seq_len(nrow(out))) {
    current_q <- vapply(qcols, function(cc) as.numeric(out[[cc]][i]), numeric(1))
    needs_repair <- as.integer(out$n_wet_gcm[i]) < min_wet_target_quantile ||
      any(!is.finite(current_q))

    if (!needs_repair) next

    yy <- as.integer(out$year[i])
    mm <- as.integer(out$month[i])

    chosen <- numeric(0)
    chosen_source <- NA_character_
    chosen_radius <- NA_integer_

    # First choice: same calendar month in nearby target years.
    for (rad in seq_len(fallback_max_year_radius)) {
      keep <- target_daily$month == mm &
        abs(as.integer(target_daily$year) - yy) <= rad
      x <- suppressWarnings(as.numeric(target_daily[[gcm_col]][keep]))
      x <- x[is.finite(x) & x > wet_threshold]

      if (length(x) >= fallback_min_pool_wet) {
        chosen <- x
        chosen_source <- paste0("same_month_pool_pm", rad)
        chosen_radius <- rad
        break
      }
    }

    # Second choice: all target years, same calendar month.
    if (length(chosen) == 0L) {
      keep <- target_daily$month == mm
      x <- suppressWarnings(as.numeric(target_daily[[gcm_col]][keep]))
      x <- x[is.finite(x) & x > wet_threshold]
      if (length(x) > 0L) {
        chosen <- x
        chosen_source <- "same_month_all_target_years"
        chosen_radius <- NA_integer_
      }
    }

    if (length(chosen) > 0L) {
      qhat <- .quantile_safe(chosen, taus)
      for (j in seq_along(qcols)) out[[qcols[j]]][i] <- qhat[j]
      out$q_fallback_gcm[i] <- TRUE
      out$q_source_gcm[i] <- chosen_source
      out$q_pool_radius_gcm[i] <- chosen_radius
      out$n_wet_pool_gcm[i] <- length(chosen)
      next
    }

    # Last resort: historical GCM same-month quantile climatology.
    if (!is.null(base_monthly)) {
      bm <- base_monthly[base_monthly$month == mm, , drop = FALSE]
      qhat <- vapply(qcols, function(cc) {
        z <- suppressWarnings(as.numeric(bm[[cc]]))
        z <- z[is.finite(z)]
        if (length(z)) stats::median(z) else NA_real_
      }, numeric(1))

      if (all(is.finite(qhat))) {
        for (j in seq_along(qcols)) out[[qcols[j]]][i] <- qhat[j]
        out$q_fallback_gcm[i] <- TRUE
        out$q_source_gcm[i] <- "base_gcm_same_month_climatology"
        out$q_pool_radius_gcm[i] <- NA_integer_
        out$n_wet_pool_gcm[i] <- if ("n_wet_gcm" %in% names(bm)) {
          as.integer(sum(bm$n_wet_gcm, na.rm = TRUE))
        } else {
          NA_integer_
        }
        next
      }
    }

    stop(sprintf(
      paste0(
        "Unable to repair target GCM wet-day quantiles for %04d-%02d. ",
        "No GCM values >= %.6g were available in same-month target pooling, ",
        "and historical GCM climatology was unavailable."
      ),
      yy, mm, wet_threshold
    ), call. = FALSE)
  }

  out
}

# =============================================================================
# STAGE 1A - FIT SHARED BASE MODEL (1990-2014)
# =============================================================================

.get_model <- function(models, kind, month, tau = NULL) {
  mask <- models$kind == kind & models$month == month
  if (!is.null(tau)) mask <- mask & !is.na(models$tau) & abs(models$tau - tau) < 1e-12
  hit <- models[mask, , drop = FALSE]
  if (nrow(hit) != 1L) {
    stop(sprintf("Expected one model for kind=%s, month=%d, tau=%s; found %d",
                 kind, month, ifelse(is.null(tau), "NULL", as.character(tau)), nrow(hit)), call. = FALSE)
  }
  hit[1, , drop = FALSE]
}

.predict_quantile_model <- function(model, x, wet_threshold) {
  x <- as.numeric(x)
  xx <- if (model$x_transform[1] == "log") log(pmax(x, wet_threshold * 1e-6)) else x
  pred <- as.numeric(model$intercept[1]) + as.numeric(model$slope[1]) * xx
  if (model$y_transform[1] == "log") pred <- exp(pred)
  as.numeric(pred)
}

.predict_occurrence_model <- function(model, p) {
  xx <- .logit(p)
  pred <- as.numeric(model$intercept[1]) + as.numeric(model$slope[1]) * xx
  .inv_logit(pred)
}

.enforce_rowwise_monotonic_quantiles <- function(df, pred_cols, floor_value) {
  out <- df
  a <- as.matrix(out[, pred_cols, drop = FALSE])
  storage.mode(a) <- "double"
  crossing <- apply(a, 1, function(r) any(diff(r) < 0, na.rm = TRUE))
  n_crossing <- sum(crossing, na.rm = TRUE)
  for (i in seq_len(nrow(a))) {
    r <- pmax(a[i, ], floor_value)
    # Preserve NA while applying cumulative maximum to finite entries in order.
    if (anyNA(r)) {
      last <- -Inf
      for (k in seq_along(r)) {
        if (is.finite(r[k])) {
          last <- max(last, r[k])
          r[k] <- last
        }
      }
    } else {
      r <- cummax(r)
    }
    a[i, ] <- r
  }
  out[, pred_cols] <- a
  list(data = out, n_crossing = as.integer(n_crossing))
}

fit_stage1_base <- function(base_daily,
                            obs_col = "chirps",
                            gcm_col = "gcm",
                            taus = c(0.10, 0.50, 0.90, 0.95, 0.99),
                            wet_threshold = 0.1,
                            min_train_years = 8L,
                            use_log_rain = TRUE,
                            enforce_monotone = TRUE) {
  start <- proc.time()[["elapsed"]]
  monthly <- paired_monthly_stats(base_daily, obs_col, gcm_col, taus, wet_threshold)
  model_rows <- list()
  k <- 0L

  for (month in 1:12) {
    bm <- monthly[monthly$month == month, , drop = FALSE]

    fit_occ <- tryCatch(
      .fit_ols(.logit(bm$wet_prop_gcm), .logit(bm$wet_prop_obs), min_train_years),
      error = function(e) stop(sprintf("Wet-day occurrence model failed for month=%d: %s", month, conditionMessage(e)), call. = FALSE)
    )
    k <- k + 1L
    model_rows[[k]] <- data.frame(
      kind = "wet_prop", month = month, tau = NA_real_,
      intercept = fit_occ$intercept, slope = fit_occ$slope, n_train = fit_occ$n,
      x_transform = "logit", y_transform = "logit", stringsAsFactors = FALSE
    )

    for (tau in taus) {
      qn <- .qname(tau)
      xraw <- as.numeric(bm[[paste0(qn, "_gcm")]])
      yraw <- as.numeric(bm[[paste0(qn, "_obs")]])
      if (isTRUE(use_log_rain)) {
        x <- log(pmax(xraw, wet_threshold * 1e-6))
        y <- log(pmax(yraw, wet_threshold * 1e-6))
        xt <- "log"; yt <- "log"
      } else {
        x <- xraw; y <- yraw; xt <- "identity"; yt <- "identity"
      }
      fq <- tryCatch(
        .fit_ols(x, y, min_train_years),
        error = function(e) stop(sprintf("Quantile model failed month=%d, tau=%s: %s", month, tau, conditionMessage(e)), call. = FALSE)
      )
      k <- k + 1L
      model_rows[[k]] <- data.frame(
        kind = "quantile", month = month, tau = as.numeric(tau),
        intercept = fq$intercept, slope = fq$slope, n_train = fq$n,
        x_transform = xt, y_transform = yt, stringsAsFactors = FALSE
      )
    }
  }

  models <- do.call(rbind, model_rows)
  fitted <- monthly
  fitted$pred_wet_prop <- NA_real_
  for (tau in taus) fitted[[paste0("pred_", .qname(tau))]] <- NA_real_

  for (month in 1:12) {
    idx <- which(fitted$month == month)
    occ <- .get_model(models, "wet_prop", month)
    fitted$pred_wet_prop[idx] <- .predict_occurrence_model(occ, fitted$wet_prop_gcm[idx])
    for (tau in taus) {
      qn <- .qname(tau)
      mod <- .get_model(models, "quantile", month, tau)
      fitted[[paste0("pred_", qn)]][idx] <- .predict_quantile_model(
        mod, fitted[[paste0(qn, "_gcm")]][idx], wet_threshold
      )
    }
  }

  fitted$pred_wet_prop <- pmin(pmax(fitted$pred_wet_prop, 1e-6), 1 - 1e-6)
  fitted$pred_dry_prop <- 1 - fitted$pred_wet_prop
  pred_cols <- vapply(taus, function(t) paste0("pred_", .qname(t)), character(1))
  if (isTRUE(enforce_monotone)) {
    tmp <- .enforce_rowwise_monotonic_quantiles(fitted, pred_cols, wet_threshold)
    fitted <- tmp$data; n_crossing <- tmp$n_crossing
  } else {
    a <- as.matrix(fitted[, pred_cols, drop = FALSE])
    n_crossing <- sum(apply(a, 1, function(r) any(diff(r) < 0, na.rm = TRUE)), na.rm = TRUE)
  }

  diagnostics <- data.frame(
    n_base_months = nrow(monthly),
    n_quantile_models = sum(models$kind == "quantile"),
    n_occurrence_models = sum(models$kind == "wet_prop"),
    n_base_months_with_quantile_crossing_before_adjustment = n_crossing,
    monotone_safeguard_applied = isTRUE(enforce_monotone),
    stringsAsFactors = FALSE
  )

  list(
    models = models,
    base_monthly = monthly,
    base_fitted = fitted,
    diagnostics = diagnostics,
    fit_seconds = .elapsed(start)
  )
}

# =============================================================================
# STAGE 1B - APPLY TO ONE SSP TARGET (2015-2025) + EVALUATE
# =============================================================================

.daily_pinball_for_conditional_wet_quantile <- function(target_daily,
                                                         projected_monthly,
                                                         obs_col,
                                                         tau,
                                                         wet_threshold) {
  qn <- .qname(tau)
  look <- projected_monthly[, c("year", "month", paste0("pred_", qn)), drop = FALSE]
  d <- merge(target_daily[, c("year", "month", obs_col), drop = FALSE],
             look, by = c("year", "month"), all.x = TRUE, sort = FALSE)
  d <- d[d[[obs_col]] > wet_threshold & !is.na(d[[obs_col]]), , drop = FALSE]
  list(loss = .pinball(d[[obs_col]], d[[paste0("pred_", qn)]], tau), n = nrow(d))
}

project_stage1_target <- function(fit,
                                  target_daily,
                                  obs_col = "chirps",
                                  gcm_col = "gcm",
                                  taus = c(0.10, 0.50, 0.90, 0.95, 0.99),
                                  wet_threshold = 0.1,
                                  enforce_monotone = TRUE,
                                  target_quantile_fallback = TRUE,
                                  min_wet_target_quantile = 1L,
                                  fallback_min_pool_wet = 10L,
                                  fallback_max_year_radius = 5L) {
  start <- proc.time()[["elapsed"]]

  # Raw target statistics. q*_gcm are NA when a target month has no wet GCM day.
  target_raw <- paired_monthly_stats(target_daily, obs_col, gcm_col, taus, wet_threshold)

  # Repair only the GCM predictor quantiles. Target CHIRPS is never used here.
  target <- .repair_sparse_target_gcm_quantiles(
    target_daily = target_daily,
    monthly = target_raw,
    gcm_col = gcm_col,
    taus = taus,
    wet_threshold = wet_threshold,
    base_monthly = fit$base_monthly,
    enabled = target_quantile_fallback,
    min_wet_target_quantile = min_wet_target_quantile,
    fallback_min_pool_wet = fallback_min_pool_wet,
    fallback_max_year_radius = fallback_max_year_radius
  )

  out <- target
  out$pred_wet_prop <- NA_real_
  for (tau in taus) out[[paste0("pred_", .qname(tau))]] <- NA_real_

  for (month in 1:12) {
    idx <- which(out$month == month)
    occ <- .get_model(fit$models, "wet_prop", month)
    out$pred_wet_prop[idx] <- .predict_occurrence_model(occ, out$wet_prop_gcm[idx])
    for (tau in taus) {
      qn <- .qname(tau)
      mod <- .get_model(fit$models, "quantile", month, tau)
      out[[paste0("pred_", qn)]][idx] <- .predict_quantile_model(
        mod, out[[paste0(qn, "_gcm")]][idx], wet_threshold
      )
    }
  }

  out$pred_wet_prop <- pmin(pmax(out$pred_wet_prop, 1e-6), 1 - 1e-6)
  out$pred_dry_prop <- 1 - out$pred_wet_prop
  pred_cols <- vapply(taus, function(tt) paste0("pred_", .qname(tt)), character(1))

  if (isTRUE(enforce_monotone)) {
    tmp <- .enforce_rowwise_monotonic_quantiles(out, pred_cols, wet_threshold)
    out <- tmp$data; n_crossing <- tmp$n_crossing
  } else {
    a <- as.matrix(out[, pred_cols, drop = FALSE])
    n_crossing <- sum(apply(a, 1, function(r) any(diff(r) < 0, na.rm = TRUE)), na.rm = TRUE)
  }

  # Fail here, with exact months, rather than later as a cryptic knot error.
  pred_mat <- as.matrix(out[, pred_cols, drop = FALSE])
  storage.mode(pred_mat) <- "double"
  bad <- !is.finite(out$pred_wet_prop) |
    apply(pred_mat, 1, function(r) any(!is.finite(r)))

  if (any(bad)) {
    keys <- sprintf("%04d-%02d", out$year[bad], out$month[bad])
    stop(sprintf(
      "Stage-1 target projection still contains non-finite values for: %s",
      paste(keys, collapse = ", ")
    ), call. = FALSE)
  }

  qmetrics <- vector("list", length(taus))
  for (i in seq_along(taus)) {
    tau <- taus[i]
    qn <- .qname(tau)
    yq <- out[[paste0(qn, "_obs")]]
    pq <- out[[paste0("pred_", qn)]]
    pin <- .daily_pinball_for_conditional_wet_quantile(target_daily, out, obs_col, tau, wet_threshold)
    qmetrics[[i]] <- data.frame(
      tau = tau,
      n_target_months = sum(is.finite(yq) & is.finite(pq)),
      n_wet_daily_for_pinball = pin$n,
      pinball_loss = pin$loss,
      MAE_monthly_quantile = .mae(yq, pq),
      RMSE_monthly_quantile = .rmse(yq, pq),
      Bias_monthly_quantile = .bias(yq, pq),
      stringsAsFactors = FALSE
    )
  }
  qmetrics <- do.call(rbind, qmetrics)

  ometrics <- data.frame(
    n_target_months = nrow(out),
    MAE_wet_prop = .mae(out$wet_prop_obs, out$pred_wet_prop),
    RMSE_wet_prop = .rmse(out$wet_prop_obs, out$pred_wet_prop),
    Bias_wet_prop = .bias(out$wet_prop_obs, out$pred_wet_prop),
    MAE_dry_prop = .mae(out$dry_prop_obs, out$pred_dry_prop),
    RMSE_dry_prop = .rmse(out$dry_prop_obs, out$pred_dry_prop),
    stringsAsFactors = FALSE
  )

  diagnostics <- data.frame(
    n_target_months = nrow(out),
    n_target_zero_wet_gcm_months = sum(out$n_wet_gcm == 0L, na.rm = TRUE),
    n_target_months_quantile_fallback = sum(out$q_fallback_gcm, na.rm = TRUE),
    target_quantile_fallback_enabled = isTRUE(target_quantile_fallback),
    min_wet_target_quantile = as.integer(min_wet_target_quantile),
    fallback_min_pool_wet = as.integer(fallback_min_pool_wet),
    fallback_max_year_radius = as.integer(fallback_max_year_radius),
    n_target_months_with_quantile_crossing_before_adjustment = n_crossing,
    monotone_safeguard_applied = isTRUE(enforce_monotone),
    stringsAsFactors = FALSE
  )

  qcols <- vapply(taus, function(tt) paste0(.qname(tt), "_gcm"), character(1))
  raw_qcols <- paste0(qcols, "_raw")
  fallback_cols <- unique(c(
    "year", "month",
    "n_days_gcm", "n_wet_gcm", "wet_prop_gcm", "dry_prop_gcm",
    "q_fallback_gcm", "q_source_gcm", "q_pool_radius_gcm", "n_wet_pool_gcm",
    raw_qcols, qcols, "pred_wet_prop", "pred_dry_prop", pred_cols
  ))
  fallback_cols <- fallback_cols[fallback_cols %in% names(out)]
  fallback_log <- out[out$q_fallback_gcm, fallback_cols, drop = FALSE]

  list(
    target_monthly_raw = target_raw,
    target_monthly = target,
    projected = out,
    fallback_log = fallback_log,
    metrics_quantiles = qmetrics,
    metrics_occurrence = ometrics,
    diagnostics = diagnostics,
    projection_seconds = .elapsed(start)
  )
}


# =============================================================================
# STAGE 2 - KOKIC RAINFALL QUANTILE MATCHING + ERROR CORRECTION
# =============================================================================

.collapse_knots <- function(xk, yk) {
  xk <- as.numeric(xk); yk <- as.numeric(yk)
  ok <- is.finite(xk) & is.finite(yk)
  xk <- xk[ok]; yk <- yk[ok]
  ord <- order(xk, seq_along(xk))
  xk <- xk[ord]; yk <- yk[ord]
  xs <- numeric(0); ys <- numeric(0)
  for (i in seq_along(xk)) {
    x <- xk[i]; y <- yk[i]
    if (length(xs) && abs(x - xs[length(xs)]) <= 1e-12) {
      ys[length(ys)] <- max(ys[length(ys)], y)
    } else {
      xs <- c(xs, x); ys <- c(ys, y)
    }
  }
  if (length(xs) < 2L) stop("Insufficient distinct quantile-matching knots", call. = FALSE)
  ys <- cummax(ys)
  list(x = xs, y = ys)
}

.piecewise_linear_with_45_tail <- function(x, xk, yk) {
  z <- .collapse_knots(xk, yk)
  .apply_prepared_map(x, z$x, z$y)
}

.prepare_map_knots <- function(xk, yk) .collapse_knots(xk, yk)

.apply_prepared_map <- function(x, xu, yu) {
  x <- as.numeric(x)
  out <- stats::approx(x = xu, y = yu, xout = x, method = "linear", rule = 2, ties = "ordered")$y
  hi <- x > xu[length(xu)]
  if (any(hi)) out[hi] <- yu[length(yu)] + (x[hi] - xu[length(xu)])
  lo <- x < xu[1]
  if (any(lo)) out[lo] <- pmax(0, yu[1] + (x[lo] - xu[1]))
  out
}

prepare_kokic_transform_plan <- function(base_obs_values,
                                         target_row,
                                         base_fit_row,
                                         taus,
                                         wet_threshold) {
  raw <- as.numeric(base_obs_values)
  if (any(!is.finite(raw))) stop("Historical source block contains missing/non-finite rainfall", call. = FALSE)
  empirical <- ifelse(raw <= wet_threshold, 0, raw)
  dry_mask <- empirical <= 0

  target_dry <- min(max(as.numeric(target_row$pred_dry_prop[1]), 1e-6), 1 - 1e-6)
  fit_dry <- min(max(as.numeric(base_fit_row$pred_dry_prop[1]), 1e-6), 1 - 1e-6)
  future_q <- vapply(taus, function(t) as.numeric(target_row[[paste0("pred_", .qname(t))]][1]), numeric(1))
  fitted_q <- vapply(taus, function(t) as.numeric(base_fit_row[[paste0("pred_", .qname(t))]][1]), numeric(1))
  future_q <- cummax(pmax(future_q, wet_threshold))
  fitted_q <- cummax(pmax(fitted_q, wet_threshold))

  p_future <- target_dry + (1 - target_dry) * taus
  p_fitted <- fit_dry + (1 - fit_dry) * taus

  x_future <- c(0, .quantile_safe(empirical, target_dry), .quantile_safe(empirical, p_future))
  x_fitted <- c(0, .quantile_safe(empirical, fit_dry), .quantile_safe(empirical, p_fitted))
  y_future_knots <- c(0, wet_threshold, future_q)
  y_fitted_knots <- c(0, wet_threshold, fitted_q)
  fk <- .prepare_map_knots(x_future, y_future_knots)
  bk <- .prepare_map_knots(x_fitted, y_fitted_knots)

  list(
    empirical = empirical,
    dry_mask = dry_mask,
    future_x = fk$x,
    future_y = fk$y,
    fitted_x = bk$x,
    fitted_y = bk$y,
    wet_threshold = wet_threshold
  )
}

apply_kokic_transform_plan <- function(plan) {
  y_work <- plan$empirical
  if (any(plan$dry_mask)) {
    y_work[plan$dry_mask] <- stats::runif(sum(plan$dry_mask), 0, plan$wet_threshold)
  }
  y_future <- .apply_prepared_map(y_work, plan$future_x, plan$future_y)
  y_fitted <- .apply_prepared_map(y_work, plan$fitted_x, plan$fitted_y)
  corrected <- y_future + y_work - y_fitted
  corrected <- ifelse(corrected > plan$wet_threshold, corrected, 0)
  pmax(corrected, 0)
}

kokic_transform_month <- function(base_obs_values,
                                  target_row,
                                  base_fit_row,
                                  taus,
                                  wet_threshold) {
  plan <- prepare_kokic_transform_plan(base_obs_values, target_row, base_fit_row, taus, wet_threshold)
  apply_kokic_transform_plan(plan)
}

.month_blocks <- function(daily) {
  keys <- interaction(daily$year, daily$month, drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(daily)), keys)
  out <- list()
  for (idx in groups) {
    g <- daily[idx, , drop = FALSE]
    g <- g[order(g$day), , drop = FALSE]
    key <- sprintf("%04d-%02d", g$year[1], g$month[1])
    out[[key]] <- g
  }
  out
}

.row_lookup <- function(monthly) {
  out <- list()
  for (i in seq_len(nrow(monthly))) {
    key <- sprintf("%04d-%02d", monthly$year[i], monthly$month[i])
    out[[key]] <- monthly[i, , drop = FALSE]
  }
  out
}

.eligible_source_blocks <- function(base_blocks, target_month, target_days, obs_col) {
  out <- list(); j <- 0L
  target_days <- as.integer(target_days)
  for (key in names(base_blocks)) {
    block <- base_blocks[[key]]
    month <- block$month[1]
    if (month != target_month) next
    if (any(is.na(block[[obs_col]]))) next
    days <- as.integer(block$day)
    if (length(days) == length(target_days) && identical(days, target_days)) {
      j <- j + 1L
      out[[j]] <- list(source_year = as.integer(block$year[1]), block = block)
    }
  }
  out
}

# Diagnose whether a historical month block can support the quantile-matching map.
# A valid piecewise-linear map requires at least two distinct empirical x-knots.
# With the Kokic dry-day treatment used here, this means that after applying
# wet_threshold there must be at least one positive/wet empirical value in
# addition to zero.
.source_block_knot_diagnostics <- function(values, wet_threshold) {
  raw <- suppressWarnings(as.numeric(values))
  finite <- is.finite(raw)
  empirical <- ifelse(finite & raw > wet_threshold, raw, 0)
  data.frame(
    n_days = length(raw),
    n_finite = sum(finite),
    n_wet = sum(empirical > 0, na.rm = TRUE),
    n_distinct_empirical = length(unique(empirical[finite])),
    n_distinct_x_support = length(unique(c(0, empirical[finite]))),
    min_empirical = if (any(finite)) min(empirical[finite], na.rm = TRUE) else NA_real_,
    max_empirical = if (any(finite)) max(empirical[finite], na.rm = TRUE) else NA_real_,
    qm_usable = any(finite) && any(empirical[finite] > 0),
    stringsAsFactors = FALSE
  )
}

# Public helper: run this on the 1990-2014 base_daily object to see which
# year-month blocks are unusable for Stage 2/3 quantile matching.
diagnose_qm_source_blocks <- function(base_daily,
                                      obs_col = "chirps",
                                      wet_threshold = 0.1) {
  blocks <- .month_blocks(base_daily)
  rows <- vector("list", length(blocks))
  i <- 0L
  for (key in names(blocks)) {
    b <- blocks[[key]]
    d <- .source_block_knot_diagnostics(b[[obs_col]], wet_threshold)
    i <- i + 1L
    rows[[i]] <- data.frame(
      source_key = key,
      source_year = as.integer(b$year[1]),
      month = as.integer(b$month[1]),
      d,
      stringsAsFactors = FALSE
    )
  }
  out <- .safe_rbind(rows)
  if (nrow(out)) {
    out <- out[order(out$month, out$source_year), , drop = FALSE]
    rownames(out) <- NULL
  }
  out
}

# Safe wrapper used by Stage 2 and Stage 3. It does NOT invent artificial
# knots. If a source block is degenerate, it is skipped and recorded.
.try_prepare_kokic_transform_plan <- function(base_obs_values,
                                              target_row,
                                              base_fit_row,
                                              taus,
                                              wet_threshold) {
  diag <- .source_block_knot_diagnostics(base_obs_values, wet_threshold)

  if (!all(is.finite(as.numeric(base_obs_values)))) {
    return(list(
      ok = FALSE,
      plan = NULL,
      reason = "Historical source block contains missing/non-finite rainfall",
      diagnostics = diag
    ))
  }

  if (!isTRUE(diag$qm_usable[1])) {
    return(list(
      ok = FALSE,
      plan = NULL,
      reason = sprintf(
        "Degenerate empirical rainfall block: n_wet=%d, n_distinct_x_support=%d",
        diag$n_wet[1], diag$n_distinct_x_support[1]
      ),
      diagnostics = diag
    ))
  }

  ans <- tryCatch(
    prepare_kokic_transform_plan(
      base_obs_values = base_obs_values,
      target_row = target_row,
      base_fit_row = base_fit_row,
      taus = taus,
      wet_threshold = wet_threshold
    ),
    error = function(e) e
  )

  if (inherits(ans, "error")) {
    return(list(
      ok = FALSE,
      plan = NULL,
      reason = conditionMessage(ans),
      diagnostics = diag
    ))
  }

  list(ok = TRUE, plan = ans, reason = NA_character_, diagnostics = diag)
}


.evaluate_stage2 <- function(point_daily, upper_taus = c(0.95, 0.99)) {
  valid <- complete.cases(point_daily[, c("obs", "qmb_stage2")])
  overall <- list(
    n_days = sum(valid),
    MAE = .mae(point_daily$obs, point_daily$qmb_stage2),
    RMSE = .rmse(point_daily$obs, point_daily$qmb_stage2),
    Bias = .bias(point_daily$obs, point_daily$qmb_stage2)
  )

  d <- point_daily
  d$year <- as.integer(format(d$date, "%Y"))
  d$month <- as.integer(format(d$date, "%m"))
  keys <- interaction(d$year, d$month, drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(d)), keys)
  rows <- list(); j <- 0L
  for (idx in groups) {
    g <- d[idx, , drop = FALSE]
    row <- list(year = g$year[1], month = g$month[1])
    for (tau in upper_taus) {
      lab <- paste0("P", round(100 * tau))
      qo <- if (any(!is.na(g$obs))) .quantile_safe(g$obs, tau) else NA_real_
      qp <- if (any(!is.na(g$qmb_stage2))) .quantile_safe(g$qmb_stage2, tau) else NA_real_
      row[[paste0(lab, "_obs")]] <- qo
      row[[paste0(lab, "_pred")]] <- qp
      row[[paste0(lab, "_error")]] <- qp - qo
      row[[paste0(lab, "_abs_error")]] <- abs(qp - qo)
    }
    j <- j + 1L
    rows[[j]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }
  upper <- do.call(rbind, rows)

  for (tau in upper_taus) {
    lab <- paste0("P", round(100 * tau))
    err <- upper[[paste0(lab, "_error")]]
    aerr <- upper[[paste0(lab, "_abs_error")]]
    overall[[paste0(lab, "_MAE")]] <- mean(aerr, na.rm = TRUE)
    overall[[paste0(lab, "_Bias")]] <- mean(err, na.rm = TRUE)
    overall[[paste0(lab, "_RMSE")]] <- sqrt(mean(err^2, na.rm = TRUE))
  }

  list(metrics = as.data.frame(overall, stringsAsFactors = FALSE), upper = upper)
}

run_stage2 <- function(base_daily,
                       target_daily,
                       fit,
                       projection,
                       obs_col = "chirps",
                       taus = c(0.10, 0.50, 0.90, 0.95, 0.99),
                       wet_threshold = 0.1,
                       random_state = 42L,
                       keep_candidate_blocks = FALSE) {
  start <- proc.time()[["elapsed"]]
  set.seed(as.integer(random_state))
  base_blocks <- .month_blocks(base_daily)
  target_blocks <- .month_blocks(target_daily)
  base_fit_lookup <- .row_lookup(fit$base_fitted)
  target_lookup <- .row_lookup(projection$projected)

  frames <- list(); f <- 0L
  plan_cache <- list()
  skipped <- list(); sk <- 0L
  n_plans_used <- 0L

  for (target_key in names(target_blocks)) {
    tg <- target_blocks[[target_key]]
    ty <- as.integer(tg$year[1]); tm <- as.integer(tg$month[1])
    if (is.null(target_lookup[[target_key]])) next

    sources <- .eligible_source_blocks(base_blocks, tm, tg$day, obs_col)
    if (length(sources) == 0L) {
      stop(sprintf("No calendar-compatible historical block for target %d-%02d", ty, tm), call. = FALSE)
    }

    plan_cache[[target_key]] <- list()

    for (src in sources) {
      sy <- src$source_year; sg <- src$block
      base_key <- sprintf("%04d-%02d", sy, tm)
      if (is.null(base_fit_lookup[[base_key]])) next

      attempt <- .try_prepare_kokic_transform_plan(
        base_obs_values = sg[[obs_col]],
        target_row = target_lookup[[target_key]],
        base_fit_row = base_fit_lookup[[base_key]],
        taus = taus,
        wet_threshold = wet_threshold
      )

      if (!isTRUE(attempt$ok)) {
        sk <- sk + 1L
        skipped[[sk]] <- data.frame(
          target_key = target_key,
          target_year = ty,
          month = tm,
          source_year = sy,
          reason = attempt$reason,
          attempt$diagnostics,
          stringsAsFactors = FALSE
        )
        next
      }

      plan <- attempt$plan
      plan_cache[[target_key]][[length(plan_cache[[target_key]]) + 1L]] <-
        list(source_year = sy, plan = plan)
      n_plans_used <- n_plans_used + 1L

      transformed <- apply_kokic_transform_plan(plan)
      f <- f + 1L
      frames[[f]] <- data.frame(
        target_date = tg$date,
        target_year = ty,
        month = tm,
        source_year = sy,
        projected = transformed,
        stringsAsFactors = FALSE
      )
    }

    # It is acceptable to skip a few degenerate historical blocks, but every
    # target month must retain at least one valid historical source block.
    if (length(plan_cache[[target_key]]) == 0L) {
      bad_years <- if (length(skipped)) {
        z <- .safe_rbind(skipped)
        z <- z[z$target_key == target_key, , drop = FALSE]
        paste(unique(z$source_year), collapse = ", ")
      } else "none recorded"

      stop(sprintf(
        paste0(
          "No valid quantile-matching source block for target %d-%02d. ",
          "%d calendar-compatible source blocks were checked, but all had degenerate/invalid knots. ",
          "Problematic source years: %s. Check base-period CHIRPS values and wet_threshold."
        ),
        ty, tm, length(sources), bad_years
      ), call. = FALSE)
    }
  }

  if (length(frames) == 0L) stop("Stage 2 generated no candidate blocks", call. = FALSE)
  candidates <- do.call(rbind, frames)
  skipped_df <- .safe_rbind(skipped)

  if (nrow(skipped_df) > 0L) {
    message(sprintf(
      "Stage 2: skipped %d degenerate source-block uses; %d valid transformation plans retained.",
      nrow(skipped_df), n_plans_used
    ))
  }

  point <- stats::aggregate(projected ~ target_date, data = candidates, FUN = median, na.rm = TRUE)
  names(point) <- c("date", "qmb_stage2")
  point$date <- as.Date(point$date, origin = "1970-01-01")
  obs <- target_daily[, c("date", obs_col), drop = FALSE]
  names(obs)[2] <- "obs"
  point <- merge(obs, point, by = "date", all.x = TRUE, sort = FALSE)
  point <- point[order(point$date), , drop = FALSE]
  rownames(point) <- NULL
  ev <- .evaluate_stage2(point)

  diagnostics <- data.frame(
    n_target_months = length(plan_cache),
    n_valid_transform_plans = n_plans_used,
    n_skipped_source_block_uses = nrow(skipped_df),
    n_unique_skipped_source_year_month = if (nrow(skipped_df))
      length(unique(sprintf("%04d-%02d", skipped_df$source_year, skipped_df$month))) else 0L,
    stringsAsFactors = FALSE
  )

  list(
    point_daily = point,
    metrics = ev$metrics,
    monthly_upper_quantiles = ev$upper,
    candidate_blocks = if (isTRUE(keep_candidate_blocks)) candidates else NULL,
    plans = plan_cache,
    diagnostics = diagnostics,
    skipped_source_blocks = skipped_df,
    seconds = .elapsed(start)
  )
}


# =============================================================================
# STAGE 3 - MONTH-BLOCK BOOTSTRAP
# =============================================================================

run_stage3 <- function(base_daily,
                       target_daily,
                       fit,
                       projection,
                       obs_col = "chirps",
                       taus = c(0.10, 0.50, 0.90, 0.95, 0.99),
                       wet_threshold = 0.1,
                       n_boot = 1000L,
                       interval = 0.95,
                       random_state = 123L,
                       keep_ensemble = FALSE,
                       precomputed_plans = NULL) {
  if (n_boot < 2L) stop("n_boot must be >= 2", call. = FALSE)
  if (!(interval > 0 && interval < 1)) stop("interval must be in (0,1)", call. = FALSE)

  start <- proc.time()[["elapsed"]]
  set.seed(as.integer(random_state))
  base_blocks <- .month_blocks(base_daily)
  target_daily <- target_daily[order(target_daily$date), , drop = FALSE]
  rownames(target_daily) <- NULL
  base_fit_lookup <- .row_lookup(fit$base_fitted)
  target_lookup <- .row_lookup(projection$projected)

  tkeys <- interaction(target_daily$year, target_daily$month, drop = TRUE, lex.order = TRUE)
  tgroups <- split(seq_len(nrow(target_daily)), tkeys)
  groups <- list(); gi <- 0L

  for (positions in tgroups) {
    tg <- target_daily[positions, , drop = FALSE]
    ty <- as.integer(tg$year[1]); tm <- as.integer(tg$month[1])
    key <- sprintf("%04d-%02d", ty, tm)
    if (is.null(target_lookup[[key]])) stop(sprintf("Missing Stage-1 target projection for %s", key), call. = FALSE)

    if (!is.null(precomputed_plans) && !is.null(precomputed_plans[[key]])) {
      plans <- precomputed_plans[[key]]
    } else {
      sources <- .eligible_source_blocks(base_blocks, tm, tg$day, obs_col)
      plans <- list()
      for (src in sources) {
        sy <- src$source_year; sg <- src$block
        base_key <- sprintf("%04d-%02d", sy, tm)
        if (is.null(base_fit_lookup[[base_key]])) next

        attempt <- .try_prepare_kokic_transform_plan(
          base_obs_values = sg[[obs_col]],
          target_row = target_lookup[[key]],
          base_fit_row = base_fit_lookup[[base_key]],
          taus = taus,
          wet_threshold = wet_threshold
        )
        if (!isTRUE(attempt$ok)) next

        plans[[length(plans) + 1L]] <- list(source_year = sy, plan = attempt$plan)
      }
    }
    if (length(plans) == 0L) {
      stop(sprintf(
        paste0(
          "No valid bootstrap source block for target %d-%02d after knot validation. ",
          "Check base-period CHIRPS values and wet_threshold; Stage 3 cannot bootstrap an undefined quantile map."
        ),
        ty, tm
      ), call. = FALSE)
    }
    gi <- gi + 1L
    groups[[gi]] <- list(key = key, positions = positions, plans = plans)
  }

  ens <- matrix(NA_real_, nrow = n_boot, ncol = nrow(target_daily))
  for (b in seq_len(n_boot)) {
    for (grp in groups) {
      pick <- sample.int(length(grp$plans), 1L)
      transformed <- apply_kokic_transform_plan(grp$plans[[pick]]$plan)
      ens[b, grp$positions] <- transformed
    }
  }

  alpha <- 1 - interval
  qcol <- function(p) apply(ens, 2, stats::quantile, probs = p, na.rm = TRUE, names = FALSE, type = 7)
  lower <- qcol(alpha / 2)
  med <- qcol(0.5)
  upper <- qcol(1 - alpha / 2)

  di <- data.frame(
    date = target_daily$date,
    obs = as.numeric(target_daily[[obs_col]]),
    lower = lower,
    median = med,
    upper = upper,
    stringsAsFactors = FALSE
  )
  valid <- complete.cases(di[, c("obs", "lower", "upper")])
  inside <- di$obs[valid] >= di$lower[valid] & di$obs[valid] <= di$upper[valid]
  coverage <- if (length(inside)) mean(inside) else NA_real_
  widths <- di$upper[valid] - di$lower[valid]
  width <- if (length(widths)) mean(widths) else NA_real_

  metrics <- data.frame(
    nominal_interval = interval,
    coverage = coverage,
    coverage_error_signed = if (is.finite(coverage)) coverage - interval else NA_real_,
    coverage_error_abs = if (is.finite(coverage)) abs(coverage - interval) else NA_real_,
    mean_interval_width = width,
    median_interval_width = if (length(widths)) median(widths) else NA_real_,
    n_days = sum(valid),
    n_boot = as.integer(n_boot),
    stringsAsFactors = FALSE
  )

  list(
    daily_interval = di,
    metrics = metrics,
    ensemble = if (isTRUE(keep_ensemble)) ens else NULL,
    seconds = .elapsed(start)
  )
}

# =============================================================================
# ONE SSP PIPELINE (uses a SHARED Stage-1 fit)
# =============================================================================

.flatten_ssp_summary <- function(region, model, ssp, fit, s1, s2, s3, n_ssp_sharing_fit) {
  row <- list(
    region = region,
    model = model,
    ssp = ssp,
    T1_fit_shared_sec = fit$fit_seconds,
    T1_projection_sec = s1$projection_seconds,
    T1_effective_amortized_sec = fit$fit_seconds / n_ssp_sharing_fit + s1$projection_seconds,
    T1_if_independent_sec = fit$fit_seconds + s1$projection_seconds,
    T2_sec = s2$seconds,
    T3_sec = s3$seconds,
    T_total_effective_amortized_sec = fit$fit_seconds / n_ssp_sharing_fit + s1$projection_seconds + s2$seconds + s3$seconds,
    T_total_if_independent_sec = fit$fit_seconds + s1$projection_seconds + s2$seconds + s3$seconds
  )

  for (i in seq_len(nrow(s1$metrics_quantiles))) {
    r <- s1$metrics_quantiles[i, , drop = FALSE]
    q <- paste0("Q", round(100 * r$tau))
    row[[paste0("S1_", q, "_Pinball")]] <- r$pinball_loss
    row[[paste0("S1_", q, "_MAE")]] <- r$MAE_monthly_quantile
    row[[paste0("S1_", q, "_RMSE")]] <- r$RMSE_monthly_quantile
    row[[paste0("S1_", q, "_Bias")]] <- r$Bias_monthly_quantile
  }

  om <- s1$metrics_occurrence[1, , drop = FALSE]
  for (cc in c("MAE_wet_prop", "RMSE_wet_prop", "Bias_wet_prop", "MAE_dry_prop", "RMSE_dry_prop")) {
    row[[paste0("S1_", cc)]] <- om[[cc]]
  }
  for (cc in names(s2$metrics)) row[[paste0("S2_", cc)]] <- s2$metrics[[cc]][1]
  for (cc in names(s3$metrics)) row[[paste0("S3_", cc)]] <- s3$metrics[[cc]][1]
  as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
}

run_one_ssp_with_shared_fit <- function(base_daily,
                                        target_daily,
                                        fit,
                                        region,
                                        model,
                                        ssp,
                                        obs_col = "chirps",
                                        gcm_col = "gcm",
                                        taus = c(0.10, 0.50, 0.90, 0.95, 0.99),
                                        wet_threshold = 0.1,
                                        target_quantile_fallback = TRUE,
                                        min_wet_target_quantile = 1L,
                                        fallback_min_pool_wet = 10L,
                                        fallback_max_year_radius = 5L,
                                        n_boot = 1000L,
                                        interval = 0.95,
                                        random_state = 42L,
                                        n_ssp_sharing_fit = 4L,
                                        keep_candidate_blocks = FALSE,
                                        keep_ensemble = FALSE) {
  s1 <- project_stage1_target(
    fit, target_daily, obs_col, gcm_col, taus, wet_threshold,
    target_quantile_fallback = target_quantile_fallback,
    min_wet_target_quantile = min_wet_target_quantile,
    fallback_min_pool_wet = fallback_min_pool_wet,
    fallback_max_year_radius = fallback_max_year_radius
  )
  s2 <- run_stage2(base_daily, target_daily, fit, s1, obs_col, taus, wet_threshold,
                   random_state = random_state + 100L,
                   keep_candidate_blocks = keep_candidate_blocks)
  s3 <- run_stage3(base_daily, target_daily, fit, s1, obs_col, taus, wet_threshold,
                   n_boot = n_boot, interval = interval, random_state = random_state + 200L,
                   keep_ensemble = keep_ensemble, precomputed_plans = s2$plans)

  runtime <- data.frame(
    component = c(
      "T1_fit_shared", "T1_projection", "T1_effective_amortized", "T1_if_independent",
      "T2_quantile_matching", "T3_block_bootstrap",
      "T_total_effective_amortized", "T_total_if_independent"
    ),
    seconds = c(
      fit$fit_seconds,
      s1$projection_seconds,
      fit$fit_seconds / n_ssp_sharing_fit + s1$projection_seconds,
      fit$fit_seconds + s1$projection_seconds,
      s2$seconds,
      s3$seconds,
      fit$fit_seconds / n_ssp_sharing_fit + s1$projection_seconds + s2$seconds + s3$seconds,
      fit$fit_seconds + s1$projection_seconds + s2$seconds + s3$seconds
    ),
    stringsAsFactors = FALSE
  )
  runtime$minutes <- runtime$seconds / 60
  summary <- .flatten_ssp_summary(region, model, ssp, fit, s1, s2, s3, n_ssp_sharing_fit)

  list(ssp = ssp, stage1 = s1, stage2 = s2, stage3 = s3, runtime = runtime, summary = summary)
}

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

save_ssp_result <- function(result, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result$stage1$projected, file.path(outdir, "stage1_projected_monthly_quantiles.csv"), row.names = FALSE)
  utils::write.csv(result$stage1$metrics_quantiles, file.path(outdir, "stage1_metrics_quantiles.csv"), row.names = FALSE)
  utils::write.csv(result$stage1$metrics_occurrence, file.path(outdir, "stage1_metrics_occurrence.csv"), row.names = FALSE)
  utils::write.csv(result$stage1$diagnostics, file.path(outdir, "stage1_projection_diagnostics.csv"), row.names = FALSE)
  if (!is.null(result$stage1$target_monthly_raw)) {
    utils::write.csv(result$stage1$target_monthly_raw, file.path(outdir, "stage1_target_monthly_raw.csv"), row.names = FALSE)
  }
  if (!is.null(result$stage1$fallback_log)) {
    utils::write.csv(result$stage1$fallback_log, file.path(outdir, "stage1_target_quantile_fallback_log.csv"), row.names = FALSE)
  }
  utils::write.csv(result$stage2$point_daily, file.path(outdir, "stage2_daily_representative.csv"), row.names = FALSE)
  utils::write.csv(result$stage2$metrics, file.path(outdir, "stage2_metrics.csv"), row.names = FALSE)
  utils::write.csv(result$stage2$monthly_upper_quantiles, file.path(outdir, "stage2_P95_P99_monthly.csv"), row.names = FALSE)
  if (!is.null(result$stage2$diagnostics)) {
    utils::write.csv(result$stage2$diagnostics, file.path(outdir, "stage2_knot_diagnostics_summary.csv"), row.names = FALSE)
  }
  if (!is.null(result$stage2$skipped_source_blocks) && nrow(result$stage2$skipped_source_blocks) > 0L) {
    utils::write.csv(result$stage2$skipped_source_blocks, file.path(outdir, "stage2_skipped_source_blocks.csv"), row.names = FALSE)
  }
  if (!is.null(result$stage2$candidate_blocks)) {
    utils::write.csv(result$stage2$candidate_blocks, file.path(outdir, "stage2_candidate_blocks.csv"), row.names = FALSE)
  }
  utils::write.csv(result$stage3$daily_interval, file.path(outdir, "stage3_daily_prediction_interval.csv"), row.names = FALSE)
  utils::write.csv(result$stage3$metrics, file.path(outdir, "stage3_metrics.csv"), row.names = FALSE)
  if (!is.null(result$stage3$ensemble)) {
    saveRDS(result$stage3$ensemble, file.path(outdir, "stage3_bootstrap_ensemble.rds"))
  }
  utils::write.csv(result$runtime, file.path(outdir, "runtime.csv"), row.names = FALSE)
  utils::write.csv(result$summary, file.path(outdir, "summary.csv"), row.names = FALSE)
}

save_shared_fit <- function(fit, outdir, base_check = NULL) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(fit$models, file.path(outdir, "stage1_shared_models.csv"), row.names = FALSE)
  utils::write.csv(fit$base_monthly, file.path(outdir, "stage1_base_monthly_quantiles.csv"), row.names = FALSE)
  utils::write.csv(fit$base_fitted, file.path(outdir, "stage1_base_fitted_quantiles.csv"), row.names = FALSE)
  utils::write.csv(fit$diagnostics, file.path(outdir, "stage1_fit_diagnostics.csv"), row.names = FALSE)
  utils::write.csv(data.frame(T1_fit_shared_sec = fit$fit_seconds,
                              T1_fit_shared_min = fit$fit_seconds / 60),
                   file.path(outdir, "stage1_fit_runtime.csv"), row.names = FALSE)
  if (!is.null(base_check)) {
    utils::write.csv(base_check, file.path(outdir, "base_consistency_across_ssp.csv"), row.names = FALSE)
  }
}

# =============================================================================
# MAIN WORKFLOW: 1 REGION + 1 GCM + 4 SSP FILES
# =============================================================================

run_one_region_one_gcm <- function(full_files,
                                   region,
                                   model,
                                   output_root,
                                   date_col = "date",
                                   obs_col = "chirps",
                                   gcm_col = "gcm",
                                   base_years = c(1990L, 2014L),
                                   target_years = c(2015L, 2025L),
                                   taus = c(0.10, 0.50, 0.90, 0.95, 0.99),
                                   wet_threshold = 0.1,
                                   target_quantile_fallback = TRUE,
                                   min_wet_target_quantile = 1L,
                                   fallback_min_pool_wet = 10L,
                                   fallback_max_year_radius = 5L,
                                   no_leap = TRUE,
                                   strict_calendar = FALSE,
                                   n_boot = 1000L,
                                   interval = 0.95,
                                   random_state = 42L,
                                   keep_candidate_blocks = FALSE,
                                   keep_ensemble = FALSE) {
  if (length(full_files) == 0L) stop("full_files is empty", call. = FALSE)
  if (is.null(names(full_files)) || any(names(full_files) == "")) {
    stop("full_files must be a named vector/list, e.g. c(ssp126='file.csv', ...)", call. = FALSE)
  }

  base_by_ssp <- list()
  target_by_ssp <- list()

  for (ssp in names(full_files)) {
    path <- as.character(full_files[[ssp]])
    d <- load_pair_csv(path, date_col, obs_col, gcm_col, no_leap)
    sp <- split_period(d, base_years, target_years)
    validate_period_coverage(sp$base, base_years, paste(region, model, ssp, "base", sep = "-"), no_leap, strict_calendar, obs_col)
    validate_period_coverage(sp$target, target_years, paste(region, model, ssp, "target", sep = "-"), no_leap, strict_calendar, obs_col)
    base_by_ssp[[ssp]] <- sp$base
    target_by_ssp[[ssp]] <- sp$target
  }

  base_check <- verify_shared_base_across_ssps(base_by_ssp, date_col, obs_col, gcm_col)
  reference_ssp <- names(full_files)[1]
  shared_base <- base_by_ssp[[reference_ssp]]

  fit <- fit_stage1_base(shared_base, obs_col, gcm_col, taus, wet_threshold)

  root <- file.path(output_root, .slug(region), .slug(model))
  shared_outdir <- file.path(root, "_shared_calibration_1990_2014")
  save_shared_fit(fit, shared_outdir, base_check)

  source_block_diag <- diagnose_qm_source_blocks(shared_base, obs_col, wet_threshold)
  utils::write.csv(
    source_block_diag,
    file.path(shared_outdir, "stage2_source_block_knot_diagnostics.csv"),
    row.names = FALSE
  )
  n_bad_source_blocks <- sum(!source_block_diag$qm_usable)
  if (n_bad_source_blocks > 0L) {
    message(sprintf(
      paste0(
        "Base-period knot check: %d historical year-month block(s) are degenerate after wet_threshold=%g. ",
        "They will be skipped in Stage 2/3 when alternative valid blocks exist."
      ),
      n_bad_source_blocks, wet_threshold
    ))
  }

  results <- list()
  summaries <- list()
  n_ssp <- length(full_files)
  i <- 0L
  for (ssp in names(target_by_ssp)) {
    i <- i + 1L
    res <- run_one_ssp_with_shared_fit(
      shared_base, target_by_ssp[[ssp]], fit,
      region = region, model = model, ssp = ssp,
      obs_col = obs_col, gcm_col = gcm_col, taus = taus,
      wet_threshold = wet_threshold,
      target_quantile_fallback = target_quantile_fallback,
      min_wet_target_quantile = min_wet_target_quantile,
      fallback_min_pool_wet = fallback_min_pool_wet,
      fallback_max_year_radius = fallback_max_year_radius,
      n_boot = n_boot, interval = interval,
      random_state = as.integer(random_state + (i - 1L) * 1000L),
      n_ssp_sharing_fit = n_ssp,
      keep_candidate_blocks = keep_candidate_blocks,
      keep_ensemble = keep_ensemble
    )
    results[[ssp]] <- res
    summaries[[i]] <- res$summary
    save_ssp_result(res, file.path(root, .slug(ssp)))
  }

  master <- do.call(rbind, summaries)
  rownames(master) <- NULL
  utils::write.csv(master, file.path(root, "QMB_summary_all_SSP.csv"), row.names = FALSE)

  list(
    region = region,
    model = model,
    fit = fit,
    ssp_results = results,
    base_consistency = base_check,
    master_summary = master
  )
}

# =============================================================================
# EASY-EDIT EXAMPLE
# =============================================================================
# In R / Google Colab with R runtime:
#
# source("qmb_kokic_single_region_gcm.R")
#
# FULL_FILES <- c(
#   ssp126 = "/content/drive/MyDrive/DISERTASI2/QMB_Kokic/DATA/Jawa_Timur/CESM2_WACCM/Jawa_Timur_CESM2_WACCM_ssp126_paired_1990_2025.csv",
#   ssp245 = "/content/drive/MyDrive/DISERTASI2/QMB_Kokic/DATA/Jawa_Timur/CESM2_WACCM/Jawa_Timur_CESM2_WACCM_ssp245_paired_1990_2025.csv",
#   ssp370 = "/content/drive/MyDrive/DISERTASI2/QMB_Kokic/DATA/Jawa_Timur/CESM2_WACCM/Jawa_Timur_CESM2_WACCM_ssp370_paired_1990_2025.csv",
#   ssp585 = "/content/drive/MyDrive/DISERTASI2/QMB_Kokic/DATA/Jawa_Timur/CESM2_WACCM/Jawa_Timur_CESM2_WACCM_ssp585_paired_1990_2025.csv"
# )
#
# TEST RUN:
# result_test <- run_one_region_one_gcm(
#   full_files = FULL_FILES,
#   region = "Jawa_Timur",
#   model = "CESM2_WACCM",
#   output_root = "/content/drive/MyDrive/DISERTASI2/QMB_Kokic/OUTPUT/TEST_R",
#   date_col = "date",
#   obs_col = "chirps",
#   gcm_col = "gcm",
#   base_years = c(1990, 2014),
#   target_years = c(2015, 2025),
#   taus = c(0.10, 0.50, 0.90, 0.95, 0.99),
#   wet_threshold = 0.1,
#   no_leap = TRUE,
#   strict_calendar = FALSE,
#   n_boot = 20,
#   interval = 0.95,
#   random_state = 42
# )
#
# print(result_test$master_summary)
# print(result_test$ssp_results$ssp126$stage1$metrics_quantiles)
# print(result_test$ssp_results$ssp126$stage2$metrics)
# print(result_test$ssp_results$ssp126$stage3$metrics)
# print(result_test$ssp_results$ssp126$runtime)
#
# FINAL RUN: change n_boot = 1000 and strict_calendar = TRUE.
# =============================================================================

# =============================================================================
# FINAL Q1-AUDIT OVERRIDES - 2026-08-25
# =============================================================================
# Key corrections:
# 1) Wet day is consistently defined as rainfall > wet_threshold.
# 2) Kokic Stage-1 defaults to ORIGINAL-SCALE linear regression (use_log_rain=FALSE).
# 3) Sparse BASE GCM quantiles are repaired only for fitted prediction/Stage-2,
#    never for Stage-1 regression training.
# 4) Regression equations are exportable as CSV.
# 5) Stage-2 extreme metrics use pooled calendar-month P95/P99 (2015-2025),
#    while legacy year-month diagnostics are retained.
# 6) Stage-3 adds the proper interval score.
# 7) Stage-2 fingerprint helper supports automatic B-invariance assertions.
# =============================================================================

QMB_KOKIC_CORE_VERSION <- "v3.0_FINAL_AUDITED_20260825"

# ----------------------------------------------------------------------------
# BASE GCM sparse-month fallback for prediction only (no training leakage)
# ----------------------------------------------------------------------------
.repair_sparse_base_gcm_quantiles <- function(
    base_daily,
    monthly,
    gcm_col,
    taus,
    wet_threshold,
    enabled = TRUE,
    min_wet_base_quantile = 1L,
    fallback_min_pool_wet = 10L,
    fallback_max_year_radius = 5L) {

  out <- monthly
  qcols <- vapply(taus, function(t) paste0(.qname(t), "_gcm"), character(1))
  raw_qcols <- paste0(qcols, "_raw")
  for (j in seq_along(qcols)) out[[raw_qcols[j]]] <- out[[qcols[j]]]

  out$q_fallback_base_gcm <- FALSE
  out$q_source_base_gcm <- "original_month"
  out$q_pool_radius_base_gcm <- 0L
  out$n_wet_pool_base_gcm <- as.integer(out$n_wet_gcm)

  if (!isTRUE(enabled)) return(out)

  min_wet_base_quantile <- max(1L, as.integer(min_wet_base_quantile))
  fallback_min_pool_wet <- max(1L, as.integer(fallback_min_pool_wet))
  fallback_max_year_radius <- max(1L, as.integer(fallback_max_year_radius))

  for (i in seq_len(nrow(out))) {
    current_q <- vapply(qcols, function(cc) as.numeric(out[[cc]][i]), numeric(1))
    nwet_i <- suppressWarnings(as.integer(out$n_wet_gcm[i]))
    needs_repair <- !is.finite(nwet_i) ||
      nwet_i < min_wet_base_quantile || any(!is.finite(current_q))
    if (!needs_repair) next

    yy <- as.integer(out$year[i]); mm <- as.integer(out$month[i])
    chosen <- numeric(0); chosen_source <- NA_character_; chosen_radius <- NA_integer_

    for (rad in seq_len(fallback_max_year_radius)) {
      keep <- base_daily$month == mm & abs(as.integer(base_daily$year) - yy) <= rad
      x <- suppressWarnings(as.numeric(base_daily[[gcm_col]][keep]))
      x <- x[is.finite(x) & x > wet_threshold]
      if (length(x) >= fallback_min_pool_wet) {
        chosen <- x
        chosen_source <- paste0("base_same_month_pool_pm", rad)
        chosen_radius <- rad
        break
      }
    }

    if (!length(chosen)) {
      keep <- base_daily$month == mm
      x <- suppressWarnings(as.numeric(base_daily[[gcm_col]][keep]))
      x <- x[is.finite(x) & x > wet_threshold]
      if (length(x)) {
        chosen <- x
        chosen_source <- "base_same_month_all_years"
        chosen_radius <- NA_integer_
      }
    }

    if (!length(chosen)) {
      stop(sprintf(
        "Unable to repair BASE GCM wet-day quantiles for %04d-%02d; no GCM values > %.6g in the same calendar month.",
        yy, mm, wet_threshold
      ), call. = FALSE)
    }

    qhat <- .quantile_safe(chosen, taus)
    if (any(!is.finite(qhat))) {
      stop(sprintf("BASE GCM fallback produced non-finite quantiles for %04d-%02d", yy, mm), call. = FALSE)
    }
    for (j in seq_along(qcols)) out[[qcols[j]]][i] <- qhat[j]
    out$q_fallback_base_gcm[i] <- TRUE
    out$q_source_base_gcm[i] <- chosen_source
    out$q_pool_radius_base_gcm[i] <- chosen_radius
    out$n_wet_pool_base_gcm[i] <- length(chosen)
  }
  out
}

# ----------------------------------------------------------------------------
# Stage-1 Kokic regression, revised for original-scale default and base fallback
# ----------------------------------------------------------------------------
fit_stage1_base <- function(base_daily,
                            obs_col = "chirps",
                            gcm_col = "gcm",
                            taus = c(0.10, 0.50, 0.90),
                            wet_threshold = 0.1,
                            min_train_years = 8L,
                            use_log_rain = FALSE,
                            enforce_monotone = TRUE,
                            base_quantile_fallback = TRUE,
                            min_wet_base_quantile = 1L,
                            base_fallback_min_pool_wet = 10L,
                            base_fallback_max_year_radius = 5L) {
  start <- proc.time()[["elapsed"]]
  monthly <- paired_monthly_stats(base_daily, obs_col, gcm_col, taus, wet_threshold)

  # IMPORTANT: monthly is the untouched regression-training table.
  # Repair is used only for historical fitted prediction required downstream.
  monthly_prediction <- .repair_sparse_base_gcm_quantiles(
    base_daily = base_daily,
    monthly = monthly,
    gcm_col = gcm_col,
    taus = taus,
    wet_threshold = wet_threshold,
    enabled = base_quantile_fallback,
    min_wet_base_quantile = min_wet_base_quantile,
    fallback_min_pool_wet = base_fallback_min_pool_wet,
    fallback_max_year_radius = base_fallback_max_year_radius
  )

  model_rows <- list(); k <- 0L
  for (month in 1:12) {
    bm <- monthly[monthly$month == month, , drop = FALSE]

    fit_occ <- tryCatch(
      .fit_ols(.logit(bm$wet_prop_gcm), .logit(bm$wet_prop_obs), min_train_years),
      error = function(e) stop(sprintf("Wet-day occurrence model failed for month=%d: %s", month, conditionMessage(e)), call. = FALSE)
    )
    k <- k + 1L
    model_rows[[k]] <- data.frame(
      kind = "wet_prop", month = month, tau = NA_real_,
      intercept = fit_occ$intercept, slope = fit_occ$slope, n_train = fit_occ$n,
      x_transform = "logit", y_transform = "logit", stringsAsFactors = FALSE
    )

    for (tau in taus) {
      qn <- .qname(tau)
      xraw <- as.numeric(bm[[paste0(qn, "_gcm")]])
      yraw <- as.numeric(bm[[paste0(qn, "_obs")]])
      if (isTRUE(use_log_rain)) {
        x <- log(pmax(xraw, wet_threshold * 1e-6))
        y <- log(pmax(yraw, wet_threshold * 1e-6))
        xt <- "log"; yt <- "log"
      } else {
        x <- xraw; y <- yraw; xt <- "identity"; yt <- "identity"
      }
      fq <- tryCatch(
        .fit_ols(x, y, min_train_years),
        error = function(e) stop(sprintf("Quantile model failed month=%d, tau=%s: %s", month, tau, conditionMessage(e)), call. = FALSE)
      )
      k <- k + 1L
      model_rows[[k]] <- data.frame(
        kind = "quantile", month = month, tau = as.numeric(tau),
        intercept = fq$intercept, slope = fq$slope, n_train = fq$n,
        x_transform = xt, y_transform = yt, stringsAsFactors = FALSE
      )
    }
  }

  models <- do.call(rbind, model_rows)
  fitted <- monthly_prediction
  fitted$pred_wet_prop <- NA_real_
  for (tau in taus) fitted[[paste0("pred_", .qname(tau))]] <- NA_real_

  for (month in 1:12) {
    idx <- which(fitted$month == month)
    occ <- .get_model(models, "wet_prop", month)
    fitted$pred_wet_prop[idx] <- .predict_occurrence_model(occ, fitted$wet_prop_gcm[idx])
    for (tau in taus) {
      qn <- .qname(tau)
      mod <- .get_model(models, "quantile", month, tau)
      fitted[[paste0("pred_", qn)]][idx] <- .predict_quantile_model(
        mod, fitted[[paste0(qn, "_gcm")]][idx], wet_threshold
      )
    }
  }

  fitted$pred_wet_prop <- pmin(pmax(fitted$pred_wet_prop, 1e-6), 1 - 1e-6)
  fitted$pred_dry_prop <- 1 - fitted$pred_wet_prop
  pred_cols <- vapply(taus, function(t) paste0("pred_", .qname(t)), character(1))
  if (isTRUE(enforce_monotone)) {
    tmp <- .enforce_rowwise_monotonic_quantiles(fitted, pred_cols, wet_threshold)
    fitted <- tmp$data; n_crossing <- tmp$n_crossing
  } else {
    a <- as.matrix(fitted[, pred_cols, drop = FALSE])
    n_crossing <- sum(apply(a, 1, function(r) any(diff(r) < 0, na.rm = TRUE)), na.rm = TRUE)
  }

  regression_equations <- .kokic_regression_equations(models)
  diagnostics <- data.frame(
    n_base_months = nrow(monthly),
    n_base_zero_wet_gcm_months = sum(monthly$n_wet_gcm == 0L, na.rm = TRUE),
    n_base_months_quantile_fallback = sum(monthly_prediction$q_fallback_base_gcm, na.rm = TRUE),
    base_quantile_fallback_enabled = isTRUE(base_quantile_fallback),
    n_quantile_models = sum(models$kind == "quantile"),
    n_occurrence_models = sum(models$kind == "wet_prop"),
    n_base_months_with_quantile_crossing_before_adjustment = n_crossing,
    monotone_safeguard_applied = isTRUE(enforce_monotone),
    quantile_regression_scale = if (isTRUE(use_log_rain)) "log-log" else "original-linear",
    wet_day_definition = paste0("rainfall > ", format(wet_threshold, scientific = FALSE)),
    stringsAsFactors = FALSE
  )

  list(
    method = "QMB_KOKIC",
    models = models,
    regression_equations = regression_equations,
    base_monthly = monthly,
    base_monthly_prediction = monthly_prediction,
    base_fallback_log = monthly_prediction[monthly_prediction$q_fallback_base_gcm, , drop = FALSE],
    base_fitted = fitted,
    diagnostics = diagnostics,
    fit_seconds = .elapsed(start),
    taus = as.numeric(taus),
    use_log_rain = isTRUE(use_log_rain)
  )
}

# ----------------------------------------------------------------------------
# Human-readable regression equations for paper/audit CSV
# ----------------------------------------------------------------------------
.kokic_regression_equations <- function(models) {
  if (is.null(models) || !nrow(models)) return(data.frame())
  rows <- vector("list", nrow(models))
  for (i in seq_len(nrow(models))) {
    r <- models[i, , drop = FALSE]
    a <- as.numeric(r$intercept[1]); b <- as.numeric(r$slope[1])
    if (r$kind[1] == "wet_prop") {
      eq_t <- sprintf("logit(p_obs) = %.12g + %.12g * logit(p_gcm)", a, b)
      eq_o <- sprintf("p_obs = inv_logit(%.12g + %.12g * logit(p_gcm))", a, b)
    } else if (r$x_transform[1] == "log" && r$y_transform[1] == "log") {
      eq_t <- sprintf("log(Q_obs) = %.12g + %.12g * log(Q_gcm)", a, b)
      eq_o <- sprintf("Q_obs = exp(%.12g) * Q_gcm^(%.12g)", a, b)
    } else {
      eq_t <- sprintf("Q_obs = %.12g + %.12g * Q_gcm", a, b)
      eq_o <- eq_t
    }
    rows[[i]] <- data.frame(
      method = "QMB_KOKIC",
      kind = as.character(r$kind[1]),
      month = as.integer(r$month[1]),
      tau = as.numeric(r$tau[1]),
      intercept = a,
      slope = b,
      n_train = as.integer(r$n_train[1]),
      x_transform = as.character(r$x_transform[1]),
      y_transform = as.character(r$y_transform[1]),
      equation_model_scale = eq_t,
      equation_original_scale = eq_o,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# ----------------------------------------------------------------------------
# Stage-2 evaluation: pooled calendar-month extremes are the primary P95/P99
# ----------------------------------------------------------------------------
.evaluate_stage2 <- function(point_daily, upper_taus = c(0.95, 0.99)) {
  valid <- complete.cases(point_daily[, c("obs", "qmb_stage2")])
  overall <- list(
    n_days = sum(valid),
    MAE = .mae(point_daily$obs, point_daily$qmb_stage2),
    RMSE = .rmse(point_daily$obs, point_daily$qmb_stage2),
    Bias = .bias(point_daily$obs, point_daily$qmb_stage2),
    extreme_metric_definition = "pooled_calendar_month_across_target_years"
  )

  d <- point_daily
  d$year <- as.integer(format(d$date, "%Y"))
  d$month <- as.integer(format(d$date, "%m"))

  # Legacy diagnostic by year-month, retained for traceability only.
  keys <- interaction(d$year, d$month, drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(d)), keys)
  legacy_rows <- list(); j <- 0L
  for (idx in groups) {
    g <- d[idx, , drop = FALSE]
    row <- list(year = g$year[1], month = g$month[1])
    for (tau in upper_taus) {
      lab <- paste0("P", round(100 * tau))
      qo <- .quantile_safe(g$obs, tau)
      qp <- .quantile_safe(g$qmb_stage2, tau)
      row[[paste0(lab, "_obs")]] <- qo
      row[[paste0(lab, "_pred")]] <- qp
      row[[paste0(lab, "_error")]] <- qp - qo
      row[[paste0(lab, "_abs_error")]] <- abs(qp - qo)
    }
    j <- j + 1L
    legacy_rows[[j]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }
  legacy <- if (length(legacy_rows)) do.call(rbind, legacy_rows) else data.frame()

  pooled_rows <- vector("list", 12L)
  for (m in 1:12) {
    g <- d[d$month == m, , drop = FALSE]
    row <- list(month = m, n_days = nrow(g))
    for (tau in upper_taus) {
      lab <- paste0("P", round(100 * tau))
      qo <- .quantile_safe(g$obs, tau)
      qp <- .quantile_safe(g$qmb_stage2, tau)
      row[[paste0(lab, "_obs")]] <- qo
      row[[paste0(lab, "_pred")]] <- qp
      row[[paste0(lab, "_error")]] <- qp - qo
      row[[paste0(lab, "_abs_error")]] <- abs(qp - qo)
    }
    pooled_rows[[m]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }
  pooled <- do.call(rbind, pooled_rows)

  for (tau in upper_taus) {
    lab <- paste0("P", round(100 * tau))
    err <- pooled[[paste0(lab, "_error")]]
    aerr <- pooled[[paste0(lab, "_abs_error")]]
    overall[[paste0(lab, "_MAE")]] <- mean(aerr, na.rm = TRUE)
    overall[[paste0(lab, "_Bias")]] <- mean(err, na.rm = TRUE)
    overall[[paste0(lab, "_RMSE")]] <- sqrt(mean(err^2, na.rm = TRUE))
  }

  list(
    metrics = as.data.frame(overall, stringsAsFactors = FALSE),
    upper = legacy,
    pooled_upper = pooled
  )
}

# Existing run_stage2() dynamically resolves .evaluate_stage2(). Wrap it only
# to expose the pooled table and a reproducibility fingerprint.
.run_stage2_pre_audit <- run_stage2
run_stage2 <- function(base_daily,
                       target_daily,
                       fit,
                       projection,
                       obs_col = "chirps",
                       taus = c(0.10, 0.50, 0.90),
                       wet_threshold = 0.1,
                       random_state = 42L,
                       keep_candidate_blocks = FALSE) {
  ans <- .run_stage2_pre_audit(
    base_daily = base_daily,
    target_daily = target_daily,
    fit = fit,
    projection = projection,
    obs_col = obs_col,
    taus = taus,
    wet_threshold = wet_threshold,
    random_state = random_state,
    keep_candidate_blocks = keep_candidate_blocks
  )
  ev <- .evaluate_stage2(ans$point_daily)
  ans$metrics <- ev$metrics
  ans$monthly_upper_quantiles <- ev$upper
  ans$pooled_upper_quantiles <- ev$pooled_upper
  ans$fingerprint <- .stage2_fingerprint(ans$point_daily)
  ans$stage2_seed <- as.integer(random_state)
  ans
}

.stage2_fingerprint <- function(point_daily) {
  if (is.null(point_daily) || !nrow(point_daily)) return(NA_character_)
  tf <- tempfile(fileext = ".txt")
  on.exit(unlink(tf), add = TRUE)
  x <- data.frame(
    date = format(as.Date(point_daily$date), "%Y-%m-%d"),
    qmb_stage2 = sprintf("%.12f", as.numeric(point_daily$qmb_stage2)),
    stringsAsFactors = FALSE
  )
  utils::write.table(x, tf, sep = "|", row.names = FALSE, col.names = TRUE, quote = FALSE)
  unname(tools::md5sum(tf))
}

# ----------------------------------------------------------------------------
# Stage-3 proper interval score (lower is better)
# ----------------------------------------------------------------------------
.interval_score <- function(y, lower, upper, interval = 0.95) {
  alpha <- 1 - interval
  y <- as.numeric(y); lower <- as.numeric(lower); upper <- as.numeric(upper)
  (upper - lower) +
    (2 / alpha) * (lower - y) * (y < lower) +
    (2 / alpha) * (y - upper) * (y > upper)
}

.run_stage3_pre_audit <- run_stage3
run_stage3 <- function(base_daily,
                       target_daily,
                       fit,
                       projection,
                       obs_col = "chirps",
                       taus = c(0.10, 0.50, 0.90),
                       wet_threshold = 0.1,
                       n_boot = 1000L,
                       interval = 0.95,
                       random_state = 123L,
                       keep_ensemble = FALSE,
                       precomputed_plans = NULL) {
  ans <- .run_stage3_pre_audit(
    base_daily = base_daily,
    target_daily = target_daily,
    fit = fit,
    projection = projection,
    obs_col = obs_col,
    taus = taus,
    wet_threshold = wet_threshold,
    n_boot = n_boot,
    interval = interval,
    random_state = random_state,
    keep_ensemble = keep_ensemble,
    precomputed_plans = precomputed_plans
  )
  di <- ans$daily_interval
  ok <- complete.cases(di[, c("obs", "lower", "upper")])
  scores <- .interval_score(di$obs[ok], di$lower[ok], di$upper[ok], interval)
  ans$metrics$mean_interval_score <- if (length(scores)) mean(scores, na.rm = TRUE) else NA_real_
  ans$metrics$median_interval_score <- if (length(scores)) stats::median(scores, na.rm = TRUE) else NA_real_
  ans$metrics$bootstrap_seed <- as.integer(random_state)
  ans
}

# ----------------------------------------------------------------------------
# Save wrappers: add equations, base fallback log, pooled extremes, fingerprint
# ----------------------------------------------------------------------------
.save_shared_fit_pre_audit <- save_shared_fit
save_shared_fit <- function(fit, outdir, base_check = NULL) {
  .save_shared_fit_pre_audit(fit, outdir, base_check)
  eq <- if (!is.null(fit$regression_equations)) fit$regression_equations else .kokic_regression_equations(fit$models)
  utils::write.csv(eq, file.path(outdir, "stage1_regression_equations.csv"), row.names = FALSE)
  if (!is.null(fit$base_monthly_prediction)) {
    utils::write.csv(fit$base_monthly_prediction,
                     file.path(outdir, "stage1_base_monthly_quantiles_prediction_repaired.csv"),
                     row.names = FALSE)
  }
  if (!is.null(fit$base_fallback_log)) {
    utils::write.csv(fit$base_fallback_log,
                     file.path(outdir, "stage1_base_quantile_fallback_log.csv"),
                     row.names = FALSE)
  }
  utils::write.csv(data.frame(core_version = QMB_KOKIC_CORE_VERSION,
                              taus = paste(fit$taus, collapse = "|"),
                              use_log_rain = fit$use_log_rain,
                              stringsAsFactors = FALSE),
                   file.path(outdir, "stage1_model_specification.csv"), row.names = FALSE)
  invisible(TRUE)
}

.save_ssp_result_pre_audit <- save_ssp_result
save_ssp_result <- function(result, outdir) {
  .save_ssp_result_pre_audit(result, outdir)
  if (!is.null(result$stage2$pooled_upper_quantiles)) {
    utils::write.csv(result$stage2$pooled_upper_quantiles,
                     file.path(outdir, "stage2_P95_P99_pooled_calendar_month.csv"),
                     row.names = FALSE)
  }
  utils::write.csv(data.frame(
    stage2_fingerprint = result$stage2$fingerprint,
    stage2_seed = result$stage2$stage2_seed,
    note = "Stage-2 is computed before Stage-3 and must be invariant to n_boot",
    stringsAsFactors = FALSE
  ), file.path(outdir, "stage2_reproducibility_fingerprint.csv"), row.names = FALSE)
  invisible(TRUE)
}

