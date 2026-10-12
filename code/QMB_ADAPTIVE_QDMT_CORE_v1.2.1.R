# =============================================================================
# QMB ADAPTIVE QDM-TAIL CUTOFF (5Q) -- PAPER / REPRODUCIBLE CORE
# Version: v1.2.1-audited (2026-09-28)
# =============================================================================
# PURPOSE
#   Extend the audited QMB Kokic 5Q implementation at Stage 2 AND Stage 3:
#     - Stage 2: QMB-LIN body below a data-selected cutoff tau*
#                full empirical QDM correction at / above tau*
#     - Stage 3: QMB month-block bootstrap is retained; each bootstrap member
#                is re-centered by the locked Stage-2 adaptive correction so
#                source-year uncertainty is preserved instead of collapsing
#                the upper-tail interval to deterministic QDM.
#
#   tau* is selected WITHOUT using 2015-2025 observations. Selection uses a
#   blocked pseudo-future cross-validation entirely inside 1990-2014.
#
# LOCKED SCIENTIFIC DESIGN
#   Calibration/source period : 1990-2014
#   Evaluation/target period  : 2015-2025
#   Stage-1 quantiles         : Q10, Q50, Q90, Q95, Q99
#   Wet-day threshold         : > 0.1 mm/day
#   QMB Stage 1               : audited Kokic core, ORIGINAL scale OLS
#   QMB Stage 2               : audited Kokic LIN + residual correction
#   QDM                       : pooled calendar-month empirical QDM
#   Candidate cutoffs         : 0.900, 0.905, ..., 0.950
#   CV                        : 5 blocked pseudo-future folds of 5 years
#   Primary selector          : robust all-fold non-inferiority to QDM at
#                               BOTH P95 and P99; default margin = 1.0 mm.
#   Selection principle       : choose the HIGHEST feasible cutoff, i.e. use
#                               QDM on as little of the distribution as possible.
#
# IMPORTANT
#   This file DOES NOT reimplement QMB Kokic. It sources the exact audited core:
#       QMB_KOKIC_CORE.R
#   Therefore Stage 1, dry-day treatment, same-calendar-month source blocks,
#   Stage-2 source-year median, 45-degree extrapolation, residual correction,
#   sparse-GCM fallback, and diagnostics remain those of the locked baseline.
#
# STAGE-3 ADAPTIVE DEFINITION
#   Let Y_LIN^(b)(t) be bootstrap member b from audited QMB-LIN Stage 3,
#       Y_LIN,2(t)   be the QMB-LIN Stage-2 representative, and
#       Y_ACT,2(t)   be the locked adaptive Stage-2 representative.
#
#   The adaptive bootstrap member is:
#       Y_ACT^(b)(t) = max[0,
#                         Y_LIN^(b)(t)
#                         + {Y_ACT,2(t) - Y_LIN,2(t)}]
#
#   Hence below the selected cutoff the ensemble is exactly QMB-LIN. At/above
#   the cutoff every QMB bootstrap member receives the same deterministic
#   Stage-2 tail correction. This preserves the source-year bootstrap anomaly
#   structure and avoids zero-width QDM tail intervals.
#
#   IMPORTANT: the Stage-2 correction is independent of B. Only the Monte Carlo
#   representation of Stage-3 uncertainty changes with n_boot.
#
# Dependencies: BASE R only + audited QMB core file.
# =============================================================================

QMB_ACT_VERSION <- "v1.1.0_STAGE3_20260907"
QMB_ACT_METHOD_ID <- "QMB_CV_ADAPTIVE_QDM_TAIL_5Q_STAGE3"
QMB_ACT_TAUS <- c(0.10, 0.50, 0.90, 0.95, 0.99)
QMB_ACT_CUTOFFS <- seq(0.90, 0.95, by = 0.005)

# -----------------------------------------------------------------------------
# 0. Load and validate the locked audited QMB core
# -----------------------------------------------------------------------------
.qmb_act_this_file <- function() {
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(f) || !nzchar(f)) return(NA_character_)
  normalizePath(f, winslash = "/", mustWork = FALSE)
}

.qmb_act_file <- .qmb_act_this_file()
.qmb_act_dir <- if (is.na(.qmb_act_file)) getwd() else dirname(.qmb_act_file)

load_qmb_audited_core <- function(core_file = file.path(.qmb_act_dir, "QMB_KOKIC_CORE.R")) {
  if (!file.exists(core_file)) {
    stop(
      paste0(
        "Audited QMB core not found:\n  ", core_file, "\n\n",
        "Place QMB_KOKIC_CORE.R in the same folder as this file, ",
        "or call load_qmb_audited_core('/full/path/to/core.R')."
      ),
      call. = FALSE
    )
  }
  source(core_file, local = .GlobalEnv)
  needed <- c(
    "prepare_daily_pair", "load_pair_csv", "split_period",
    "verify_shared_base_across_ssps", "fit_stage1_base",
    "project_stage1_target", "run_stage2", "run_stage3"
  )
  miss <- needed[!vapply(needed, exists, logical(1), mode = "function", inherits = TRUE)]
  if (length(miss)) {
    stop("Audited core loaded, but required functions are missing: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# 1. Small utilities kept local to this adaptive extension
# -----------------------------------------------------------------------------
.act_num <- function(x) suppressWarnings(as.numeric(x))

.act_q <- function(x, p) {
  x <- .act_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, type = 7, names = FALSE, na.rm = TRUE))
}

.act_mae <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat)
  if (!any(ok)) return(NA_real_)
  mean(abs(yhat[ok] - y[ok]))
}

.act_rmse <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((yhat[ok] - y[ok])^2))
}

.act_bias <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat)
  if (!any(ok)) return(NA_real_)
  mean(yhat[ok] - y[ok])
}

.act_safe_rbind <- function(x) {
  if (!length(x)) return(data.frame())
  do.call(rbind, x)
}

.act_norm_gcm <- function(x) {
  z <- toupper(as.character(x))
  z <- gsub("[- ]", "_", z)
  z <- gsub("__+", "_", z)
  z
}

.act_norm_ssp <- function(x) {
  z <- tolower(as.character(x))
  z <- gsub("[^0-9]", "", z)
  paste0("ssp", z)
}

# Default research grid routing used in the current dissertation pipeline.
# Override LOCATION_ID explicitly in a manifest if a future experiment uses
# truly independent location-1/location-2 files for every GCM.
default_location_id <- function(regime, gcm) {
  gg <- .act_norm_gcm(gcm)
  suffix <- ifelse(gg == "MRI_ESM2_0", 2L, 1L)
  paste(as.character(regime), suffix)
}

# -----------------------------------------------------------------------------
# 2. Empirical monthly QDM -- exactly the benchmark used in the final tests
# -----------------------------------------------------------------------------
# For each calendar month:
#   p_t = average_rank(GCM_target_t) / (n_target_month + 1)
#   Qg  = Q_GCM,base(p_t)
#   Qo  = Q_OBS,base(p_t)
#   if Qg > wet_threshold:
#       Y_QDM = Qo * (GCM_target / Qg)
#   else:
#       Y_QDM = Qo + (GCM_target - Qg)
#   then truncate at zero.
#
# Observations in target_daily are NOT used to calculate qdm_corrected; they are
# carried only for evaluation after prediction.
run_empirical_qdm <- function(base_daily,
                              target_daily,
                              obs_col = "chirps",
                              gcm_col = "gcm",
                              wet_threshold = 0.1) {
  b <- base_daily[order(base_daily$date), , drop = FALSE]
  t <- target_daily[order(target_daily$date), , drop = FALSE]

  out <- data.frame(
    date = t$date,
    year = t$year,
    month = t$month,
    day = t$day,
    obs = .act_num(t[[obs_col]]),
    gcm = .act_num(t[[gcm_col]]),
    p_future_qdm = NA_real_,
    q_gcm_hist = NA_real_,
    q_obs_hist = NA_real_,
    qdm_full = NA_real_,
    stringsAsFactors = FALSE
  )

  for (mm in 1:12) {
    ib <- which(b$month == mm)
    it <- which(t$month == mm)
    if (!length(it)) next

    bg <- .act_num(b[[gcm_col]][ib])
    bo <- .act_num(b[[obs_col]][ib])
    tg <- .act_num(t[[gcm_col]][it])

    ok_bg <- is.finite(bg)
    ok_bo <- is.finite(bo)
    if (!any(ok_bg) || !any(ok_bo)) {
      stop(sprintf("QDM: insufficient base data in calendar month %d", mm), call. = FALSE)
    }

    p <- rep(NA_real_, length(tg))
    ok_tg <- is.finite(tg)
    if (any(ok_tg)) {
      p[ok_tg] <- rank(tg[ok_tg], ties.method = "average") / (sum(ok_tg) + 1)
    }

    qg <- rep(NA_real_, length(tg))
    qo <- rep(NA_real_, length(tg))
    if (any(ok_tg)) {
      qg[ok_tg] <- as.numeric(stats::quantile(
        bg[ok_bg], probs = p[ok_tg], type = 7, names = FALSE, na.rm = TRUE
      ))
      qo[ok_tg] <- as.numeric(stats::quantile(
        bo[ok_bo], probs = p[ok_tg], type = 7, names = FALSE, na.rm = TRUE
      ))
    }

    corr <- rep(NA_real_, length(tg))
    mult <- ok_tg & is.finite(qg) & is.finite(qo) & qg > wet_threshold
    add  <- ok_tg & is.finite(qg) & is.finite(qo) & !mult
    corr[mult] <- qo[mult] * (tg[mult] / qg[mult])
    corr[add]  <- qo[add] + (tg[add] - qg[add])
    corr <- pmax(corr, 0)

    out$p_future_qdm[it] <- p
    out$q_gcm_hist[it] <- qg
    out$q_obs_hist[it] <- qo
    out$qdm_full[it] <- corr
  }

  out
}

# -----------------------------------------------------------------------------
# 3. Exact audited QMB-LIN Stage 2 wrapper
# -----------------------------------------------------------------------------
run_qmb_lin_5q_stage2 <- function(source_daily,
                                  target_daily,
                                  obs_col = "chirps",
                                  gcm_col = "gcm",
                                  wet_threshold = 0.1,
                                  taus = QMB_ACT_TAUS,
                                  random_state = 42L,
                                  target_quantile_fallback = TRUE,
                                  min_wet_target_quantile = 1L,
                                  fallback_min_pool_wet = 10L,
                                  fallback_max_year_radius = 5L,
                                  keep_plans = FALSE) {
  # Explicit use_log_rain=FALSE is deliberate: this is the locked final audited
  # baseline used in the dissertation experiments.
  fit <- fit_stage1_base(
    base_daily = source_daily,
    obs_col = obs_col,
    gcm_col = gcm_col,
    taus = taus,
    wet_threshold = wet_threshold,
    use_log_rain = FALSE
  )

  proj <- project_stage1_target(
    fit = fit,
    target_daily = target_daily,
    obs_col = obs_col,
    gcm_col = gcm_col,
    taus = taus,
    wet_threshold = wet_threshold,
    target_quantile_fallback = target_quantile_fallback,
    min_wet_target_quantile = min_wet_target_quantile,
    fallback_min_pool_wet = fallback_min_pool_wet,
    fallback_max_year_radius = fallback_max_year_radius
  )

  s2 <- run_stage2(
    base_daily = source_daily,
    target_daily = target_daily,
    fit = fit,
    projection = proj,
    obs_col = obs_col,
    taus = taus,
    wet_threshold = wet_threshold,
    random_state = as.integer(random_state) + 100L,
    keep_candidate_blocks = FALSE
  )

  d <- s2$point_daily
  names(d)[names(d) == "qmb_stage2"] <- "qmb_kokic5q_stage2"

  list(
    fit = fit,
    projection = proj,
    stage2 = s2,
    daily = d,
    plans = if (isTRUE(keep_plans)) s2$plans else NULL
  )
}

# Same wrapper when a full 1990-2014 Stage-1 fit is already available.
run_qmb_lin_5q_stage2_with_fit <- function(base_daily,
                                           target_daily,
                                           fit,
                                           obs_col = "chirps",
                                           gcm_col = "gcm",
                                           wet_threshold = 0.1,
                                           taus = QMB_ACT_TAUS,
                                           random_state = 42L,
                                           target_quantile_fallback = TRUE,
                                           min_wet_target_quantile = 1L,
                                           fallback_min_pool_wet = 10L,
                                           fallback_max_year_radius = 5L) {
  proj <- project_stage1_target(
    fit = fit,
    target_daily = target_daily,
    obs_col = obs_col,
    gcm_col = gcm_col,
    taus = taus,
    wet_threshold = wet_threshold,
    target_quantile_fallback = target_quantile_fallback,
    min_wet_target_quantile = min_wet_target_quantile,
    fallback_min_pool_wet = fallback_min_pool_wet,
    fallback_max_year_radius = fallback_max_year_radius
  )
  s2 <- run_stage2(
    base_daily = base_daily,
    target_daily = target_daily,
    fit = fit,
    projection = proj,
    obs_col = obs_col,
    taus = taus,
    wet_threshold = wet_threshold,
    random_state = as.integer(random_state) + 100L,
    keep_candidate_blocks = FALSE
  )
  d <- s2$point_daily
  names(d)[names(d) == "qmb_stage2"] <- "qmb_kokic5q_stage2"
  list(projection = proj, stage2 = s2, daily = d)
}

# -----------------------------------------------------------------------------
# 4. Hybrid QMB body + hard QDM tail at a month-specific cutoff
# -----------------------------------------------------------------------------
apply_qmb_qdm_tail <- function(qmb_daily,
                               qdm_daily,
                               cutoff_by_month,
                               output_col = "qmb_adaptive_stage2") {
  q <- qmb_daily[, c("date", "obs", "qmb_kokic5q_stage2"), drop = FALSE]
  d <- qdm_daily[, c("date", "year", "month", "day", "gcm", "p_future_qdm",
                     "q_gcm_hist", "q_obs_hist", "qdm_full"), drop = FALSE]
  z <- merge(d, q, by = "date", all.x = TRUE, sort = FALSE)
  z <- z[order(z$date), , drop = FALSE]

  cut <- if (is.data.frame(cutoff_by_month)) {
    if (!all(c("MONTH", "SELECTED_CUTOFF") %in% names(cutoff_by_month))) {
      stop("cutoff_by_month data.frame requires MONTH and SELECTED_CUTOFF", call. = FALSE)
    }
    setNames(.act_num(cutoff_by_month$SELECTED_CUTOFF), as.character(as.integer(cutoff_by_month$MONTH)))
  } else {
    x <- .act_num(cutoff_by_month)
    if (length(x) == 1L) x <- rep(x, 12L)
    if (length(x) != 12L) stop("cutoff_by_month must contain 12 monthly values", call. = FALSE)
    setNames(x, as.character(1:12))
  }

  z$selected_cutoff <- cut[as.character(z$month)]
  if (any(!is.finite(z$selected_cutoff))) {
    stop("Missing/non-finite monthly cutoff when applying adaptive tail", call. = FALSE)
  }
  z$use_qdm_tail <- is.finite(z$p_future_qdm) & z$p_future_qdm >= z$selected_cutoff
  z[[output_col]] <- z$qmb_kokic5q_stage2
  z[[output_col]][z$use_qdm_tail] <- z$qdm_full[z$use_qdm_tail]
  z
}

# -----------------------------------------------------------------------------
# 5. Calendar-month metrics used in CV and final evaluation
# -----------------------------------------------------------------------------
calendar_month_metrics <- function(daily,
                                   pred_col,
                                   obs_col = "obs",
                                   probs = c(0.95, 0.99)) {
  if (!all(c("date", "month", obs_col, pred_col) %in% names(daily))) {
    stop("calendar_month_metrics: required columns are missing", call. = FALSE)
  }
  rows <- vector("list", 12L)
  k <- 0L
  for (mm in 1:12) {
    g <- daily[daily$month == mm, , drop = FALSE]
    if (!nrow(g)) next
    y <- .act_num(g[[obs_col]])
    yh <- .act_num(g[[pred_col]])
    ok <- is.finite(y) & is.finite(yh)
    y <- y[ok]; yh <- yh[ok]
    if (!length(y)) next
    p95o <- .act_q(y, probs[1]); p95p <- .act_q(yh, probs[1])
    p99o <- .act_q(y, probs[2]); p99p <- .act_q(yh, probs[2])
    k <- k + 1L
    rows[[k]] <- data.frame(
      MONTH = mm,
      N_DAYS = length(y),
      MAE = .act_mae(y, yh),
      RMSE = .act_rmse(y, yh),
      BIAS = .act_bias(y, yh),
      P95_OBS = p95o,
      P95_PRED = p95p,
      P95_ERROR = p95p - p95o,
      P95_ABS_ERROR = abs(p95p - p95o),
      P99_OBS = p99o,
      P99_PRED = p99p,
      P99_ERROR = p99p - p99o,
      P99_ABS_ERROR = abs(p99p - p99o),
      stringsAsFactors = FALSE
    )
  }
  rows <- rows[seq_len(k)]
  .act_safe_rbind(rows)
}

# -----------------------------------------------------------------------------
# 6. Five blocked pseudo-future folds inside 1990-2014
# -----------------------------------------------------------------------------
make_blocked_cv_folds <- function(base_years = c(1990L, 2014L), block_size = 5L) {
  yy <- seq.int(as.integer(base_years[1]), as.integer(base_years[2]))
  if (length(yy) %% block_size != 0L) {
    stop("For the locked design, number of base years must be divisible by block_size", call. = FALSE)
  }
  split(yy, ceiling(seq_along(yy) / block_size))
}

# Run CV for one historical region-GCM series.
# historical_daily must contain ONLY the 1990-2014 base period (or at least be
# filterable to it). The held-out fold is pseudo-target; all other years remain
# source years. No 2015-2025 data enters this function.
cv_one_gcm <- function(historical_daily,
                       regime,
                       location_id,
                       gcm,
                       obs_col = "chirps",
                       gcm_col = "gcm",
                       base_years = c(1990L, 2014L),
                       folds = make_blocked_cv_folds(base_years, 5L),
                       cutoffs = QMB_ACT_CUTOFFS,
                       wet_threshold = 0.1,
                       random_state = 42L) {
  hist <- historical_daily[
    historical_daily$year >= base_years[1] & historical_daily$year <= base_years[2],
    , drop = FALSE
  ]
  if (!nrow(hist)) stop("cv_one_gcm: no historical rows", call. = FALSE)

  detail <- list(); audit <- list(); k <- 0L; a <- 0L

  for (fi in seq_along(folds)) {
    valid_years <- as.integer(folds[[fi]])
    source <- hist[!hist$year %in% valid_years, , drop = FALSE]
    pseudo <- hist[ hist$year %in% valid_years, , drop = FALSE]

    seed_fold <- as.integer(random_state) + fi * 10000L + sum(utf8ToInt(as.character(gcm)))
    t0 <- proc.time()[["elapsed"]]

    res <- run_qmb_lin_5q_stage2(
      source_daily = source,
      target_daily = pseudo,
      obs_col = obs_col,
      gcm_col = gcm_col,
      wet_threshold = wet_threshold,
      taus = QMB_ACT_TAUS,
      random_state = seed_fold
    )
    qdm <- run_empirical_qdm(source, pseudo, obs_col, gcm_col, wet_threshold)

    qmb <- merge(
      qdm[, c("date", "year", "month", "day", "obs", "gcm", "p_future_qdm",
              "q_gcm_hist", "q_obs_hist", "qdm_full")],
      res$daily[, c("date", "qmb_kokic5q_stage2")],
      by = "date", all.x = TRUE, sort = FALSE
    )
    qmb <- qmb[order(qmb$date), , drop = FALSE]

    m_qdm <- calendar_month_metrics(qmb, "qdm_full", "obs")
    m_qmb <- calendar_month_metrics(qmb, "qmb_kokic5q_stage2", "obs")

    for (cc in cutoffs) {
      h <- qmb
      h$hybrid <- h$qmb_kokic5q_stage2
      use <- is.finite(h$p_future_qdm) & h$p_future_qdm >= cc
      h$hybrid[use] <- h$qdm_full[use]
      mh <- calendar_month_metrics(h, "hybrid", "obs")

      for (mm in 1:12) {
        hh <- mh[mh$MONTH == mm, , drop = FALSE]
        qq <- m_qdm[m_qdm$MONTH == mm, , drop = FALSE]
        bb <- m_qmb[m_qmb$MONTH == mm, , drop = FALSE]
        if (!nrow(hh) || !nrow(qq) || !nrow(bb)) next
        k <- k + 1L
        detail[[k]] <- data.frame(
          REGIME = regime,
          LOCATION_ID = location_id,
          GCM = .act_norm_gcm(gcm),
          FOLD = fi,
          VALID_YEARS = paste(range(valid_years), collapse = "-"),
          MONTH = mm,
          CUTOFF = cc,
          N_DAYS = hh$N_DAYS[1],
          MAE = hh$MAE[1],
          RMSE = hh$RMSE[1],
          BIAS = hh$BIAS[1],
          P95_ABS_ERROR = hh$P95_ABS_ERROR[1],
          P99_ABS_ERROR = hh$P99_ABS_ERROR[1],
          QDM_P95_ABS_ERROR = qq$P95_ABS_ERROR[1],
          QDM_P99_ABS_ERROR = qq$P99_ABS_ERROR[1],
          QMB_P95_ABS_ERROR = bb$P95_ABS_ERROR[1],
          QMB_P99_ABS_ERROR = bb$P99_ABS_ERROR[1],
          stringsAsFactors = FALSE
        )
      }
    }

    a <- a + 1L
    audit[[a]] <- data.frame(
      REGIME = regime,
      LOCATION_ID = location_id,
      GCM = .act_norm_gcm(gcm),
      FOLD = fi,
      VALID_YEARS = paste(range(valid_years), collapse = "-"),
      N_SOURCE_DAYS = nrow(source),
      N_PSEUDOTARGET_DAYS = nrow(pseudo),
      ELAPSED_SECONDS = proc.time()[["elapsed"]] - t0,
      STATUS = "SUCCESS",
      stringsAsFactors = FALSE
    )
  }

  list(detail = .act_safe_rbind(detail), audit = .act_safe_rbind(audit))
}

# -----------------------------------------------------------------------------
# 7. Aggregate CV scores and select cutoff per LOCATION x MONTH
# -----------------------------------------------------------------------------
aggregate_cv_scores <- function(cv_detail) {
  if (!nrow(cv_detail)) stop("aggregate_cv_scores: empty cv_detail", call. = FALSE)

  # GCM x month x cutoff summary across folds.
  gkeys <- interaction(
    cv_detail$REGIME, cv_detail$LOCATION_ID, cv_detail$GCM,
    cv_detail$MONTH, cv_detail$CUTOFF,
    drop = TRUE, lex.order = TRUE
  )
  gg <- split(seq_len(nrow(cv_detail)), gkeys)
  gcm_rows <- lapply(gg, function(idx) {
    x <- cv_detail[idx, , drop = FALSE]
    data.frame(
      REGIME = x$REGIME[1], LOCATION_ID = x$LOCATION_ID[1], GCM = x$GCM[1],
      MONTH = x$MONTH[1], CUTOFF = x$CUTOFF[1],
      CV_MAE = mean(x$MAE, na.rm = TRUE),
      CV_RMSE = mean(x$RMSE, na.rm = TRUE),
      CV_P95_MAE = mean(x$P95_ABS_ERROR, na.rm = TRUE),
      CV_P99_MAE = mean(x$P99_ABS_ERROR, na.rm = TRUE),
      QDM_P95_MAE = mean(x$QDM_P95_ABS_ERROR, na.rm = TRUE),
      QDM_P99_MAE = mean(x$QDM_P99_ABS_ERROR, na.rm = TRUE),
      QMB_P95_MAE = mean(x$QMB_P95_ABS_ERROR, na.rm = TRUE),
      QMB_P99_MAE = mean(x$QMB_P99_ABS_ERROR, na.rm = TRUE),
      N_FOLDS = nrow(x), stringsAsFactors = FALSE
    )
  })

  # Location x month x cutoff summary across all GCM-fold components.
  lkeys <- interaction(
    cv_detail$REGIME, cv_detail$LOCATION_ID, cv_detail$MONTH, cv_detail$CUTOFF,
    drop = TRUE, lex.order = TRUE
  )
  ll <- split(seq_len(nrow(cv_detail)), lkeys)
  loc_rows <- lapply(ll, function(idx) {
    x <- cv_detail[idx, , drop = FALSE]
    data.frame(
      REGIME = x$REGIME[1], LOCATION_ID = x$LOCATION_ID[1],
      MONTH = x$MONTH[1], CUTOFF = x$CUTOFF[1],
      CV_MAE = mean(x$MAE, na.rm = TRUE),
      CV_RMSE = mean(x$RMSE, na.rm = TRUE),
      CV_P95_MAE = mean(x$P95_ABS_ERROR, na.rm = TRUE),
      CV_P99_MAE = mean(x$P99_ABS_ERROR, na.rm = TRUE),
      QDM_P95_MAE = mean(x$QDM_P95_ABS_ERROR, na.rm = TRUE),
      QDM_P99_MAE = mean(x$QDM_P99_ABS_ERROR, na.rm = TRUE),
      QMB_P95_MAE = mean(x$QMB_P95_ABS_ERROR, na.rm = TRUE),
      QMB_P99_MAE = mean(x$QMB_P99_ABS_ERROR, na.rm = TRUE),
      N_FOLD_GCM = nrow(x), stringsAsFactors = FALSE
    )
  })

  list(
    gcm_month_cutoff = .act_safe_rbind(gcm_rows),
    location_month_cutoff = .act_safe_rbind(loc_rows)
  )
}

select_cutoff_robust_allfolds <- function(cv_detail,
                                          tolerance_mm = 1.0,
                                          cutoffs = QMB_ACT_CUTOFFS) {
  combos <- unique(cv_detail[, c("REGIME", "LOCATION_ID", "MONTH"), drop = FALSE])
  combos <- combos[order(combos$REGIME, combos$LOCATION_ID, combos$MONTH), , drop = FALSE]
  out <- vector("list", nrow(combos))

  for (i in seq_len(nrow(combos))) {
    reg <- combos$REGIME[i]; loc <- combos$LOCATION_ID[i]; mm <- combos$MONTH[i]
    x <- cv_detail[
      cv_detail$REGIME == reg & cv_detail$LOCATION_ID == loc & cv_detail$MONTH == mm,
      , drop = FALSE
    ]

    cand <- vector("list", length(cutoffs))
    for (j in seq_along(cutoffs)) {
      cc <- cutoffs[j]
      z <- x[abs(x$CUTOFF - cc) < 1e-12, , drop = FALSE]
      ex95 <- z$P95_ABS_ERROR - z$QDM_P95_ABS_ERROR
      ex99 <- z$P99_ABS_ERROR - z$QDM_P99_ABS_ERROR
      worst <- max(c(ex95, ex99), na.rm = TRUE)
      feasible <- all(ex95 <= tolerance_mm & ex99 <= tolerance_mm, na.rm = TRUE)
      cand[[j]] <- data.frame(
        CUTOFF = cc,
        FEASIBLE = feasible,
        WORST_EXCESS = worst,
        CV_MAE = mean(z$MAE, na.rm = TRUE),
        CV_RMSE = mean(z$RMSE, na.rm = TRUE),
        CV_P95 = mean(z$P95_ABS_ERROR, na.rm = TRUE),
        CV_P99 = mean(z$P99_ABS_ERROR, na.rm = TRUE),
        QDM_P95 = mean(z$QDM_P95_ABS_ERROR, na.rm = TRUE),
        QDM_P99 = mean(z$QDM_P99_ABS_ERROR, na.rm = TRUE),
        N_COMPONENTS = nrow(z), stringsAsFactors = FALSE
      )
    }
    cand <- .act_safe_rbind(cand)

    if (any(cand$FEASIBLE)) {
      feasible <- cand[cand$FEASIBLE, , drop = FALSE]
      chosen <- feasible[which.max(feasible$CUTOFF), , drop = FALSE]
      status <- "FEASIBLE_ALL_FOLDS"
    } else {
      # Deterministic fallback: minimize worst tail excess; then RMSE; then use
      # the highest cutoff among remaining exact ties.
      best_w <- min(cand$WORST_EXCESS, na.rm = TRUE)
      z <- cand[abs(cand$WORST_EXCESS - best_w) < 1e-12, , drop = FALSE]
      best_r <- min(z$CV_RMSE, na.rm = TRUE)
      z <- z[abs(z$CV_RMSE - best_r) < 1e-12, , drop = FALSE]
      chosen <- z[which.max(z$CUTOFF), , drop = FALSE]
      status <- "FALLBACK_MIN_WORST_EXCESS"
    }

    out[[i]] <- data.frame(
      REGIME = reg,
      LOCATION_ID = loc,
      MONTH = mm,
      SELECTED_CUTOFF = chosen$CUTOFF[1],
      STATUS = status,
      WORST_EXCESS = chosen$WORST_EXCESS[1],
      CV_MAE = chosen$CV_MAE[1],
      CV_RMSE = chosen$CV_RMSE[1],
      CV_P95 = chosen$CV_P95[1],
      CV_P99 = chosen$CV_P99[1],
      QDM_P95 = chosen$QDM_P95[1],
      QDM_P99 = chosen$QDM_P99[1],
      TOLERANCE_MM = tolerance_mm,
      stringsAsFactors = FALSE
    )
  }

  ans <- .act_safe_rbind(out)
  ans$MONTH_NAME <- month.abb[ans$MONTH]
  ans
}


# -----------------------------------------------------------------------------
# 8. Stage-3 adaptive bootstrap helpers
# -----------------------------------------------------------------------------
# Stage 3 is generated from the SAME QMB-LIN source-year bootstrap ensemble
# for both methods. Adaptive-CV differs only by the locked Stage-2 correction.
#
# This gives a paired comparison:
#   QMB-LIN Stage 3       : Y_LIN^(b)
#   Adaptive-CV Stage 3  : Y_LIN^(b) + (Y_ACT,2 - Y_LIN,2)
#
# Because the same bootstrap draws are used, differences between the two
# methods are caused by the adaptive interpolation/tail treatment rather than
# different random resampling sequences.

.act_summarize_stage3_ensemble <- function(ensemble,
                                           target_daily,
                                           obs_col = "chirps",
                                           interval = 0.95) {
  if (!is.matrix(ensemble)) ensemble <- as.matrix(ensemble)
  if (nrow(ensemble) < 2L) stop("Stage3 ensemble requires >=2 bootstrap members", call. = FALSE)
  if (ncol(ensemble) != nrow(target_daily)) {
    stop("Stage3 ensemble column count does not match target_daily rows", call. = FALSE)
  }
  if (!(interval > 0 && interval < 1)) stop("interval must be in (0,1)", call. = FALSE)

  alpha <- 1 - interval
  qcol <- function(p) {
    apply(
      ensemble, 2, stats::quantile,
      probs = p, na.rm = TRUE, names = FALSE, type = 7
    )
  }

  lower <- qcol(alpha / 2)
  med   <- qcol(0.5)
  upper <- qcol(1 - alpha / 2)

  td <- target_daily[order(target_daily$date), , drop = FALSE]
  rownames(td) <- NULL

  daily_interval <- data.frame(
    date = td$date,
    year = td$year,
    month = td$month,
    day = td$day,
    obs = .act_num(td[[obs_col]]),
    lower = lower,
    median = med,
    upper = upper,
    stringsAsFactors = FALSE
  )

  valid <- complete.cases(daily_interval[, c("obs", "lower", "median", "upper")])
  inside <- daily_interval$obs[valid] >= daily_interval$lower[valid] &
    daily_interval$obs[valid] <= daily_interval$upper[valid]
  widths <- daily_interval$upper[valid] - daily_interval$lower[valid]

  e <- daily_interval$median[valid] - daily_interval$obs[valid]
  accuracy <- data.frame(
    n_days = sum(valid),
    sum_abs_error = sum(abs(e)),
    sum_sq_error = sum(e^2),
    sum_error = sum(e),
    MAE = if (length(e)) mean(abs(e)) else NA_real_,
    RMSE = if (length(e)) sqrt(mean(e^2)) else NA_real_,
    Bias = if (length(e)) mean(e) else NA_real_,
    stringsAsFactors = FALSE
  )

  interval_metrics <- data.frame(
    nominal_interval = interval,
    coverage = if (length(inside)) mean(inside) else NA_real_,
    n_inside = if (length(inside)) sum(inside) else 0L,
    n_interval_days = length(inside),
    coverage_error_signed =
      if (length(inside)) mean(inside) - interval else NA_real_,
    coverage_error_abs =
      if (length(inside)) abs(mean(inside) - interval) else NA_real_,
    mean_interval_width = if (length(widths)) mean(widths) else NA_real_,
    median_interval_width = if (length(widths)) stats::median(widths) else NA_real_,
    sum_interval_width = if (length(widths)) sum(widths) else 0,
    n_boot = nrow(ensemble),
    stringsAsFactors = FALSE
  )

  list(
    daily_interval = daily_interval,
    accuracy = accuracy,
    interval_metrics = interval_metrics
  )
}

.act_recenter_adaptive_ensemble <- function(qmb_ensemble,
                                            qmb_stage2_daily,
                                            adaptive_stage2_daily,
                                            target_daily,
                                            obs_col = "chirps",
                                            interval = 0.95,
                                            qmb_col = "qmb_kokic5q_stage2",
                                            adaptive_col = "qmb_adaptive_stage2",
                                            clip_zero = TRUE) {
  q2 <- qmb_stage2_daily[, c("date", qmb_col), drop = FALSE]
  a2_cols <- c("date", adaptive_col, "selected_cutoff", "use_qdm_tail",
               "p_future_qdm", "qdm_full")
  a2_cols <- a2_cols[a2_cols %in% names(adaptive_stage2_daily)]
  a2 <- adaptive_stage2_daily[, a2_cols, drop = FALSE]

  center <- merge(q2, a2, by = "date", all = FALSE, sort = FALSE)
  center <- center[order(center$date), , drop = FALSE]
  td <- target_daily[order(target_daily$date), , drop = FALSE]

  if (nrow(center) != nrow(td) || !identical(as.Date(center$date), as.Date(td$date))) {
    stop("Stage3 recentering date alignment failed", call. = FALSE)
  }

  shift <- .act_num(center[[adaptive_col]]) - .act_num(center[[qmb_col]])
  if (any(!is.finite(shift))) {
    stop("Non-finite Stage2 adaptive shift encountered before Stage3", call. = FALSE)
  }

  ens <- sweep(as.matrix(qmb_ensemble), 2, shift, FUN = "+")
  if (isTRUE(clip_zero)) ens <- pmax(ens, 0)

  out <- .act_summarize_stage3_ensemble(
    ensemble = ens,
    target_daily = td,
    obs_col = obs_col,
    interval = interval
  )

  # Add the locked adaptive metadata to the daily Stage-3 output.
  meta_cols <- setdiff(names(center), c("date", qmb_col, adaptive_col))
  if (length(meta_cols)) {
    out$daily_interval <- merge(
      out$daily_interval,
      center[, c("date", meta_cols), drop = FALSE],
      by = "date", all.x = TRUE, sort = FALSE
    )
    out$daily_interval <- out$daily_interval[order(out$daily_interval$date), , drop = FALSE]
  }

  out$shift <- shift
  out$ensemble <- ens
  out
}

.act_probabilistic_calibration <- function(ensemble, target_daily, obs_col = "chirps",
                                           levels = c(0.50, 0.80, 0.90, 0.95)) {
  if (!is.matrix(ensemble)) ensemble <- as.matrix(ensemble)
  td <- target_daily[order(target_daily$date), , drop = FALSE]
  y <- .act_num(td[[obs_col]])
  ok <- is.finite(y)
  ens <- ensemble[, ok, drop = FALSE]
  y <- y[ok]
  if (!length(y)) return(data.frame())
  rows <- lapply(levels, function(level) {
    a <- 1 - level
    lo <- apply(ens, 2, stats::quantile, probs = a/2, type = 7, na.rm = TRUE, names = FALSE)
    hi <- apply(ens, 2, stats::quantile, probs = 1-a/2, type = 7, na.rm = TRUE, names = FALSE)
    inside <- y >= lo & y <= hi
    width <- hi-lo
    # Gneiting-Raftery central interval score; lower is better.
    iscore <- width + (2/a)*(lo-y)*(y < lo) + (2/a)*(y-hi)*(y > hi)
    data.frame(SCOPE="ALL", NOMINAL=level, N=length(y), COVERAGE=mean(inside),
               COVERAGE_ERROR=mean(inside)-level, ABS_COVERAGE_ERROR=abs(mean(inside)-level),
               MEAN_WIDTH=mean(width), INTERVAL_SCORE=mean(iscore), stringsAsFactors=FALSE)
  })
  out <- .act_safe_rbind(rows)
  # Upper-tail calibration: observations at/above their pooled P90/P95/P99.
  for (q in c(.90,.95,.99)) {
    thr <- as.numeric(stats::quantile(y, q, type=7, na.rm=TRUE))
    ix <- which(y >= thr)
    if (length(ix)) {
      e <- ens[,ix,drop=FALSE]; yy <- y[ix]
      a <- .05
      lo <- apply(e,2,stats::quantile,probs=a/2,type=7,na.rm=TRUE,names=FALSE)
      hi <- apply(e,2,stats::quantile,probs=1-a/2,type=7,na.rm=TRUE,names=FALSE)
      inside <- yy>=lo & yy<=hi; width <- hi-lo
      iscore <- width + (2/a)*(lo-yy)*(yy<lo) + (2/a)*(yy-hi)*(yy>hi)
      out <- rbind(out, data.frame(SCOPE=paste0("OBS_GE_P",round(100*q)), NOMINAL=.95,
        N=length(ix), COVERAGE=mean(inside), COVERAGE_ERROR=mean(inside)-.95,
        ABS_COVERAGE_ERROR=abs(mean(inside)-.95), MEAN_WIDTH=mean(width),
        INTERVAL_SCORE=mean(iscore), stringsAsFactors=FALSE))
    }
  }
  out
}

.act_stage3_monthly_metrics <- function(daily_interval) {
  x <- daily_interval
  x$stage3_median <- x$median
  calendar_month_metrics(
    daily = x,
    pred_col = "stage3_median",
    obs_col = "obs"
  )
}

.act_stage3_overall_from_rows <- function(acc_rows, monthly_rows, interval_rows) {
  if (!nrow(acc_rows)) return(data.frame())

  keys <- unique(acc_rows[, c("METHOD", "N_BOOT"), drop = FALSE])
  out <- vector("list", nrow(keys))

  for (i in seq_len(nrow(keys))) {
    meth <- keys$METHOD[i]
    B <- keys$N_BOOT[i]

    aa <- acc_rows[
      acc_rows$METHOD == meth & acc_rows$N_BOOT == B,
      , drop = FALSE
    ]
    mm <- monthly_rows[
      monthly_rows$METHOD == meth & monthly_rows$N_BOOT == B,
      , drop = FALSE
    ]
    ii <- interval_rows[
      interval_rows$METHOD == meth & interval_rows$N_BOOT == B,
      , drop = FALSE
    ]

    n <- sum(aa$n_days, na.rm = TRUE)
    ni <- sum(ii$n_interval_days, na.rm = TRUE)

    out[[i]] <- data.frame(
      METHOD = meth,
      N_BOOT = B,
      DAILY_MAE =
        if (n > 0) sum(aa$sum_abs_error, na.rm = TRUE) / n else NA_real_,
      DAILY_RMSE =
        if (n > 0) sqrt(sum(aa$sum_sq_error, na.rm = TRUE) / n) else NA_real_,
      DAILY_BIAS =
        if (n > 0) sum(aa$sum_error, na.rm = TRUE) / n else NA_real_,
      MEAN_MONTHLY_P95_ABS_ERROR =
        mean(mm$P95_ABS_ERROR, na.rm = TRUE),
      MEAN_MONTHLY_P99_ABS_ERROR =
        mean(mm$P99_ABS_ERROR, na.rm = TRUE),
      COVERAGE =
        if (ni > 0) sum(ii$n_inside, na.rm = TRUE) / ni else NA_real_,
      MEAN_INTERVAL_WIDTH =
        if (ni > 0) sum(ii$sum_interval_width, na.rm = TRUE) / ni else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  .act_safe_rbind(out)
}


# -----------------------------------------------------------------------------
# 8. Manifest validation
# -----------------------------------------------------------------------------
# Required columns:
#   REGIME, LOCATION_ID, GCM, SSP, FILE
# Each row points to one paired daily CSV containing date, chirps, gcm.
validate_manifest <- function(manifest) {
  req <- c("REGIME", "LOCATION_ID", "GCM", "SSP", "FILE")
  miss <- setdiff(req, names(manifest))
  if (length(miss)) stop("Manifest missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  m <- manifest[, req, drop = FALSE]
  m$REGIME <- as.character(m$REGIME)
  m$LOCATION_ID <- as.character(m$LOCATION_ID)
  m$GCM <- .act_norm_gcm(m$GCM)
  m$SSP <- .act_norm_ssp(m$SSP)
  m$FILE <- path.expand(as.character(m$FILE))
  if (any(!file.exists(m$FILE))) {
    stop("Manifest contains missing files:\n", paste(m$FILE[!file.exists(m$FILE)], collapse = "\n"), call. = FALSE)
  }
  key <- paste(m$REGIME, m$LOCATION_ID, m$GCM, m$SSP, sep = "|")
  if (anyDuplicated(key)) stop("Manifest has duplicate REGIME/LOCATION/GCM/SSP combinations", call. = FALSE)
  m[order(m$REGIME, m$LOCATION_ID, m$GCM, m$SSP), , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 9. High-level end-to-end runner
# -----------------------------------------------------------------------------
run_qmb_adaptive_qdm_tail <- function(manifest,
                                      output_root,
                                      core_file = file.path(.qmb_act_dir, "QMB_KOKIC_CORE.R"),
                                      date_col = "date",
                                      obs_col = "chirps",
                                      gcm_col = "gcm",
                                      base_years = c(1990L, 2014L),
                                      target_years = c(2015L, 2025L),
                                      no_leap = TRUE,
                                      wet_threshold = 0.1,
                                      cutoffs = QMB_ACT_CUTOFFS,
                                      cv_tolerance_mm = 1.0,
                                      also_run_strict_05mm = TRUE,
                                      random_state = 42L,
                                      run_stage3_bootstrap = TRUE,
                                      b_values = c(30L, 100L, 250L, 500L, 1000L, 2000L),
                                      stage3_interval = 0.95,
                                      save_stage3_daily_all_b = FALSE,
                                      keep_stage3_ensemble_rds = FALSE,
                                      save_rds = TRUE) {
  load_qmb_audited_core(core_file)
  m <- validate_manifest(manifest)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  writeLines(c(
    paste0("METHOD_ID=", QMB_ACT_METHOD_ID),
    paste0("VERSION=", QMB_ACT_VERSION),
    paste0("BASE_YEARS=", paste(base_years, collapse = "-")),
    paste0("TARGET_YEARS=", paste(target_years, collapse = "-")),
    paste0("TAUS=", paste(QMB_ACT_TAUS, collapse = ",")),
    paste0("CUTOFFS=", paste(sprintf("%.3f", cutoffs), collapse = ",")),
    paste0("CV_TOLERANCE_MM=", cv_tolerance_mm),
    paste0("WET_THRESHOLD=", wet_threshold),
    "SELECTOR=ROBUST_ALL_FOLDS_HIGHEST_FEASIBLE_CUTOFF",
    "TARGET_2015_2025_USED_FOR_CUTOFF_SELECTION=NO",
    paste0("RUN_STAGE3_BOOTSTRAP=", isTRUE(run_stage3_bootstrap)),
    paste0("STAGE3_B_VALUES=", paste(as.integer(b_values), collapse = ",")),
    paste0("STAGE3_INTERVAL=", stage3_interval),
    "STAGE3_ADAPTIVE=PAIRED_QMB_BOOTSTRAP_PLUS_LOCKED_STAGE2_RECENTER_SHIFT"
  ), file.path(output_root, "METHOD_LOCK.txt"))
  utils::write.csv(m, file.path(output_root, "INPUT_MANIFEST.csv"), row.names = FALSE)

  # Load all paired files once.
  daily <- vector("list", nrow(m))
  for (i in seq_len(nrow(m))) {
    daily[[i]] <- load_pair_csv(
      m$FILE[i], date_col = date_col, obs_col = obs_col,
      gcm_col = gcm_col, no_leap = no_leap
    )
  }

  # ---------------------------------------------------------------------------
  # A. Historical CV. Use ONE historical series per regime-location-GCM,
  #    because SSP historical data must be identical.
  # ---------------------------------------------------------------------------
  gcm_keys <- unique(m[, c("REGIME", "LOCATION_ID", "GCM"), drop = FALSE])
  cv_details <- list(); cv_audits <- list(); cd <- 0L; ca <- 0L
  shared_base_audit <- list(); ba <- 0L

  for (ii in seq_len(nrow(gcm_keys))) {
    reg <- gcm_keys$REGIME[ii]; loc <- gcm_keys$LOCATION_ID[ii]; gg <- gcm_keys$GCM[ii]
    ids <- which(m$REGIME == reg & m$LOCATION_ID == loc & m$GCM == gg)

    # Verify that 1990-2014 is identical across SSPs before sharing a fit/CV.
    bases <- list()
    for (j in ids) {
      sp <- split_period(daily[[j]], base_years, target_years)
      bases[[m$SSP[j]]] <- sp$base
    }
    aud <- verify_shared_base_across_ssps(
      bases, date_col = "date", obs_col = obs_col, gcm_col = gcm_col
    )
    aud$REGIME <- reg; aud$LOCATION_ID <- loc; aud$GCM <- gg
    ba <- ba + 1L; shared_base_audit[[ba]] <- aud

    representative <- bases[[1]]
    cv <- cv_one_gcm(
      historical_daily = representative,
      regime = reg,
      location_id = loc,
      gcm = gg,
      obs_col = obs_col,
      gcm_col = gcm_col,
      base_years = base_years,
      folds = make_blocked_cv_folds(base_years, 5L),
      cutoffs = cutoffs,
      wet_threshold = wet_threshold,
      random_state = as.integer(random_state) + ii * 100000L
    )
    cd <- cd + 1L; cv_details[[cd]] <- cv$detail
    ca <- ca + 1L; cv_audits[[ca]] <- cv$audit
  }

  cv_detail <- .act_safe_rbind(cv_details)
  cv_audit <- .act_safe_rbind(cv_audits)
  base_audit <- .act_safe_rbind(shared_base_audit)
  cv_scores <- aggregate_cv_scores(cv_detail)
  selected <- select_cutoff_robust_allfolds(cv_detail, cv_tolerance_mm, cutoffs)
  selected_strict <- if (isTRUE(also_run_strict_05mm)) {
    select_cutoff_robust_allfolds(cv_detail, 0.5, cutoffs)
  } else NULL

  utils::write.csv(cv_detail, file.path(output_root, "CV_DETAIL_ALL_GCM_FOLD_MONTH_CUTOFF.csv"), row.names = FALSE)
  utils::write.csv(cv_audit, file.path(output_root, "CV_RUN_AUDIT.csv"), row.names = FALSE)
  utils::write.csv(base_audit, file.path(output_root, "SHARED_BASE_SSP_AUDIT.csv"), row.names = FALSE)
  utils::write.csv(cv_scores$gcm_month_cutoff, file.path(output_root, "CV_SCORES_GCM_MONTH_CUTOFF.csv"), row.names = FALSE)
  utils::write.csv(cv_scores$location_month_cutoff, file.path(output_root, "CV_SCORES_LOCATION_MONTH_CUTOFF.csv"), row.names = FALSE)
  utils::write.csv(selected, file.path(output_root, "SELECTED_CUTOFF_LOCATION_MONTH_ROBUST.csv"), row.names = FALSE)
  if (!is.null(selected_strict)) {
    utils::write.csv(selected_strict, file.path(output_root, "SELECTED_CUTOFF_LOCATION_MONTH_STRICT_05MM.csv"), row.names = FALSE)
  }

  # ---------------------------------------------------------------------------
  # B. Final 1990-2014 fit and locked application to 2015-2025.
  # ---------------------------------------------------------------------------
  final_daily_all <- list(); final_month_all <- list(); fd <- 0L; fm <- 0L
  fits <- list()

  # Stage-3 masters. QMB-LIN and Adaptive-CV use exactly the same bootstrap
  # source-year draws. We run once at max(B) and summarize nested prefixes.
  stage3_accuracy_all <- list(); s3a <- 0L
  stage3_interval_all <- list(); s3i <- 0L
  stage3_month_all <- list(); s3m <- 0L
  stage3_calibration_all <- list(); s3c <- 0L
  b_values <- sort(unique(as.integer(b_values)))
  if (isTRUE(run_stage3_bootstrap)) {
    if (!length(b_values) || any(!is.finite(b_values)) || any(b_values < 2L)) {
      stop("b_values must contain integers >= 2", call. = FALSE)
    }
    if (!(stage3_interval > 0 && stage3_interval < 1)) {
      stop("stage3_interval must be in (0,1)", call. = FALSE)
    }
  }

  for (ii in seq_len(nrow(gcm_keys))) {
    reg <- gcm_keys$REGIME[ii]; loc <- gcm_keys$LOCATION_ID[ii]; gg <- gcm_keys$GCM[ii]
    ids <- which(m$REGIME == reg & m$LOCATION_ID == loc & m$GCM == gg)

    sp0 <- split_period(daily[[ids[1]]], base_years, target_years)
    base <- sp0$base
    fit_key <- paste(reg, loc, gg, sep = "|")
    fits[[fit_key]] <- fit_stage1_base(
      base_daily = base,
      obs_col = obs_col,
      gcm_col = gcm_col,
      taus = QMB_ACT_TAUS,
      wet_threshold = wet_threshold,
      use_log_rain = FALSE
    )

    cut_loc <- selected[selected$REGIME == reg & selected$LOCATION_ID == loc, c("MONTH", "SELECTED_CUTOFF"), drop = FALSE]
    if (nrow(cut_loc) != 12L) {
      stop(sprintf("Expected 12 selected cutoffs for %s / %s; found %d", reg, loc, nrow(cut_loc)), call. = FALSE)
    }

    for (j in ids) {
      ssp <- m$SSP[j]
      sp <- split_period(daily[[j]], base_years, target_years)
      target <- sp$target

      qmb <- run_qmb_lin_5q_stage2_with_fit(
        base_daily = base,
        target_daily = target,
        fit = fits[[fit_key]],
        obs_col = obs_col,
        gcm_col = gcm_col,
        wet_threshold = wet_threshold,
        taus = QMB_ACT_TAUS,
        random_state = as.integer(random_state) + ii * 100000L + j * 1000L
      )
      qdm <- run_empirical_qdm(base, target, obs_col, gcm_col, wet_threshold)
      adaptive <- apply_qmb_qdm_tail(qmb$daily, qdm, cut_loc, "qmb_adaptive_stage2")

      # Auditable comparator columns.
      adaptive$qmb_hard90_stage2 <- adaptive$qmb_kokic5q_stage2
      use90 <- is.finite(adaptive$p_future_qdm) & adaptive$p_future_qdm >= 0.90
      adaptive$qmb_hard90_stage2[use90] <- adaptive$qdm_full[use90]

      adaptive$qmb_fixed94_stage2 <- adaptive$qmb_kokic5q_stage2
      use94 <- is.finite(adaptive$p_future_qdm) & adaptive$p_future_qdm >= 0.94
      adaptive$qmb_fixed94_stage2[use94] <- adaptive$qdm_full[use94]

      adaptive$REGIME <- reg
      adaptive$LOCATION_ID <- loc
      adaptive$GCM <- gg
      adaptive$SSP <- ssp
      adaptive$METHOD_ID <- QMB_ACT_METHOD_ID
      adaptive$METHOD_VERSION <- QMB_ACT_VERSION

      outdir <- file.path(output_root, "FINAL_STAGE2", reg, loc, gg, ssp)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(adaptive, file.path(outdir, "stage2_daily_qmb_adaptive_qdm_tail.csv"), row.names = FALSE)
      utils::write.csv(qmb$projection$projected, file.path(outdir, "stage1_projected_5q.csv"), row.names = FALSE)

      methods <- c(
        QMB = "qmb_kokic5q_stage2",
        HARD90 = "qmb_hard90_stage2",
        FIXED94 = "qmb_fixed94_stage2",
        ADAPTIVE_CV = "qmb_adaptive_stage2",
        QDM = "qdm_full"
      )
      for (mn in names(methods)) {
        mm <- calendar_month_metrics(adaptive, methods[[mn]], "obs")
        mm$REGIME <- reg; mm$LOCATION_ID <- loc; mm$GCM <- gg; mm$SSP <- ssp
        mm$METHOD <- mn
        mm$CV_CUTOFF <- if (mn == "ADAPTIVE_CV") {
          cut_loc$SELECTED_CUTOFF[match(mm$MONTH, cut_loc$MONTH)]
        } else NA_real_
        fm <- fm + 1L; final_month_all[[fm]] <- mm
      }

      fd <- fd + 1L
      final_daily_all[[fd]] <- adaptive[, c(
        "date", "year", "month", "day", "obs", "gcm", "p_future_qdm",
        "q_gcm_hist", "q_obs_hist", "qdm_full", "qmb_kokic5q_stage2",
        "qmb_hard90_stage2", "qmb_fixed94_stage2", "selected_cutoff",
        "use_qdm_tail", "qmb_adaptive_stage2", "REGIME", "LOCATION_ID",
        "GCM", "SSP", "METHOD_ID", "METHOD_VERSION"
      ), drop = FALSE]

      # -----------------------------------------------------------------------
      # C. Stage 3 — paired source-year bootstrap comparison.
      # -----------------------------------------------------------------------
      if (isTRUE(run_stage3_bootstrap)) {
        max_b <- max(b_values)

        # Audited QMB-LIN Stage 3. The precomputed Stage-2 plans guarantee that
        # Stage 2 is not recomputed inside bootstrap.
        qmb_s3_max <- run_stage3(
          base_daily = base,
          target_daily = target,
          fit = fits[[fit_key]],
          projection = qmb$projection,
          obs_col = obs_col,
          taus = QMB_ACT_TAUS,
          wet_threshold = wet_threshold,
          n_boot = max_b,
          interval = stage3_interval,
          random_state = as.integer(random_state) + ii * 100000L + j * 1000L + 200L,
          keep_ensemble = TRUE,
          precomputed_plans = qmb$stage2$plans
        )

        ens_max <- qmb_s3_max$ensemble
        if (is.null(ens_max)) {
          stop("QMB Stage3 did not return ensemble although keep_ensemble=TRUE", call. = FALSE)
        }

        # Optional large ensemble audit. Saved only for max(B).
        if (isTRUE(keep_stage3_ensemble_rds)) {
          s3_ens_dir <- file.path(
            output_root, "FINAL_STAGE3", reg, loc, gg, ssp,
            sprintf("B_%04d", max_b)
          )
          dir.create(s3_ens_dir, recursive = TRUE, showWarnings = FALSE)
          saveRDS(
            ens_max,
            file.path(s3_ens_dir, "qmb_lin_stage3_ensemble_maxB.rds")
          )
        }

        for (B in b_values) {
          ens_lin <- ens_max[seq_len(B), , drop = FALSE]

          lin_s3 <- .act_summarize_stage3_ensemble(
            ensemble = ens_lin,
            target_daily = target,
            obs_col = obs_col,
            interval = stage3_interval
          )
          lin_s3$ensemble <- ens_lin

          # All four Stage-2 deterministic methods are propagated through the SAME
          # paired QMB source-year ensemble. HARD90, FIXED94 and ACT therefore differ
          # from QMB-LIN only by their locked Stage-2 centre shift.
          hard90_s3 <- .act_recenter_adaptive_ensemble(
            qmb_ensemble = ens_lin, qmb_stage2_daily = qmb$daily,
            adaptive_stage2_daily = adaptive, target_daily = target,
            obs_col = obs_col, interval = stage3_interval,
            qmb_col = "qmb_kokic5q_stage2", adaptive_col = "qmb_hard90_stage2",
            clip_zero = TRUE
          )
          fixed94_s3 <- .act_recenter_adaptive_ensemble(
            qmb_ensemble = ens_lin, qmb_stage2_daily = qmb$daily,
            adaptive_stage2_daily = adaptive, target_daily = target,
            obs_col = obs_col, interval = stage3_interval,
            qmb_col = "qmb_kokic5q_stage2", adaptive_col = "qmb_fixed94_stage2",
            clip_zero = TRUE
          )
          act_s3 <- .act_recenter_adaptive_ensemble(
            qmb_ensemble = ens_lin, qmb_stage2_daily = qmb$daily,
            adaptive_stage2_daily = adaptive, target_daily = target,
            obs_col = obs_col, interval = stage3_interval,
            qmb_col = "qmb_kokic5q_stage2", adaptive_col = "qmb_adaptive_stage2",
            clip_zero = TRUE
          )
          s3_objects <- list(
            QMB_LIN_STAGE3 = lin_s3,
            HARD90_STAGE3 = hard90_s3,
            FIXED94_STAGE3 = fixed94_s3,
            ADAPTIVE_CV_STAGE3 = act_s3
          )

          # Metadata + accuracy, interval and probabilistic-calibration rows.
          for (meth in names(s3_objects)) {
            obj <- s3_objects[[meth]]

            aa <- obj$accuracy
            aa$REGIME <- reg; aa$LOCATION_ID <- loc
            aa$GCM <- gg; aa$SSP <- ssp
            aa$METHOD <- meth; aa$N_BOOT <- B
            s3a <- s3a + 1L
            stage3_accuracy_all[[s3a]] <- aa

            im <- obj$interval_metrics
            im$REGIME <- reg; im$LOCATION_ID <- loc
            im$GCM <- gg; im$SSP <- ssp
            im$METHOD <- meth; im$N_BOOT <- B
            s3i <- s3i + 1L
            stage3_interval_all[[s3i]] <- im

            mm3 <- .act_stage3_monthly_metrics(obj$daily_interval)
            mm3$REGIME <- reg; mm3$LOCATION_ID <- loc
            mm3$GCM <- gg; mm3$SSP <- ssp
            mm3$METHOD <- meth; mm3$N_BOOT <- B
            mm3$CV_CUTOFF <- if (meth == "ADAPTIVE_CV_STAGE3") {
              cut_loc$SELECTED_CUTOFF[match(mm3$MONTH, cut_loc$MONTH)]
            } else NA_real_
            s3m <- s3m + 1L
            stage3_month_all[[s3m]] <- mm3

            cal <- .act_probabilistic_calibration(obj$ensemble, target, obs_col = obs_col)
            cal$REGIME <- reg; cal$LOCATION_ID <- loc
            cal$GCM <- gg; cal$SSP <- ssp
            cal$METHOD <- meth; cal$N_BOOT <- B
            s3c <- s3c + 1L
            stage3_calibration_all[[s3c]] <- cal
          }

          # Daily intervals: save every B only if requested. By default save
          # only max(B), which is usually what is used as the final paper series.
          if (isTRUE(save_stage3_daily_all_b) || B == max_b) {
            bdir <- file.path(
              output_root, "FINAL_STAGE3", reg, loc, gg, ssp,
              sprintf("B_%04d", B)
            )
            dir.create(bdir, recursive = TRUE, showWarnings = FALSE)

            utils::write.csv(
              lin_s3$daily_interval,
              file.path(bdir, "qmb_lin_stage3_daily_interval.csv"),
              row.names = FALSE
            )
            utils::write.csv(hard90_s3$daily_interval,
              file.path(bdir, "hard90_stage3_daily_interval.csv"), row.names = FALSE)
            utils::write.csv(fixed94_s3$daily_interval,
              file.path(bdir, "fixed94_stage3_daily_interval.csv"), row.names = FALSE)
            utils::write.csv(act_s3$daily_interval,
              file.path(bdir, "adaptive_cv_stage3_daily_interval.csv"), row.names = FALSE)
          }

          # Do not retain the recentered ensemble between B values.
          hard90_s3$ensemble <- NULL; fixed94_s3$ensemble <- NULL; act_s3$ensemble <- NULL
          rm(ens_lin, lin_s3, hard90_s3, fixed94_s3, act_s3, s3_objects)
          invisible(gc(FALSE))
        }

        rm(ens_max, qmb_s3_max)
        invisible(gc(FALSE))
      }
    }
  }

  final_daily <- .act_safe_rbind(final_daily_all)
  final_month <- .act_safe_rbind(final_month_all)
  utils::write.csv(final_month, file.path(output_root, "FINAL_MONTHLY_METRICS_5METHODS.csv"), row.names = FALSE)

  # Overall comparison computed from the complete final daily collection.
  overall <- list(); oo <- 0L
  method_cols <- c(
    QMB = "qmb_kokic5q_stage2",
    HARD90 = "qmb_hard90_stage2",
    FIXED94 = "qmb_fixed94_stage2",
    ADAPTIVE_CV = "qmb_adaptive_stage2",
    QDM = "qdm_full"
  )
  for (mn in names(method_cols)) {
    col <- method_cols[[mn]]
    oo <- oo + 1L
    z <- final_month[final_month$METHOD == mn, , drop = FALSE]
    overall[[oo]] <- data.frame(
      METHOD = mn,
      DAILY_MAE = .act_mae(final_daily$obs, final_daily[[col]]),
      DAILY_RMSE = .act_rmse(final_daily$obs, final_daily[[col]]),
      DAILY_BIAS = .act_bias(final_daily$obs, final_daily[[col]]),
      MEAN_MONTHLY_P95_ABS_ERROR = mean(z$P95_ABS_ERROR, na.rm = TRUE),
      MEAN_MONTHLY_P99_ABS_ERROR = mean(z$P99_ABS_ERROR, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  overall <- .act_safe_rbind(overall)
  utils::write.csv(overall, file.path(output_root, "FINAL_OVERALL_COMPARISON.csv"), row.names = FALSE)

  # ---------------------------------------------------------------------------
  # D. Final Stage-3 master outputs and fair QMB-LIN vs Adaptive comparison.
  # ---------------------------------------------------------------------------
  stage3_accuracy <- data.frame()
  stage3_interval_metrics <- data.frame()
  stage3_monthly <- data.frame()
  stage3_calibration <- data.frame()
  stage3_overall <- data.frame()

  if (isTRUE(run_stage3_bootstrap)) {
    stage3_accuracy <- .act_safe_rbind(stage3_accuracy_all)
    stage3_interval_metrics <- .act_safe_rbind(stage3_interval_all)
    stage3_monthly <- .act_safe_rbind(stage3_month_all)
    stage3_calibration <- .act_safe_rbind(stage3_calibration_all)

    stage3_overall <- .act_stage3_overall_from_rows(
      stage3_accuracy, stage3_monthly, stage3_interval_metrics
    )

    # Add deterministic full-QDM as reference for each B. QDM has no bootstrap
    # interval; only accuracy/tail columns are meaningful.
    qdm_month <- final_month[final_month$METHOD == "QDM", , drop = FALSE]
    qdm_rows <- vector("list", length(b_values))
    for (ib in seq_along(b_values)) {
      B <- b_values[ib]
      qdm_rows[[ib]] <- data.frame(
        METHOD = "QDM_DETERMINISTIC",
        N_BOOT = B,
        DAILY_MAE = .act_mae(final_daily$obs, final_daily$qdm_full),
        DAILY_RMSE = .act_rmse(final_daily$obs, final_daily$qdm_full),
        DAILY_BIAS = .act_bias(final_daily$obs, final_daily$qdm_full),
        MEAN_MONTHLY_P95_ABS_ERROR =
          mean(qdm_month$P95_ABS_ERROR, na.rm = TRUE),
        MEAN_MONTHLY_P99_ABS_ERROR =
          mean(qdm_month$P99_ABS_ERROR, na.rm = TRUE),
        COVERAGE = NA_real_,
        MEAN_INTERVAL_WIDTH = NA_real_,
        stringsAsFactors = FALSE
      )
    }
    stage3_overall <- rbind(stage3_overall, .act_safe_rbind(qdm_rows))
    stage3_overall <- stage3_overall[
      order(stage3_overall$N_BOOT, stage3_overall$METHOD),
      , drop = FALSE
    ]

    utils::write.csv(
      stage3_accuracy,
      file.path(output_root, "FINAL_STAGE3_ACCURACY_BY_COMBINATION.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      stage3_interval_metrics,
      file.path(output_root, "FINAL_STAGE3_INTERVAL_METRICS_BY_COMBINATION.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      stage3_monthly,
      file.path(output_root, "FINAL_STAGE3_MONTHLY_METRICS.csv"),
      row.names = FALSE
    )
    utils::write.csv(stage3_calibration,
      file.path(output_root, "MASTER_STAGE3_PROBABILISTIC_CALIBRATION.csv"), row.names = FALSE)
    utils::write.csv(
      stage3_overall,
      file.path(output_root, "FINAL_STAGE3_OVERALL_COMPARISON.csv"),
      row.names = FALSE
    )
  }

  # Paper-facing cutoff map.
  cutoff_map <- reshape(
    selected[, c("LOCATION_ID", "MONTH", "SELECTED_CUTOFF")],
    idvar = "LOCATION_ID", timevar = "MONTH", direction = "wide"
  )
  names(cutoff_map) <- sub("SELECTED_CUTOFF\\.", "M", names(cutoff_map))
  utils::write.csv(cutoff_map, file.path(output_root, "PAPER_CUTOFF_MAP_6LOCATION_X_12MONTH.csv"), row.names = FALSE)

  if (isTRUE(save_rds)) {
    saveRDS(
      list(
        method_id = QMB_ACT_METHOD_ID,
        version = QMB_ACT_VERSION,
        manifest = m,
        cv_detail = cv_detail,
        cv_audit = cv_audit,
        selected_cutoff = selected,
        selected_cutoff_strict_05mm = selected_strict,
        final_monthly_metrics = final_month,
        final_overall_comparison = overall,
        stage3_accuracy = stage3_accuracy,
        stage3_interval_metrics = stage3_interval_metrics,
        stage3_monthly_metrics = stage3_monthly,
        stage3_probabilistic_calibration = stage3_calibration,
        stage3_overall_comparison = stage3_overall
      ),
      file.path(output_root, "QMB_ADAPTIVE_QDM_TAIL_MASTER.rds")
    )
  }

  invisible(list(
    manifest = m,
    cv_detail = cv_detail,
    cv_audit = cv_audit,
    selected_cutoff = selected,
    selected_cutoff_strict_05mm = selected_strict,
    final_monthly_metrics = final_month,
    final_overall_comparison = overall,
    stage3_accuracy = stage3_accuracy,
    stage3_interval_metrics = stage3_interval_metrics,
    stage3_monthly_metrics = stage3_monthly,
    stage3_probabilistic_calibration = stage3_calibration,
    stage3_overall_comparison = stage3_overall,
    output_root = normalizePath(output_root, winslash = "/", mustWork = FALSE)
  ))
}

# =============================================================================
# END CORE
# =============================================================================
