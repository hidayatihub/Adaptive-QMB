# =============================================================================
# PLOT_CODE_v1.0.1.R
# ACT-QMB major-revision publication figure generator
# =============================================================================
#
# PURPOSE
#   Generate reproducible, publication-ready figures directly from the outputs
#   of RUN_FULL_ANALYSIS_v1.2.4.R and RUN_INFERENCE_v1.2.1.R.
#
# OUTPUT FORMATS
#   PDF  : vector
#   TIFF : 600 dpi, LZW compression
#   PNG  : 300 dpi preview
#
# IMPORTANT
#   - No scientific values are hard-coded.
#   - RAW GCM and QDM are deterministic references.
#   - Only QMB_LIN, HARD90, FIXED94 and ACT_QMB have Stage-3 uncertainty.
#   - All figures are generated from actual CSV output.
# =============================================================================

options(stringsAsFactors = FALSE)

PLOT_CODE_VERSION <- "v1.0.1"

# -----------------------------------------------------------------------------
# 0. Resolve repository root
# -----------------------------------------------------------------------------
.get_this_script <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", ca, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[1]),
                                        winslash = "/", mustWork = FALSE))
  nf <- sys.nframe()
  if (nf > 0L) {
    for (i in rev(seq_len(nf))) {
      f <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
      if (!is.null(f) && length(f) == 1L && nzchar(f))
        return(normalizePath(f, winslash = "/", mustWork = FALSE))
    }
  }
  NA_character_
}

THIS_SCRIPT <- .get_this_script()
SCRIPT_DIR <- if (is.na(THIS_SCRIPT)) normalizePath(getwd(), winslash="/",
                                                     mustWork=FALSE) else dirname(THIS_SCRIPT)
REPO_ROOT <- if (tolower(basename(SCRIPT_DIR)) == "code") dirname(SCRIPT_DIR) else SCRIPT_DIR
REPO_ROOT <- normalizePath(REPO_ROOT, winslash="/", mustWork=FALSE)

FULL_DIR <- file.path(REPO_ROOT, "results", "full_analysis")
INF_DIR  <- file.path(REPO_ROOT, "results", "major_revision_inference")
FIG_ROOT <- file.path(REPO_ROOT, "results", "figures")
PDF_DIR  <- file.path(FIG_ROOT, "pdf")
TIFF_DIR <- file.path(FIG_ROOT, "tiff_600dpi")
PNG_DIR  <- file.path(FIG_ROOT, "png_preview")
dir.create(PDF_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TIFF_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(PNG_DIR, recursive=TRUE, showWarnings=FALSE)

# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------
required <- c("ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Install required package(s): ", paste(missing, collapse=", "), call.=FALSE)

library(ggplot2)

has_scales <- requireNamespace("scales", quietly=TRUE)

# -----------------------------------------------------------------------------
# 2. Locked visual identity
# Color-blind-friendly palette; method colors remain identical across figures.
# -----------------------------------------------------------------------------
METHOD_COLORS <- c(
  RAW_GCM = "#6B7280",
  QDM     = "#0072B2",
  QMB_LIN = "#E69F00",
  HARD90  = "#009E73",
  FIXED94 = "#CC79A7",
  ACT_QMB = "#D55E00"
)
METHOD_ORDER <- names(METHOD_COLORS)
METHOD_LABELS <- c(
  RAW_GCM="RAW GCM", QDM="QDM", QMB_LIN="QMB-LIN",
  HARD90="Hard P90", FIXED94="Fixed P94", ACT_QMB="ACT-QMB"
)

theme_act <- function(base_size=11) {
  theme_minimal(base_size=base_size) +
    theme(
      plot.title=element_text(face="bold", size=base_size+2),
      plot.subtitle=element_text(size=base_size),
      axis.title=element_text(face="bold"),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_line(linewidth=0.25, colour="#D9D9D9"),
      strip.text=element_text(face="bold"),
      legend.title=element_blank(),
      legend.position="bottom",
      plot.caption=element_text(size=base_size-2, colour="#555555"),
      plot.margin=margin(8, 10, 8, 8)
    )
}

norm_method <- function(x) {
  z <- toupper(gsub("[^A-Z0-9]+","_",as.character(x)))
  z[z %in% c("RAW","GCM","RAWGCM")] <- "RAW_GCM"
  z[z %in% c("QMB","QMB_KOKIC","QMB_KOKIC5Q","QMB_LIN_STAGE3")] <- "QMB_LIN"
  z[z %in% c("ADAPTIVE_CV","ADAPTIVE_CV_STAGE3","ACT","ACTQMB")] <- "ACT_QMB"
  z[z %in% c("HARD90_STAGE3")] <- "HARD90"
  z[z %in% c("FIXED94_STAGE3")] <- "FIXED94"
  z
}

read_csv0 <- function(path) {
  if (!file.exists(path)) return(NULL)
  utils::read.csv(path, stringsAsFactors=FALSE, check.names=FALSE)
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) return(NA_character_)
  hit[1]
}

save_pub <- function(p, stem, width=8.2, height=5.4) {
  ggsave(file.path(PDF_DIR, paste0(stem,".pdf")), p,
         width=width, height=height, units="in", device="pdf")
  ggsave(file.path(TIFF_DIR, paste0(stem,".tiff")), p,
         width=width, height=height, units="in", dpi=600,
         device="tiff", compression="lzw")
  ggsave(file.path(PNG_DIR, paste0(stem,".png")), p,
         width=width, height=height, units="in", dpi=300,
         device="png")
}

month_labels <- month.abb

# -----------------------------------------------------------------------------
# 3. Input discovery
# -----------------------------------------------------------------------------
cutoff_file <- first_existing(c(
  file.path(FULL_DIR,"SELECTED_CUTOFF_LOCATION_MONTH_ROBUST.csv"),
  file.path(FULL_DIR,"PAPER_CUTOFF_MAP_6LOCATION_X_12MONTH.csv")
))
cell_file <- file.path(INF_DIR,"CELL_METRICS_6ARM_ALL_LAYERS.csv")
tail_file <- file.path(INF_DIR,"MONTHLY_TAIL_6ARM_ALL_LAYERS.csv")
pair_file <- file.path(INF_DIR,"PAIRED_DIFFERENCES_ACT_VS_COMPARATORS.csv")
boot_file <- file.path(INF_DIR,"HIERARCHICAL_MOVING_YEAR_BLOCK_BOOTSTRAP_CI.csv")
ni_file <- file.path(INF_DIR,"ACT_VS_QDM_NONINFERIORITY_SENSITIVITY.csv")
cal_file <- file.path(INF_DIR,"STAGE3_PROBABILISTIC_CALIBRATION_4METHODS_ALL_B.csv")
daily_file <- file.path(INF_DIR,"CANONICAL_DAILY_6ARM_ALL_LAYERS.csv")

CELL <- read_csv0(cell_file)
TAIL <- read_csv0(tail_file)
PAIR <- read_csv0(pair_file)
BOOT <- read_csv0(boot_file)
NI <- read_csv0(ni_file)
CAL <- read_csv0(cal_file)
DAILY <- read_csv0(daily_file)

if (!is.null(CELL)) {
  CELL$METHOD <- norm_method(CELL$METHOD)
  CELL$METHOD <- factor(CELL$METHOD, levels=METHOD_ORDER)
}
if (!is.null(TAIL)) {
  TAIL$METHOD <- norm_method(TAIL$METHOD)
  TAIL$METHOD <- factor(TAIL$METHOD, levels=METHOD_ORDER)
}
if (!is.null(CAL) && "METHOD" %in% names(CAL)) {
  CAL$METHOD <- norm_method(CAL$METHOD)
  CAL$METHOD <- factor(CAL$METHOD, levels=METHOD_ORDER)
}

# -----------------------------------------------------------------------------
# FIGURE 1 — CV-selected cutoff heatmap (6 locations x 12 months)
# -----------------------------------------------------------------------------
if (!is.na(cutoff_file)) {
  C <- read_csv0(cutoff_file)

  if (all(c("LOCATION_ID","MONTH","SELECTED_CUTOFF") %in% names(C))) {
    H <- C
  } else {
    # Wide paper map fallback: M1...M12
    mcols <- grep("^M([1-9]|1[0-2])$", names(C), value=TRUE)
    if (length(mcols)==12L) {
      H <- do.call(rbind, lapply(seq_len(nrow(C)), function(i) {
        data.frame(
          REGIME=if ("REGIME" %in% names(C)) C$REGIME[i] else NA,
          LOCATION_ID=C$LOCATION_ID[i],
          MONTH=as.integer(sub("^M","",mcols)),
          SELECTED_CUTOFF=as.numeric(C[i,mcols]),
          stringsAsFactors=FALSE
        )
      }))
    } else H <- NULL
  }

  if (!is.null(H) && nrow(H)) {
    H$MONTH <- factor(H$MONTH, levels=1:12, labels=month_labels)
    H$LOCATION_ID <- factor(H$LOCATION_ID, levels=unique(H$LOCATION_ID))
    H$label <- paste0("P", format(round(100*as.numeric(H$SELECTED_CUTOFF),1),
                                  nsmall=1, trim=TRUE))
    p <- ggplot(H, aes(MONTH, LOCATION_ID, fill=SELECTED_CUTOFF)) +
      geom_tile(colour="white", linewidth=0.7) +
      geom_text(aes(label=label), size=3.1, fontface="bold") +
      scale_fill_gradient2(low="#2166AC", mid="#FFFFBF", high="#B2182B",
                           midpoint=0.94, name="Cutoff") +
      labs(title="Cross-validated ACT-QMB tail cutoff",
           subtitle="Historical blocked CV; P94 is the fixed-cutoff reference",
           x="Calendar month", y="Climate-regime location",
           caption="Cell labels show the selected percentile cutoff.") +
      theme_act()
    save_pub(p,"Fig_CV_cutoff_heatmap",9.2,4.8)

    H$DELTA_P94 <- H$SELECTED_CUTOFF - 0.94
    p2 <- ggplot(H, aes(MONTH, LOCATION_ID, fill=DELTA_P94)) +
      geom_tile(colour="white", linewidth=0.7) +
      geom_text(aes(label=sprintf("%+.1f",100*DELTA_P94)), size=3.0, fontface="bold") +
      scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",
                           midpoint=0, name=expression(Delta*" cutoff")) +
      labs(title="Departure of CV-selected cutoff from Fixed P94",
           subtitle="Values are percentage-point departures from P94",
           x="Calendar month", y="Climate-regime location") +
      theme_act()
    save_pub(p2,"Fig_CV_cutoff_departure_from_P94",9.2,4.8)
  }
}

# -----------------------------------------------------------------------------
# FIGURE 2 — Six-arm Stage-2 performance boxplots
# -----------------------------------------------------------------------------
if (!is.null(CELL)) {
  S2 <- CELL[toupper(as.character(CELL$EVALUATION_LAYER))=="STAGE2",,drop=FALSE]
  metrics <- intersect(c("DAILY_MAE","DAILY_RMSE","MONTHLY_TOTAL_MAE","WET_OCC_MAE"), names(S2))
  if (nrow(S2) && length(metrics)) {
    L <- do.call(rbind,lapply(metrics,function(m)
      data.frame(METHOD=S2$METHOD, METRIC=m, VALUE=as.numeric(S2[[m]]))))
    L$METRIC <- factor(L$METRIC, levels=metrics,
      labels=c(DAILY_MAE="Daily MAE",DAILY_RMSE="Daily RMSE",
               MONTHLY_TOTAL_MAE="Monthly-total MAE",WET_OCC_MAE="Wet-day occurrence MAE")[metrics])
    p <- ggplot(L,aes(METHOD,VALUE,fill=METHOD)) +
      geom_boxplot(width=.68,outlier.alpha=.35,linewidth=.45) +
      scale_fill_manual(values=METHOD_COLORS,labels=METHOD_LABELS,drop=FALSE) +
      scale_x_discrete(labels=METHOD_LABELS) +
      facet_wrap(~METRIC,scales="free_y",nrow=1) +
      labs(title="Deterministic Stage-2 performance across matched experiments",
           x=NULL,y="Error") +
      theme_act() + theme(axis.text.x=element_text(angle=35,hjust=1),legend.position="none")
    save_pub(p,"Fig_Stage2_six_arm_boxplots",12.5,4.6)
  }
}

# -----------------------------------------------------------------------------
# Supplementary — Stage-2 correlation and KGE
# -----------------------------------------------------------------------------
if (!is.null(CELL) && all(c("DAILY_CORRELATION","DAILY_KGE") %in% names(CELL))) {
  S2a <- CELL[toupper(as.character(CELL$EVALUATION_LAYER))=="STAGE2",,drop=FALSE]
  if(nrow(S2a)) {
    L <- rbind(
      data.frame(METHOD=S2a$METHOD,METRIC="Pearson correlation",VALUE=as.numeric(S2a$DAILY_CORRELATION)),
      data.frame(METHOD=S2a$METHOD,METRIC="KGE",VALUE=as.numeric(S2a$DAILY_KGE))
    )
    p <- ggplot(L,aes(METHOD,VALUE,fill=METHOD)) +
      geom_boxplot(width=.68,outlier.alpha=.35,linewidth=.45) +
      scale_fill_manual(values=METHOD_COLORS,labels=METHOD_LABELS,drop=FALSE) +
      scale_x_discrete(labels=METHOD_LABELS) +
      facet_wrap(~METRIC,scales="free_y",nrow=1) +
      labs(title="Stage-2 association and hydrological efficiency metrics",
           x=NULL,y="Metric value") +
      theme_act() + theme(axis.text.x=element_text(angle=35,hjust=1),legend.position="none")
    save_pub(p,"FigS_Stage2_Correlation_KGE",8.8,4.6)
  }
}

# -----------------------------------------------------------------------------
# FIGURE 3 — Stage-3 median convergence with B
# -----------------------------------------------------------------------------
if (!is.null(CELL)) {
  S3 <- CELL[toupper(as.character(CELL$EVALUATION_LAYER))=="STAGE3_MEDIAN",,drop=FALSE]
  if (nrow(S3) && all(c("N_BOOT","DAILY_MAE") %in% names(S3))) {
    A <- aggregate(DAILY_MAE~METHOD+N_BOOT,S3,mean,na.rm=TRUE)
    p <- ggplot(A,aes(N_BOOT,DAILY_MAE,colour=METHOD,group=METHOD)) +
      geom_line(linewidth=.9) + geom_point(size=2.2) +
      scale_colour_manual(values=METHOD_COLORS,labels=METHOD_LABELS,drop=FALSE) +
      scale_x_continuous(breaks=c(30,100,250,500,1000,2000)) +
      labs(title="Stage-3 median convergence across bootstrap ensemble size",
           subtitle="RAW GCM and QDM remain deterministic matched references",
           x="Bootstrap ensemble size (B)",y="Mean daily MAE (mm/day)") +
      theme_act()
    save_pub(p,"Fig_Stage3_convergence_DailyMAE",8.5,5.3)
  }
}

# -----------------------------------------------------------------------------
# FIGURE 4 — Bulk-tail trade-off / Pareto view at B=2000
# -----------------------------------------------------------------------------
if (!is.null(CELL) && !is.null(TAIL)) {
  S3 <- CELL[toupper(as.character(CELL$EVALUATION_LAYER))=="STAGE3_MEDIAN" &
             suppressWarnings(as.numeric(CELL$N_BOOT))==2000,,drop=FALSE]
  T3 <- TAIL[toupper(as.character(TAIL$EVALUATION_LAYER))=="STAGE3_MEDIAN" &
             suppressWarnings(as.numeric(TAIL$N_BOOT))==2000,,drop=FALSE]
  if (nrow(S3) && nrow(T3) && "P95_AE" %in% names(T3)) {
    a <- aggregate(DAILY_MAE~METHOD,S3,mean,na.rm=TRUE)
    b <- aggregate(P95_AE~METHOD,T3,mean,na.rm=TRUE)
    z <- merge(a,b,by="METHOD")
    p <- ggplot(z,aes(DAILY_MAE,P95_AE,colour=METHOD,label=METHOD)) +
      geom_point(size=4) +
      geom_text(nudge_y=.03,check_overlap=TRUE,size=3.5) +
      scale_colour_manual(values=METHOD_COLORS,labels=METHOD_LABELS,drop=FALSE) +
      labs(title="Bulk–tail performance trade-off at Stage 3 (B = 2000)",
           subtitle="Lower-left indicates smaller daily and P95 errors",
           x="Daily MAE (mm/day)",y="P95 absolute error (mm/day)") +
      theme_act() + theme(legend.position="none")
    save_pub(p,"Fig_Bulk_tail_Pareto_B2000",7.3,5.5)
  }
}

# -----------------------------------------------------------------------------
# FIGURE 5 — Paired ACT differences vs comparators
# -----------------------------------------------------------------------------
if (!is.null(PAIR)) {
  # Flexible discovery of the paired-difference columns.
  diffcols <- grep("(^DELTA_|_DIFF$|DIFFERENCE)", names(PAIR), value=TRUE, ignore.case=TRUE)
  compcol <- intersect(c("COMPARATOR","REFERENCE_METHOD","METHOD_B"),names(PAIR))
  if (length(diffcols) && length(compcol)) {
    keep <- diffcols[seq_len(min(length(diffcols),6L))]
    L <- do.call(rbind,lapply(keep,function(m)
      data.frame(COMPARATOR=PAIR[[compcol[1]]],METRIC=m,DELTA=as.numeric(PAIR[[m]]))))
    p <- ggplot(L,aes(COMPARATOR,DELTA,fill=COMPARATOR)) +
      geom_hline(yintercept=0,linetype=2,colour="#333333") +
      geom_boxplot(outlier.alpha=.25) +
      facet_wrap(~METRIC,scales="free_y") +
      labs(title="Paired performance differences: ACT-QMB minus comparator",
           subtitle="For error metrics, values below zero favor ACT-QMB",
           x=NULL,y=expression(Delta*" metric")) +
      theme_act() + theme(axis.text.x=element_text(angle=35,hjust=1),legend.position="none")
    save_pub(p,"Fig_ACT_paired_differences",11.5,6.4)
  }
}

# -----------------------------------------------------------------------------
# FIGURE 6 — Hierarchical block-bootstrap confidence intervals
# -----------------------------------------------------------------------------
if (!is.null(BOOT)) {
  lo <- intersect(c("CI_LOW","LOWER","CI_LOWER","LCL"),names(BOOT))
  hi <- intersect(c("CI_HIGH","UPPER","CI_UPPER","UCL"),names(BOOT))
  est <- intersect(c("ESTIMATE","MEAN_DIFF","DELTA","POINT_ESTIMATE"),names(BOOT))
  metric <- intersect(c("METRIC","OUTCOME"),names(BOOT))
  comp <- intersect(c("COMPARATOR","REFERENCE_METHOD"),names(BOOT))
  if (length(lo)&&length(hi)&&length(est)&&length(metric)&&length(comp)) {
    B <- BOOT
    B$key <- paste(B[[comp[1]]],B[[metric[1]]],sep=" | ")
    p <- ggplot(B,aes(x=.data[[est[1]]],y=reorder(key,.data[[est[1]]]))) +
      geom_vline(xintercept=0,linetype=2,colour="#555555") +
      geom_errorbarh(aes(xmin=.data[[lo[1]]],xmax=.data[[hi[1]]]),height=.18) +
      geom_point(size=2.4) +
      labs(title="Hierarchical moving-year-block bootstrap confidence intervals",
           subtitle="ACT-QMB minus comparator; negative error differences favor ACT-QMB",
           x="Paired difference with 95% CI",y=NULL) +
      theme_act()
    save_pub(p,"Fig_Hierarchical_bootstrap_CI",9.0,7.0)
  }
}

# -----------------------------------------------------------------------------
# FIGURE 7 — ACT vs QDM non-inferiority
# -----------------------------------------------------------------------------
if (!is.null(NI)) {
  margin <- intersect(c("MARGIN_MM","DELTA_MARGIN","MARGIN"),names(NI))
  hi <- intersect(c("CI_HIGH","CI_UPPER","UCL"),names(NI))
  est <- intersect(c("ESTIMATE","MEAN_DIFF","DELTA","POINT_ESTIMATE"),names(NI))
  metric <- intersect(c("METRIC","OUTCOME"),names(NI))
  if (length(margin)&&length(hi)&&length(est)&&length(metric)) {
    N <- NI
    N$MARGIN <- as.numeric(N[[margin[1]]])
    N$EST <- as.numeric(N[[est[1]]])
    N$HI <- as.numeric(N[[hi[1]]])
    p <- ggplot(N,aes(MARGIN,EST,colour=.data[[metric[1]]],group=.data[[metric[1]]])) +
      geom_line(linewidth=.8) + geom_point(size=2.5) +
      geom_errorbar(aes(ymin=EST,ymax=HI),width=.03) +
      geom_abline(intercept=0,slope=1,linetype=2,colour="#555555") +
      labs(title="Non-inferiority of ACT-QMB versus QDM in the upper tail",
           subtitle="Upper 95% CI below the margin supports non-inferiority",
           x=expression("Non-inferiority margin "*delta*" (mm)"),
           y="ACT-QMB − QDM absolute-error difference (mm)") +
      theme_act()
    save_pub(p,"Fig_ACT_QDM_noninferiority",7.8,5.4)
  }
}

# -----------------------------------------------------------------------------
# FIGURE 8 — Probabilistic calibration vs B
# -----------------------------------------------------------------------------
if (!is.null(CAL) && nrow(CAL)) {
  covcol <- intersect(c("COVERAGE","coverage"),names(CAL))
  bcol <- intersect(c("N_BOOT","B"),names(CAL))
  nomcol <- intersect(c("NOMINAL","NOMINAL_INTERVAL","nominal"),names(CAL))
  strcol <- intersect(c("SCOPE","STRATUM","TAIL_STRATUM","RANGE"),names(CAL))
  if (length(covcol)&&length(bcol)&&length(nomcol)) {
    CAL$COVERAGE_PLOT <- as.numeric(CAL[[covcol[1]]])
    CAL$B_PLOT <- as.numeric(CAL[[bcol[1]]])
    CAL$NOMINAL_PLOT <- as.numeric(CAL[[nomcol[1]]])
    CAL$SCOPE_PLOT <- if(length(strcol)) as.character(CAL[[strcol[1]]]) else "ALL"

    # Full-distribution calibration: show each nominal level separately.
    A <- CAL[CAL$SCOPE_PLOT=="ALL",,drop=FALSE]
    if(nrow(A)) {
      p <- ggplot(A,aes(B_PLOT,COVERAGE_PLOT,colour=METHOD,group=METHOD)) +
        geom_line(linewidth=.85) + geom_point(size=2) +
        geom_abline(intercept=0,slope=1,linetype=2,colour="#444444") +
        scale_colour_manual(values=METHOD_COLORS,labels=METHOD_LABELS,drop=FALSE) +
        scale_x_continuous(breaks=c(30,100,250,500,1000,2000)) +
        coord_cartesian(ylim=c(0,1)) +
        facet_wrap(~NOMINAL_PLOT,labeller=label_both) +
        labs(title="Stage-3 probabilistic calibration across ensemble size",
             subtitle="Each panel compares empirical coverage with its own nominal level",
             x="Bootstrap ensemble size (B)",y="Empirical coverage") +
        theme_act()
      # Add correct horizontal nominal line separately per facet.
      ref <- unique(A[,c("NOMINAL_PLOT"),drop=FALSE])
      p <- p + geom_hline(data=ref,aes(yintercept=NOMINAL_PLOT),
                          inherit.aes=FALSE,linetype=2,colour="#444444")
      save_pub(p,"Fig_Probabilistic_coverage_vs_B",10.0,7.0)
    }

    # Upper-tail calibration uses the locked 95% predictive interval.
    T <- CAL[CAL$SCOPE_PLOT!="ALL" & abs(CAL$NOMINAL_PLOT-.95)<1e-12,,drop=FALSE]
    if(nrow(T)) {
      p2 <- ggplot(T,aes(B_PLOT,COVERAGE_PLOT,colour=METHOD,group=METHOD)) +
        geom_hline(yintercept=.95,linetype=2,colour="#444444") +
        geom_line(linewidth=.85) + geom_point(size=2) +
        scale_colour_manual(values=METHOD_COLORS,labels=METHOD_LABELS,drop=FALSE) +
        scale_x_continuous(breaks=c(30,100,250,500,1000,2000)) +
        coord_cartesian(ylim=c(0,1)) +
        facet_wrap(~SCOPE_PLOT) +
        labs(title="Upper-tail Stage-3 probabilistic calibration",
             subtitle="Dashed line = nominal 95% coverage",
             x="Bootstrap ensemble size (B)",y="Empirical coverage") +
        theme_act()
      save_pub(p2,"Fig_Probabilistic_upper_tail_coverage_vs_B",10.0,6.0)
    }
  }
}

# -----------------------------------------------------------------------------
# FIGURE 9 — Calibration heatmap by regime / method / tail scope
# -----------------------------------------------------------------------------
if (!is.null(CAL) && all(c("REGIME","METHOD") %in% names(CAL))) {
  covcol <- intersect(c("COVERAGE","coverage"),names(CAL))
  strcol <- intersect(c("SCOPE","STRATUM","TAIL_STRATUM","RANGE"),names(CAL))
  nomcol <- intersect(c("NOMINAL","NOMINAL_INTERVAL","nominal"),names(CAL))
  if (length(covcol)&&length(strcol)&&length(nomcol)) {
    Z <- CAL
    Z <- Z[suppressWarnings(as.numeric(Z$N_BOOT))==
           max(suppressWarnings(as.numeric(Z$N_BOOT)),na.rm=TRUE),,drop=FALSE]
    # For a clean comparable heatmap, use the 95% interval in ALL and tail scopes.
    Z <- Z[abs(as.numeric(Z[[nomcol[1]]])-.95)<1e-12,,drop=FALSE]
    if(nrow(Z)) {
      f <- stats::as.formula(paste(covcol[1],"~ REGIME+METHOD+",strcol[1]))
      A <- aggregate(f,Z,mean,na.rm=TRUE)
      A$ERR <- abs(as.numeric(A[[covcol[1]]])-.95)
      p <- ggplot(A,aes(METHOD,REGIME,fill=ERR)) +
        geom_tile(colour="white",linewidth=.6) +
        geom_text(aes(label=sprintf("%.3f",as.numeric(.data[[covcol[1]]]))),size=3) +
        scale_fill_gradient(low="#F7FBFF",high="#CB181D",name="|Coverage−0.95|") +
        scale_x_discrete(labels=METHOD_LABELS) +
        labs(title="Probabilistic calibration by rainfall regime and tail scope",
             subtitle="B = 2000; cell labels are empirical 95% coverage",
             x=NULL,y=NULL) + theme_act() +
        facet_wrap(stats::as.formula(paste("~",strcol[1])))
      save_pub(p,"Fig_Calibration_heatmap_regime_tail",10.5,6.3)
    }
  }
}

# -----------------------------------------------------------------------------
# FIGURE 10 — Interval sharpness vs calibration (B=2000)
# -----------------------------------------------------------------------------
if (!is.null(CAL)) {
  covcol <- intersect(c("COVERAGE","coverage"),names(CAL))
  widthcol <- intersect(c("MEAN_WIDTH","MEAN_INTERVAL_WIDTH","mean_interval_width","WIDTH"),names(CAL))
  nomcol <- intersect(c("NOMINAL","NOMINAL_INTERVAL","nominal"),names(CAL))
  strcol <- intersect(c("SCOPE","STRATUM","TAIL_STRATUM","RANGE"),names(CAL))
  if (length(covcol)&&length(widthcol)&&length(nomcol)&&"N_BOOT"%in%names(CAL)) {
    Z <- CAL[suppressWarnings(as.numeric(CAL$N_BOOT))==2000 &
             abs(as.numeric(CAL[[nomcol[1]]])-.95)<1e-12,,drop=FALSE]
    if(nrow(Z)) {
      p <- ggplot(Z,aes(x=.data[[widthcol[1]]],y=.data[[covcol[1]]],colour=METHOD)) +
        geom_hline(yintercept=.95,linetype=2,colour="#444444") +
        geom_point(size=3,alpha=.8) +
        scale_colour_manual(values=METHOD_COLORS,labels=METHOD_LABELS,drop=FALSE) +
        labs(title="Stage-3 calibration–sharpness trade-off (B = 2000)",
             subtitle="95% predictive intervals; narrower is sharper if calibration is retained",
             x="Mean predictive interval width (mm/day)",y="Empirical 95% coverage") +
        theme_act()
      if(length(strcol)) p <- p + facet_wrap(stats::as.formula(paste("~",strcol[1])))
      save_pub(p,"Fig_Calibration_sharpness_B2000",9.5,6.0)
    }
  }
}

# -----------------------------------------------------------------------------
# Supplementary figures — metric heatmaps by GCM and regime
# -----------------------------------------------------------------------------
if (!is.null(CELL) && all(c("GCM","REGIME","DAILY_MAE") %in% names(CELL))) {
  Z <- CELL[toupper(as.character(CELL$EVALUATION_LAYER))=="STAGE3_MEDIAN" &
            suppressWarnings(as.numeric(CELL$N_BOOT))==2000,,drop=FALSE]
  if(nrow(Z)) {
    A <- aggregate(DAILY_MAE~REGIME+GCM+METHOD,Z,mean,na.rm=TRUE)
    p <- ggplot(A,aes(GCM,REGIME,fill=DAILY_MAE)) +
      geom_tile(colour="white",linewidth=.55) +
      geom_text(aes(label=sprintf("%.2f",DAILY_MAE)),size=2.8) +
      scale_fill_gradient(low="#FFF7EC",high="#7F0000",name="Daily MAE") +
      facet_wrap(~METHOD,labeller=as_labeller(METHOD_LABELS)) +
      labs(title="Daily MAE by climate regime and GCM (Stage 3, B = 2000)",
           x=NULL,y=NULL) +
      theme_act() + theme(axis.text.x=element_text(angle=35,hjust=1))
    save_pub(p,"FigS_DailyMAE_heatmap_regime_GCM",12.0,6.0)
  }
}

# -----------------------------------------------------------------------------
# Manifest
# -----------------------------------------------------------------------------
manifest <- data.frame(
  plot_code_version=PLOT_CODE_VERSION,
  generated_at=as.character(Sys.time()),
  repository_root=REPO_ROOT,
  full_analysis_dir=FULL_DIR,
  inference_dir=INF_DIR,
  figure_root=FIG_ROOT,
  stringsAsFactors=FALSE
)
utils::write.csv(manifest,file.path(FIG_ROOT,"FIGURE_MANIFEST.csv"),row.names=FALSE)

cat("\n============================================================\n")
cat("PLOT_CODE_v1.0.1 completed.\n")
cat("Publication figures:", FIG_ROOT, "\n")
cat("PDF vector          :", PDF_DIR, "\n")
cat("TIFF 600 dpi        :", TIFF_DIR, "\n")
cat("PNG previews        :", PNG_DIR, "\n")
cat("============================================================\n")
