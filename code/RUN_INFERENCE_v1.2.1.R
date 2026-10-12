# =============================================================================
# RUN_INFERENCE_v1.2.1.R
# ACT-QMB major-revision inferential/common-comparison analysis
# Release: 2026-09-30
# Metric dependency: COMMON_METRIC_CORE_v1.1.0.R
# =============================================================================
# Requires RUN_FULL_ANALYSIS_v1.2.4.R to finish first.
# Design:
#   Evaluation layers: STAGE2 and STAGE3_MEDIAN at B=30,100,250,500,1000,2000
#   Six matched arms per layer: RAW_GCM, QDM, QMB_LIN, HARD90, FIXED94, ACT_QMB
#   RAW/QDM are deterministic benchmarks repeated unchanged at each Stage-3 B.
#   Only QMB_LIN/HARD90/FIXED94/ACT_QMB are bootstrapped in Stage 3.
#   Inference: paired ACT-minus-comparator differences with hierarchical bootstrap:
#      regime-GCM clusters -> retain all SSPs -> moving 2-year evaluation-year blocks.
# =============================================================================
options(stringsAsFactors=FALSE)
.get_script <- function(){a<-commandArgs(FALSE);h<-grep('^--file=',a,value=TRUE);if(length(h))return(normalizePath(sub('^--file=','',h[1]),winslash='/',mustWork=FALSE));NA_character_}
S <- .get_script(); SD <- if(is.na(S)) normalizePath(getwd(),winslash='/',mustWork=FALSE) else dirname(S)
ROOT <- if(tolower(basename(SD))=='code') dirname(SD) else SD
FULL <- file.path(ROOT,'results','full_analysis')
DATA <- file.path(ROOT,'data','adaptive')
OUT <- file.path(ROOT,'results','major_revision_inference'); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
COMMON_METRIC_FILE <- file.path(ROOT,'code','COMMON_METRIC_CORE_v1.1.0.R')
if(!file.exists(COMMON_METRIC_FILE)) stop('Missing ',COMMON_METRIC_FILE,
  '. Place COMMON_METRIC_CORE_v1.1.0.R in code/ before running inference.')
source(COMMON_METRIC_FILE, local=FALSE)
if(!exists('COMMON_METRIC_CORE_VERSION') ||
   !identical(COMMON_METRIC_CORE_VERSION,'v1.1.0'))
  stop('RUN_INFERENCE_v1.2.1.R requires COMMON_METRIC_CORE_v1.1.0.R')

B_VALUES <- CM_B_VALUES
N_BOOT_CI <- 2000L; SEED <- 20260928L; BLOCK_YEARS <- 2L
NI_MARGINS_MM <- c(0.5,1.0,1.5)
WET_THRESHOLD <- CM_WET_THRESHOLD

# v1.1.3 execution controls ------------------------------------------------------
REUSE_EXISTING <- TRUE
FORCE_REBUILD_DESCRIPTIVE <- TRUE
CHECKPOINT_EVERY <- 100L
PROGRESS_EVERY <- 25L
STAGE3_PILOT_R <- 10L
STRICT_EXPECTED_SSP <- 4L
# IMPORTANT: retain the v1.1.3 checkpoint directory intentionally.
# v1.1.4 is checkpoint-compatible with v1.1.3, so completed bootstrap
# replicates (especially BOOT_STAGE2.rds) must not be recomputed.
CHECKPOINT_DIR <- file.path(OUT, "checkpoints_v1.2.1")
dir.create(CHECKPOINT_DIR, recursive=TRUE, showWarnings=FALSE)

cat("\nRUN_INFERENCE_v1.2.1 — all-B Stage-3 discovery fix / audited inference\n")
cat("Bootstrap replicates :", N_BOOT_CI, "\n")
cat("Moving-year block    :", BLOCK_YEARS, "years\n")
cat("Reuse existing       :", REUSE_EXISTING, "\n\n")
rbind_fill <- function(z){z<-z[vapply(z,function(x)is.data.frame(x)&&nrow(x)>0,logical(1))];if(!length(z))return(data.frame());n<-unique(unlist(lapply(z,names)));do.call(rbind,lapply(z,function(x){for(k in setdiff(n,names(x)))x[[k]]<-NA;x[,n,drop=FALSE]}))}
num <- function(x)suppressWarnings(as.numeric(as.character(x)))
canon_gcm <- function(x){
  y<-toupper(gsub("[- .]","_",as.character(x))); y<-gsub("_+","_",y)
  out<-rep(NA_character_,length(y))
  out[grepl("CESM2_WACCM",y,fixed=TRUE)]<-"CESM2_WACCM"
  out[grepl("CMCC_ESM2",y,fixed=TRUE)]<-"CMCC_ESM2"
  out[grepl("MRI_ESM2_0",y,fixed=TRUE)]<-"MRI_ESM2_0"
  out[grepl("TAIESM1",y,fixed=TRUE)]<-"TaiESM1"
  out
}
default_location_id <- function(regime,gcm){paste(regime,ifelse(canon_gcm(gcm)=="MRI_ESM2_0",2L,1L))}
parse_meta <- function(path){
  z<-gsub("\\\\","/",path); zu<-toupper(gsub("[- ]","_",z))
  g<-if(grepl("CESM2_WACCM",zu,fixed=TRUE))"CESM2_WACCM" else if(grepl("CMCC_ESM2",zu,fixed=TRUE))"CMCC_ESM2" else if(grepl("MRI_ESM2_0",zu,fixed=TRUE))"MRI_ESM2_0" else if(grepl("TAIESM1",zu,fixed=TRUE))"TaiESM1" else NA_character_
  s<-regmatches(zu,regexpr("SSP(126|245|370|585)",zu,perl=TRUE)); if(!length(s)||!nzchar(s))s<-NA_character_ else s<-tolower(s)
  r<-if(grepl("MONSOONAL",zu))"Monsoonal" else if(grepl("EQUATORIAL",zu))"Equatorial" else if(grepl("LOCAL",zu))"Local" else NA_character_
  loc<-regmatches(z,regexpr("(Monsoonal|Equatorial|Local)[ _-]*[12]",z,ignore.case=TRUE,perl=TRUE))
  if(!length(loc)||!nzchar(loc))loc<-default_location_id(r,g) else {loc<-gsub("[_-]+"," ",loc);loc<-trimws(loc)}
  data.frame(REGIME=r,LOCATION_ID=loc,GCM=g,SSP=s,stringsAsFactors=FALSE)
}
# ---- Locked metric adapters: COMMON_METRIC_CORE_v1.1.0 is authoritative ----
cm_small_daily <- function(date,o,p){
  d <- data.frame(
    date=as.Date(date), obs=num(o), prediction=num(p),
    stringsAsFactors=FALSE
  )
  keep <- !is.na(d$date) & format(d$date,'%m-%d')!='02-29' &
    as.integer(format(d$date,'%Y'))>=CM_TARGET_YEARS[1] &
    as.integer(format(d$date,'%Y'))<=CM_TARGET_YEARS[2]
  d <- d[keep,,drop=FALSE]
  d$year <- as.integer(format(d$date,'%Y'))
  d$month <- as.integer(format(d$date,'%m'))
  d
}
metric_daily <- function(o,p){
  c(DAILY_MAE=cm_mae(o,p),
    DAILY_RMSE=cm_rmse(o,p),
    DAILY_BIAS=cm_bias(o,p),
    DAILY_CORRELATION=cm_cor(o,p),
    DAILY_KGE=cm_kge(o,p))
}
monthly_total_mae <- function(date,o,p){
  d <- cm_small_daily(date,o,p)
  z <- cm_monthly_totals(d,min_completeness=0.90)
  ok <- z$valid_month & is.finite(z$Total_obs) & is.finite(z$Total_prediction)
  cm_mae(z$Total_obs[ok],z$Total_prediction[ok])
}
wet_occ_mae <- function(date,o,p){
  d <- cm_small_daily(date,o,p)
  z <- cm_monthly_occurrence(d,wet_threshold=CM_WET_THRESHOLD)
  cm_mae(z$WetProp_obs,z$WetProp_prediction)
}
tail_by_month <- function(date,o,p){
  d <- cm_small_daily(date,o,p)
  z <- cm_pooled_calendar_extremes(d,probs=CM_EXTREME_PROBS)
  data.frame(
    MONTH=z$month,
    P90_AE=z$P90_abs_error,
    P95_AE=z$P95_abs_error,
    P99_AE=z$P99_abs_error,
    stringsAsFactors=FALSE
  )
}

# ---- Global paired-difference helper ---------------------------------------------
# v1.1.4 bug fix:
# In v1.1.3 this helper was defined only inside the descriptive BUILD branch.
# Therefore it did not exist when the four descriptive CSVs were REUSED, causing
# an error after a completed bootstrap. It is now defined unconditionally.
pair_metric_table <- function(df, metric, cp) {
  a <- df[df$METHOD == 'ACT_QMB', , drop=FALSE]
  b <- df[df$METHOD == cp, , drop=FALSE]

  by <- intersect(
    c('REGIME','LOCATION_ID','GCM','SSP',
      'EVALUATION_LAYER','N_BOOT','MONTH'),
    names(df)
  )

  z <- merge(a, b, by=by, suffixes=c('_ACT','_CMP'))

  act_col <- paste0(metric, '_ACT')
  cmp_col <- paste0(metric, '_CMP')

  if (!all(c(act_col, cmp_col) %in% names(z))) {
    stop(
      "pair_metric_table(): metric columns not found for metric=", metric,
      ", comparator=", cp, call.=FALSE
    )
  }

  data.frame(
    z[, by, drop=FALSE],
    COMPARATOR=cp,
    METRIC=metric,
    ACT=num(z[[act_col]]),
    COMPARATOR_VALUE=num(z[[cmp_col]]),
    DIFF_ACT_MINUS_COMPARATOR=
      num(z[[act_col]]) - num(z[[cmp_col]]),
    stringsAsFactors=FALSE
  )
}

# ---- v1.2.0 fail-fast canonical audit -------------------------------------------
audit_stage2_48 <- function(d){
  q<-d[d$EVALUATION_LAYER=="STAGE2" & is.na(d$N_BOOT),,drop=FALSE]
  q$GCM<-canon_gcm(q$GCM)
  expected<-c("RAW_GCM","QDM","QMB_LIN","HARD90","FIXED94","ACT_QMB")
  cat("\nCANONICAL BUILD AUDIT\n")
  for(me in expected){
    x<-q[q$METHOD==me,,drop=FALSE]
    cells<-unique(x[,c("REGIME","GCM","SSP")])
    key<-interaction(x$REGIME,x$GCM,x$SSP,drop=TRUE,lex.order=TRUE)
    nd<-if(nrow(x))vapply(split(x$date,key),length,integer(1))else integer()
    ntai<-sum(cells$GCM=="TaiESM1",na.rm=TRUE)
    cat(sprintf("  %-10s : cells=%d | TaiESM1=%d | days=%s\n",me,nrow(cells),ntai,if(length(nd))paste(range(nd),collapse="-")else"NA"))
    if(nrow(cells)!=48L||ntai!=12L||length(nd)!=48L||any(nd!=4015L))stop("CANONICAL BUILD AUDIT FAILED for ",me,". Bootstrap blocked.",call.=FALSE)
  }
  cat("CANONICAL BUILD AUDIT: PASS — 6 methods x 48 cells x 4015 days\n\n")
  invisible(TRUE)
}

# ---- Resume/reuse descriptive layer ---------------------------------------------
canonical_file <- file.path(OUT,'CANONICAL_DAILY_6ARM_ALL_LAYERS.csv')
cell_file <- file.path(OUT,'CELL_METRICS_6ARM_ALL_LAYERS.csv')
month_file <- file.path(OUT,'MONTHLY_TAIL_6ARM_ALL_LAYERS.csv')
pairs_file <- file.path(OUT,'PAIRED_DIFFERENCES_ACT_VS_COMPARATORS.csv')

descriptive_ready <- REUSE_EXISTING && !FORCE_REBUILD_DESCRIPTIVE &&
  all(file.exists(c(canonical_file,cell_file,month_file,pairs_file)))

if (descriptive_ready) {
  cat("[1/6] Canonical daily ............... REUSE\n")
  ALL <- read.csv(canonical_file, check.names=FALSE)
  ALL$date <- as.Date(ALL$date)
  ALL <- cm_prepare_major_revision_schema(ALL)

  cat("[2/6] Cell/monthly metrics ........... REUSE\n")
  CELL <- read.csv(cell_file, check.names=FALSE)
  MONTH <- read.csv(month_file, check.names=FALSE)

  cat("[3/6] Paired differences ............. REUSE\n")
  PAIRS <- read.csv(pairs_file, check.names=FALSE)

  expected <- c('RAW_GCM','QDM','QMB_LIN','HARD90','FIXED94','ACT_QMB')
  if(!all(expected %in% unique(ALL$METHOD)))
    stop('Reused canonical dataset is not a complete six-arm dataset.')
  cm_assert_benchmark_invariance(ALL)
} else {
  cat("[1/6] Building canonical daily dataset ...\n")
  # ---- Load Stage 2 five corrected/reference series --------------------------------
  s2files<-list.files(file.path(FULL,'FINAL_STAGE2'),pattern='stage2_daily_qmb_adaptive_qdm_tail\\.csv$',recursive=TRUE,full.names=TRUE)
  if(!length(s2files))stop('No FINAL_STAGE2 daily files found. Run RUN_FULL_ANALYSIS_v1.2.4.R first.')
  s2meta<-rbind_fill(lapply(s2files,function(f){m<-parse_meta(f);m$FILE<-f;m}))
  cat("\nSOURCE FILE AUDIT\n  Stage-2 files :",nrow(s2meta),"\n  TaiESM1 files :",sum(s2meta$GCM=="TaiESM1",na.rm=TRUE),"\n")
  if(nrow(s2meta)!=48L||sum(s2meta$GCM=="TaiESM1",na.rm=TRUE)!=12L||any(is.na(s2meta$GCM)))stop("SOURCE FILE AUDIT FAILED: expected 48 Stage-2 files including 12 TaiESM1.",call.=FALSE)
  S2<-rbind_fill(lapply(s2files,function(f){x<-read.csv(f,check.names=FALSE);m<-parse_meta(f);map<-c(QMB_LIN='qmb_kokic5q_stage2',HARD90='qmb_hard90_stage2',FIXED94='qmb_fixed94_stage2',ACT_QMB='qmb_adaptive_stage2',QDM='qdm_full');miss<-setdiff(c('date','obs',unname(map)),names(x));if(length(miss))stop('Missing Stage2 columns in ',f,': ',paste(miss,collapse=', '));rbind_fill(lapply(names(map),function(me)data.frame(m,METHOD=me,EVALUATION_LAYER='STAGE2',N_BOOT=NA_integer_,IS_BOOTSTRAPPED=FALSE,date=x$date,obs=x$obs,prediction=x[[map[me]]])))}))
  S2$IS_BOOTSTRAPPED[S2$METHOD=='QMB_LIN']<-FALSE # Stage-2 itself is deterministic.

  # ---- RAW from the exact paired inputs used by the analysis ----------------------
  infiles<-list.files(DATA,pattern='\\.csv$',recursive=TRUE,full.names=TRUE); infiles<-infiles[grepl('paired',basename(infiles),ignore.case=TRUE)]
  RAW<-rbind_fill(lapply(infiles,function(f){x<-read.csv(f,check.names=FALSE); if(!all(c('date','chirps','gcm')%in%names(x)))return(data.frame());d<-as.Date(x$date);keep<-as.integer(format(d,'%Y'))>=2015 & as.integer(format(d,'%Y'))<=2025 & format(d,'%m-%d')!='02-29';m<-parse_meta(f);data.frame(m,METHOD='RAW_GCM',EVALUATION_LAYER='STAGE2',N_BOOT=NA_integer_,IS_BOOTSTRAPPED=FALSE,date=as.character(d[keep]),obs=num(x$chirps[keep]),prediction=num(x$gcm[keep]))}))
  # Stage2 has six arms.
  S2ALL<-rbind_fill(list(RAW,S2))

  # ---- Stage 3 four stochastic methods --------------------------------------------
  # IMPORTANT v1.2.1 FIX:
  # list.files(pattern=...) applies the regex to the BASENAME, not to the full
  # recursive path. The previous runner incorrectly placed "B_0030/" inside
  # pattern=, so all four stochastic methods were discovered as zero files.
  # We now discover by filename first, then filter the normalized FULL PATH by B.
  fnmap<-c(
    QMB_LIN='qmb_lin_stage3_daily_interval.csv',
    HARD90='hard90_stage3_daily_interval.csv',
    FIXED94='fixed94_stage3_daily_interval.csv',
    ACT_QMB='adaptive_cv_stage3_daily_interval.csv'
  )

  stage3_root <- file.path(FULL,'FINAL_STAGE3')
  all_s3_daily <- list.files(
    stage3_root,
    pattern='_stage3_daily_interval\\.csv$',
    recursive=TRUE,
    full.names=TRUE
  )
  all_s3_daily_norm <- gsub('\\\\','/',all_s3_daily)

  # Fail-fast source audit BEFORE constructing the canonical Stage-3 table.
  cat('\nSTAGE-3 SOURCE DISCOVERY AUDIT\n')
  s3_audit <- data.frame()
  s3_parts <- list()
  kk <- 0L

  for(B in B_VALUES){
    btoken <- paste0('/B_',sprintf('%04d',B),'/')
    for(me in names(fnmap)){
      fs <- all_s3_daily[
        grepl(btoken,all_s3_daily_norm,fixed=TRUE) &
        basename(all_s3_daily)==unname(fnmap[me])
      ]

      meta_fs <- if(length(fs)) rbind_fill(lapply(fs,parse_meta)) else data.frame()
      n_cells <- if(nrow(meta_fs)) nrow(unique(meta_fs[,c('REGIME','GCM','SSP')])) else 0L
      n_tai <- if(nrow(meta_fs)) sum(unique(meta_fs[,c('REGIME','GCM','SSP')])$GCM=='TaiESM1',na.rm=TRUE) else 0L

      cat(sprintf('  B=%-4d %-8s : files=%2d | cells=%2d | TaiESM1=%d %s\n',
                  B,me,length(fs),n_cells,n_tai,
                  if(length(fs)==48L && n_cells==48L && n_tai==12L) 'PASS' else 'FAIL'))

      s3_audit <- rbind(
        s3_audit,
        data.frame(N_BOOT=B,METHOD=me,N_FILES=length(fs),
                   N_CELLS=n_cells,N_TAIESM1_CELLS=n_tai,
                   stringsAsFactors=FALSE)
      )

      if(length(fs)!=48L || n_cells!=48L || n_tai!=12L){
        stop(
          'STAGE-3 SOURCE DISCOVERY AUDIT FAILED for B=',B,' method=',me,
          ': expected 48 daily files / 48 cells / 12 TaiESM1 cells; found ',
          length(fs),' / ',n_cells,' / ',n_tai,
          '. Confirm RUN_FULL_ANALYSIS_v1.2.4.R finished with STAGE-3 DAILY ALL-B AUDIT: PASS.',
          call.=FALSE
        )
      }

      kk <- kk + 1L
      s3_parts[[kk]] <- rbind_fill(lapply(fs,function(f){
        x<-read.csv(f,check.names=FALSE)
        req<-c('date','obs','median','lower','upper')
        miss<-setdiff(req,names(x))
        if(length(miss)) stop('Missing Stage-3 columns in ',f,': ',paste(miss,collapse=', '),call.=FALSE)
        if(nrow(x)!=4015L) stop('Stage-3 file does not contain 4015 evaluation days: ',f,call.=FALSE)
        m<-parse_meta(f)
        data.frame(
          m,METHOD=me,EVALUATION_LAYER='STAGE3_MEDIAN',N_BOOT=B,
          IS_BOOTSTRAPPED=TRUE,date=x$date,obs=x$obs,
          prediction=x$median,lower=x$lower,upper=x$upper
        )
      }))
    }
  }

  write.csv(s3_audit,file.path(OUT,'STAGE3_SOURCE_DISCOVERY_AUDIT.csv'),row.names=FALSE)
  S3STOCH <- rbind_fill(s3_parts)
  cat('STAGE-3 SOURCE DISCOVERY AUDIT: PASS\n\n')
  # Repeat RAW and QDM unchanged at every Stage-3 B strictly for matched benchmarking.
  bench0<-S2ALL[S2ALL$METHOD%in%c('RAW_GCM','QDM'),c('REGIME','LOCATION_ID','GCM','SSP','METHOD','date','obs','prediction')]
  S3BENCH<-rbind_fill(lapply(B_VALUES,function(B){z<-bench0;z$EVALUATION_LAYER<-'STAGE3_MEDIAN';z$N_BOOT<-B;z$IS_BOOTSTRAPPED<-FALSE;z}))
  ALL<-rbind_fill(list(S2ALL,S3STOCH,S3BENCH))
  ALL$date<-as.Date(ALL$date)
  ALL$GCM<-canon_gcm(ALL$GCM)
  ALL <- cm_prepare_major_revision_schema(ALL)
  cm_assert_benchmark_invariance(ALL)
  audit_stage2_48(ALL)
  # Audit invariants: six methods per evaluation layer; RAW/QDM invariant over B.
  expected<-c('RAW_GCM','QDM','QMB_LIN','HARD90','FIXED94','ACT_QMB')
  if(!all(expected%in%unique(ALL$METHOD)))stop('Six-arm canonical dataset incomplete.')
  for(me in c('RAW_GCM','QDM'))for(B in B_VALUES){a<-ALL[ALL$METHOD==me & ALL$EVALUATION_LAYER=='STAGE2',c('REGIME','GCM','SSP','date','prediction')];b<-ALL[ALL$METHOD==me & ALL$EVALUATION_LAYER=='STAGE3_MEDIAN' & ALL$N_BOOT==B,c('REGIME','GCM','SSP','date','prediction')];names(a)[5]<-'a';names(b)[5]<-'b';z<-merge(a,b,by=c('REGIME','GCM','SSP','date'));if(nrow(z)==0||any(abs(z$a-z$b)>1e-12,na.rm=TRUE))stop(me,' deterministic benchmark changed at B=',B)}
  write.csv(ALL,file.path(OUT,'CANONICAL_DAILY_6ARM_ALL_LAYERS.csv'),row.names=FALSE)

  # ---- Descriptive metrics ---------------------------------------------------------
  grpkey<-interaction(ALL$REGIME,ALL$LOCATION_ID,ALL$GCM,ALL$SSP,ALL$METHOD,ALL$EVALUATION_LAYER,ifelse(is.na(ALL$N_BOOT),'NA',ALL$N_BOOT),drop=TRUE,lex.order=TRUE)
  CELL<-rbind_fill(lapply(split(seq_len(nrow(ALL)),grpkey),function(ix){d<-ALL[ix,];id<-d[1,c('REGIME','LOCATION_ID','GCM','SSP','METHOD','EVALUATION_LAYER','N_BOOT','IS_BOOTSTRAPPED')];o<-num(d$obs);p<-num(d$prediction);mm<-metric_daily(o,p);data.frame(id,
    DAILY_MAE=mm['DAILY_MAE'],
    DAILY_RMSE=mm['DAILY_RMSE'],
    DAILY_BIAS=mm['DAILY_BIAS'],
    DAILY_CORRELATION=mm['DAILY_CORRELATION'],
    DAILY_KGE=mm['DAILY_KGE'],
    MONTHLY_TOTAL_MAE=monthly_total_mae(d$date,o,p),
    WET_OCC_MAE=wet_occ_mae(d$date,o,p))}))
  MONTH<-rbind_fill(lapply(split(seq_len(nrow(ALL)),grpkey),function(ix){d<-ALL[ix,];id<-d[1,c('REGIME','LOCATION_ID','GCM','SSP','METHOD','EVALUATION_LAYER','N_BOOT','IS_BOOTSTRAPPED')];mm<-tail_by_month(d$date,num(d$obs),num(d$prediction));if(nrow(mm))cbind(id[rep(1,nrow(mm)),],mm)else data.frame()}))
  write.csv(CELL,file.path(OUT,'CELL_METRICS_6ARM_ALL_LAYERS.csv'),row.names=FALSE);write.csv(MONTH,file.path(OUT,'MONTHLY_TAIL_6ARM_ALL_LAYERS.csv'),row.names=FALSE)

  # ---- Paired descriptive differences ACT - comparator ----------------------------
  PAIRS<-rbind_fill(c(lapply(c('DAILY_MAE','DAILY_RMSE','MONTHLY_TOTAL_MAE','WET_OCC_MAE'),function(m)rbind_fill(lapply(c('FIXED94','HARD90','QDM','QMB_LIN','RAW_GCM'),function(cp)pair_metric_table(CELL,m,cp)))),lapply(c('P90_AE','P95_AE','P99_AE'),function(m)rbind_fill(lapply(c('FIXED94','HARD90','QDM','QMB_LIN','RAW_GCM'),function(cp)pair_metric_table(MONTH,m,cp))))))
  write.csv(PAIRS,file.path(OUT,'PAIRED_DIFFERENCES_ACT_VS_COMPARATORS.csv'),row.names=FALSE)

  cat("[1/6] Canonical daily ............... BUILT\n")
  cat("[2/6] Cell/monthly metrics ........... BUILT\n")
  cat("[3/6] Paired differences ............. BUILT\n")
}

# ---- Optimized hierarchical moving-year-block bootstrap -------------------------
# Scientific design is unchanged from v1.1.2:
#   1) resample regime-GCM clusters with replacement;
#   2) retain all four SSPs in each selected cluster;
#   3) within each SSP, resample evaluation years in moving 2-year blocks;
#   4) recalculate metrics from daily values.
#
# Optimization in v1.1.3:
#   - one hierarchy/year-block draw is shared by ALL comparators and ALL metrics
#     within a layer, instead of repeating the same resampling separately for
#     each comparator x metric combination;
#   - daily data are split/indexed once per layer;
#   - checkpoint is written every CHECKPOINT_EVERY replicates and can be resumed.

make_year_sequence <- function(years,L=2L){
  u <- sort(unique(years))
  starts <- u[u <= max(u)-L+1L]
  out <- integer()
  while(length(out) < length(u)){
    st <- sample(starts,1)
    out <- c(out,seq.int(st,st+L-1L))
  }
  out[seq_along(u)]
}

# Re-index sampled source years to unique pseudo-year slots. This preserves
# repeated sampled years as separate temporal units for monthly calculations.
sample_slot_index <- function(year_vec,L=2L){
  ys <- make_year_sequence(year_vec,L)
  idx <- unlist(lapply(ys,function(y) which(year_vec==y)),use.names=FALSE)
  slot <- unlist(lapply(seq_along(ys),function(k) rep.int(k,sum(year_vec==ys[k]))),
                 use.names=FALSE)
  list(idx=idx,slot=slot,source_year=ys)
}

calc_metric_bundle <- function(date,o,p){
  o <- num(o); p <- num(p)
  ans <- c(
    DAILY_MAE=cm_mae(o,p),
    DAILY_RMSE=cm_rmse(o,p)
  )

  d <- data.frame(date=as.Date(date),obs=o,prediction=p,
                  stringsAsFactors=FALSE)
  d$year <- as.integer(format(d$date,'%Y'))
  d$month <- as.integer(format(d$date,'%m'))

  mt <- cm_monthly_totals(d,min_completeness=0.90)
  ok <- mt$valid_month & is.finite(mt$Total_obs) & is.finite(mt$Total_prediction)
  ans['MONTHLY_TOTAL_MAE'] <- cm_mae(mt$Total_obs[ok],mt$Total_prediction[ok])

  wo <- cm_monthly_occurrence(d,wet_threshold=CM_WET_THRESHOLD)
  ans['WET_OCC_MAE'] <- cm_mae(wo$WetProp_obs,wo$WetProp_prediction)

  ex <- cm_pooled_calendar_extremes(d,probs=CM_EXTREME_PROBS)
  ans['P90_AE'] <- mean(ex$P90_abs_error,na.rm=TRUE)
  ans['P95_AE'] <- mean(ex$P95_abs_error,na.rm=TRUE)
  ans['P99_AE'] <- mean(ex$P99_abs_error,na.rm=TRUE)
  ans
}

comparators <- c('FIXED94','HARD90','QDM','QMB_LIN','RAW_GCM')
metrics <- c('DAILY_MAE','DAILY_RMSE','MONTHLY_TOTAL_MAE','WET_OCC_MAE',
             'P90_AE','P95_AE','P99_AE')
layers <- rbind(
  data.frame(EVALUATION_LAYER='STAGE2',N_BOOT=NA_integer_),
  data.frame(EVALUATION_LAYER='STAGE3_MEDIAN',N_BOOT=B_VALUES)
)

layer_tag <- function(layer,B)
  if(is.na(B)) "STAGE2" else paste0("STAGE3_B",sprintf("%04d",as.integer(B)))

prepare_layer <- function(layer,B){
  z <- ALL[
    ALL$EVALUATION_LAYER==layer &
      (if(is.na(B)) is.na(ALL$N_BOOT) else ALL$N_BOOT==B),
    c('REGIME','LOCATION_ID','GCM','SSP','METHOD','date','obs','prediction'),
    drop=FALSE
  ]
  if(!nrow(z)) stop("No canonical data for ",layer," B=",B)

  z$date <- as.Date(z$date)
  z$REGIME <- as.character(z$REGIME)
  z$LOCATION_ID <- as.character(z$LOCATION_ID)
  z$GCM <- as.character(z$GCM)
  z$SSP <- as.character(z$SSP)
  z$METHOD <- as.character(z$METHOD)

  expected_methods <- c('ACT_QMB',comparators)
  z <- z[z$METHOD %in% expected_methods,,drop=FALSE]

  # ---------------------------------------------------------------------------
  # v1.1.6 canonical pairing rule
  #
  # LOCATION_ID is audited, but is NOT used as a merge key. In this project the
  # physical grid/location is deterministically assigned by GCM (MRI uses the
  # paired location-2 grid; the other GCMs use location-1). Earlier runners may
  # encode the same physical pairing with different LOCATION_ID labels. Using
  # LOCATION_ID as a merge key therefore silently discarded MRI and reduced the
  # paired support from 12 to 9 regime-GCM clusters.
  #
  # Safety: before dropping LOCATION_ID from the merge key we require exactly
  # one location label per METHOD x REGIME x GCM x SSP and later verify that
  # observations are identical across all six methods for every paired date.
  # ---------------------------------------------------------------------------
  loc_key <- interaction(z$METHOD,z$REGIME,z$GCM,z$SSP,
                         drop=TRUE,lex.order=TRUE)
  loc_n <- vapply(split(z$LOCATION_ID,loc_key),
                  function(v) length(unique(v[!is.na(v)])), integer(1))
  if(any(loc_n != 1L)){
    bad <- names(loc_n)[loc_n != 1L]
    stop(
      "Canonical pairing audit failed: more than one LOCATION_ID occurs within ",
      "METHOD x REGIME x GCM x SSP for: ",
      paste(head(bad,20),collapse=", "),
      call.=FALSE
    )
  }

  # Audit design support BEFORE six-arm matching.
  support <- unique(z[,c('METHOD','REGIME','GCM','SSP')])
  support_n <- aggregate(SSP ~ METHOD + REGIME + GCM,
                         data=support,
                         FUN=function(v) length(unique(v)))
  if(any(support_n$SSP != STRICT_EXPECTED_SSP)){
    bad <- support_n[support_n$SSP != STRICT_EXPECTED_SSP,,drop=FALSE]
    stop(
      "Canonical support audit failed: every method must contain four SSPs ",
      "within each available regime-GCM cell. First failures: ",
      paste(apply(head(bad,12),1,paste,collapse="|"),collapse=", "),
      call.=FALSE
    )
  }

  method_clusters <- lapply(expected_methods,function(me)
    sort(unique(paste(z$REGIME[z$METHOD==me],z$GCM[z$METHOD==me],sep='|'))))
  names(method_clusters) <- expected_methods
  nmc <- vapply(method_clusters,length,integer(1))
  if(any(nmc != 12L)){
    detail <- paste(names(nmc),nmc,sep="=",collapse=", ")
    stop(
      "Canonical support audit failed before matching: expected 12 regime-GCM ",
      "clusters for EACH method; found ",detail,
      ". This indicates a missing method/GCM input and must not be bypassed.",
      call.=FALSE
    )
  }

  canonical_id <- c('REGIME','GCM','SSP','date')

  # Build a six-arm wide table without LOCATION_ID as a matching key.
  # Carry method-specific observation columns first so exact pairing can be
  # audited rather than assumed.
  wide <- NULL
  for(me in expected_methods){
    q <- z[z$METHOD==me,
           c(canonical_id,'LOCATION_ID','obs','prediction'),
           drop=FALSE]
    # Exactly one row per canonical method/date is required.
    dup <- duplicated(q[,canonical_id]) | duplicated(q[,canonical_id],fromLast=TRUE)
    if(any(dup)){
      stop("Canonical pairing audit failed: duplicate ",me,
           " rows for REGIME/GCM/SSP/date.",call.=FALSE)
    }
    names(q)[names(q)=='LOCATION_ID'] <- paste0('LOCATION_ID__',me)
    names(q)[names(q)=='obs'] <- paste0('OBS__',me)
    names(q)[names(q)=='prediction'] <- me
    wide <- if(is.null(wide)) q else
      merge(wide,q,by=canonical_id,all=FALSE,sort=FALSE)
  }

  if(!nrow(wide))
    stop("Canonical six-arm pairing produced zero matched daily rows.",call.=FALSE)

  # Verify observations are numerically identical across all six arms.
  obs_cols <- paste0('OBS__',expected_methods)
  ref_obs <- num(wide[[obs_cols[1]]])
  obs_ok <- vapply(obs_cols,function(nm){
    x <- num(wide[[nm]])
    all((is.na(x)&is.na(ref_obs)) |
          (is.finite(x)&is.finite(ref_obs)&abs(x-ref_obs)<=1e-10))
  },logical(1))
  if(!all(obs_ok)){
    stop(
      "Canonical pairing audit failed: CHIRPS observations differ across methods ",
      "after REGIME/GCM/SSP/date matching. LOCATION_ID cannot safely be ignored. ",
      "Failing methods: ",paste(expected_methods[!obs_ok],collapse=", "),
      call.=FALSE
    )
  }

  # Verify every canonical experimental cell has complete evaluation dates.
  cell_key <- interaction(wide$REGIME,wide$GCM,wide$SSP,
                          drop=TRUE,lex.order=TRUE)
  nd <- vapply(split(wide$date,cell_key),length,integer(1))
  if(length(nd) != 48L || any(nd != 4015L)){
    stop(
      "Canonical pairing audit failed: expected 48 cells x 4015 evaluation days; ",
      "found ",length(nd)," cells; day-count range ",
      min(nd),"–",max(nd),".",
      call.=FALSE
    )
  }

  wide$obs <- ref_obs
  wide$YEAR <- as.integer(format(wide$date,'%Y'))
  wide$MD <- format(wide$date,'%m-%d')
  wide$CLUSTER <- as.character(paste(wide$REGIME,wide$GCM,sep='|'))

  cat("      canonical pairing audit: 48 cells x 4015 days = PASS\n")
  cat("      six-arm support audit: 12 regime-GCM clusters x 4 SSP = PASS\n")

  # Keep only fields required downstream. LOCATION_ID has already been audited.
  keep <- c('REGIME','GCM','SSP','date','obs',
            expected_methods,'YEAR','MD','CLUSTER')
  wide[,keep,drop=FALSE]
}

# Calculate original point estimate with equal weighting:
# SSP within cluster, then regime-GCM clusters.
point_estimates <- function(wide){
  clusters <- unique(wide$CLUSTER)
  arr <- array(NA_real_,
    dim=c(length(clusters),length(comparators),length(metrics)),
    dimnames=list(clusters,comparators,metrics))
  for(cl in clusters){
    x <- wide[wide$CLUSTER==cl,,drop=FALSE]
    ssps <- unique(x$SSP)
    sarr <- array(NA_real_,
      dim=c(length(ssps),length(comparators),length(metrics)),
      dimnames=list(ssps,comparators,metrics))
    for(s in ssps){
      q <- x[x$SSP==s,,drop=FALSE]
      ma <- calc_metric_bundle(q$date,q$obs,q$ACT_QMB)
      for(cp in comparators){
        mc <- calc_metric_bundle(q$date,q$obs,q[[cp]])
        sarr[s,cp,] <- ma[metrics]-mc[metrics]
      }
    }
    arr[cl,,] <- apply(sarr,c(2,3),mean,na.rm=TRUE)
  }
  apply(arr,c(2,3),mean,na.rm=TRUE)
}

bootstrap_layer <- function(layer,B,R=N_BOOT_CI){
  tag <- layer_tag(layer,B)
  ck <- file.path(CHECKPOINT_DIR,paste0("BOOT_",tag,".rds"))
  cat("\n[4/6] Hierarchical bootstrap — ",tag,"\n",sep="")

  wide <- prepare_layer(layer,B)
  clusters <- unique(as.character(wide$CLUSTER))
  split_cs <- split(
    wide,
    interaction(as.character(wide$CLUSTER), as.character(wide$SSP),
                drop=TRUE, lex.order=TRUE)
  )
  cluster_ssp_keys <- lapply(clusters,function(cl){
    names(split_cs)[vapply(
      split_cs,
      function(x) as.character(x$CLUSTER[1]) == as.character(cl),
      logical(1)
    )]
  })
  names(cluster_ssp_keys) <- clusters

  # Structural audit before any expensive bootstrap work.
  n_keys <- lengths(cluster_ssp_keys)
  if(any(n_keys==0L)){
    stop(
      "Bootstrap structural audit failed: empty cluster/SSP mapping for: ",
      paste(names(n_keys)[n_keys==0L],collapse=", "),
      call.=FALSE
    )
  }
  if(any(n_keys!=STRICT_EXPECTED_SSP)){
    bad <- paste0(names(n_keys)[n_keys!=STRICT_EXPECTED_SSP],
                  " [nSSP=",n_keys[n_keys!=STRICT_EXPECTED_SSP],"]")
    stop(
      "Bootstrap structural audit failed: expected ",
      STRICT_EXPECTED_SSP," SSP per regime-GCM cluster; found ",
      paste(bad,collapse=", "),
      call.=FALSE
    )
  }
  if(length(clusters)!=12L){
    stop("Bootstrap structural audit failed: expected 12 regime-GCM clusters; found ",
         length(clusters),".",call.=FALSE)
  }
  cat("      bootstrap hierarchy audit: 12 clusters x 4 SSP = PASS\n")

  est <- point_estimates(wide)
  if(any(!is.finite(est))){
    stop("Point-estimate audit failed: non-finite ACT-minus-comparator estimate.",
         call.=FALSE)
  }

  bs <- array(NA_real_,
    dim=c(R,length(comparators),length(metrics)),
    dimnames=list(NULL,comparators,metrics))
  start_r <- 1L

  DESIGN_SIGNATURE <- "V1.2.1_48CELL_ALLB_PATHSAFE_STAGE3"
  if(REUSE_EXISTING && file.exists(ck)){
    old <- tryCatch(readRDS(ck),error=function(e) NULL)
    signature_ok <- if(tag=="STAGE2") TRUE else
      (!is.null(old$design_signature) &&
       identical(old$design_signature,DESIGN_SIGNATURE))
    valid_shape <- is.list(old) && identical(old$tag,tag) &&
      identical(old$R,as.integer(R)) && !is.null(old$bs) &&
      identical(dim(old$bs),dim(bs)) && signature_ok

    if(valid_shape){
      finite_rows <- apply(old$bs,1,function(v) all(is.finite(v)))
      done_idx <- which(finite_rows)

      # A valid checkpoint must be a contiguous prefix of completed replicates.
      contiguous <- !length(done_idx) ||
        identical(done_idx,seq_len(max(done_idx)))

      if(contiguous && length(done_idx)){
        bs <- old$bs
        start_r <- max(done_idx)+1L
        cat("      checkpoint found; completed ",max(done_idx),"/",R,"\n",sep="")
      } else if(length(done_idx)==0L && tag!="STAGE2"){
        cat("      invalid/non-finite Stage-3 checkpoint detected — IGNORE and rebuild\n")
      } else if(!contiguous){
        cat("      non-contiguous checkpoint detected — IGNORE and rebuild\n")
      }
    } else {
      cat("      incompatible checkpoint detected — IGNORE and rebuild\n")
    }
  }
  if(start_r>R){
    cat("      checkpoint already complete — REUSE; constructing CI summary\n")
  } else {
    set.seed(SEED + sum(utf8ToInt(tag)))
    # Advance RNG deterministically when resuming by replaying only the random
    # draws needed to reach the next replicate. This is cheap relative to metrics.
    if(start_r>1L){
      for(rr in seq_len(start_r-1L)){
        sel <- sample(clusters,length(clusters),replace=TRUE)
        for(cl in sel){
          for(key in cluster_ssp_keys[[cl]]){
            q <- split_cs[[key]]
            invisible(make_year_sequence(q$YEAR,BLOCK_YEARS))
          }
        }
      }
    }

    # Stage-3 fail-fast pilot: validate the full resampling/metric path before
    # committing to 2,000 replicates. Pilot values are diagnostic only and are
    # not inserted into the final bootstrap array.
    if(layer=="STAGE3_MEDIAN" && start_r==1L){
      cat("      Stage-3 internal pilot: ",STAGE3_PILOT_R," replicates ... ",sep="")
      pilot_ok <- TRUE
      pilot_t0 <- proc.time()[3]
      for(pr in seq_len(STAGE3_PILOT_R)){
        selp <- sample(clusters,length(clusters),replace=TRUE)
        parr <- array(NA_real_,
          dim=c(length(selp),length(comparators),length(metrics)))
        for(j in seq_along(selp)){
          cl <- selp[j]
          keys <- cluster_ssp_keys[[cl]]
          sarr <- array(NA_real_,
            dim=c(length(keys),length(comparators),length(metrics)))
          for(k in seq_along(keys)){
            q <- split_cs[[keys[k]]]
            smp <- sample_slot_index(q$YEAR,BLOCK_YEARS)
            qq <- q[smp$idx,,drop=FALSE]
            pseudo_year <- CM_TARGET_YEARS[1] + smp$slot - 1L
            pdate <- as.Date(paste0(pseudo_year,"-",qq$MD))
            ma <- calc_metric_bundle(pdate,qq$obs,qq$ACT_QMB)
            for(ci in seq_along(comparators)){
              cp <- comparators[ci]
              mc <- calc_metric_bundle(pdate,qq$obs,qq[[cp]])
              sarr[k,ci,] <- ma[metrics]-mc[metrics]
            }
          }
          parr[j,,] <- apply(sarr,c(2,3),mean,na.rm=TRUE)
        }
        pv <- apply(parr,c(2,3),mean,na.rm=TRUE)
        if(any(!is.finite(pv))){
          pilot_ok <- FALSE
          break
        }
      }
      pilot_elapsed <- proc.time()[3]-pilot_t0
      if(!pilot_ok){
        stop("Stage-3 internal pilot failed: non-finite bootstrap metrics.",
             call.=FALSE)
      }
      cat(sprintf("PASS (%.2f min; projected %.1f min for %d replicates)\n",
                  pilot_elapsed/60,
                  (pilot_elapsed/max(STAGE3_PILOT_R,1))*R/60,R))
      # Reset RNG so the diagnostic pilot does not alter the reproducible
      # production bootstrap sequence.
      set.seed(SEED + sum(utf8ToInt(tag)))
    }

    t0 <- proc.time()[3]
    for(r in seq.int(start_r,R)){
      sel <- sample(clusters,length(clusters),replace=TRUE)
      cl_arr <- array(NA_real_,
        dim=c(length(sel),length(comparators),length(metrics)))

      for(j in seq_along(sel)){
        cl <- sel[j]
        keys <- cluster_ssp_keys[[cl]]
        s_arr <- array(NA_real_,
          dim=c(length(keys),length(comparators),length(metrics)))

        for(k in seq_along(keys)){
          q <- split_cs[[keys[k]]]
          smp <- sample_slot_index(q$YEAR,BLOCK_YEARS)
          qq <- q[smp$idx,,drop=FALSE]

          # pseudo dates retain month-day but make repeated sampled source years
          # distinct temporal slots.
          pseudo_year <- CM_TARGET_YEARS[1] + smp$slot - 1L
          pdate <- as.Date(paste0(pseudo_year,"-",qq$MD))

          ma <- calc_metric_bundle(pdate,qq$obs,qq$ACT_QMB)
          for(ci in seq_along(comparators)){
            cp <- comparators[ci]
            mc <- calc_metric_bundle(pdate,qq$obs,qq[[cp]])
            s_arr[k,ci,] <- ma[metrics]-mc[metrics]
          }
        }
        cl_arr[j,,] <- apply(s_arr,c(2,3),mean,na.rm=TRUE)
      }
      bs[r,,] <- apply(cl_arr,c(2,3),mean,na.rm=TRUE)
      if(any(!is.finite(bs[r,,]))){
        stop("Bootstrap replicate ",r," produced non-finite metrics for ",tag,
             ". Check cluster/SSP/year resampling.",call.=FALSE)
      }

      if(r %% PROGRESS_EVERY==0L || r==R){
        elapsed <- proc.time()[3]-t0
        rate <- (r-start_r+1)/max(elapsed,1e-9)
        eta <- (R-r)/max(rate,1e-9)/60
        cat(sprintf("      %4d/%d  elapsed %.1f min  ETA %.1f min\n",
                    r,R,elapsed/60,eta))
      }
      if(r %% CHECKPOINT_EVERY==0L || r==R){
        saveRDS(list(tag=tag,R=as.integer(R),bs=bs,
                     design_signature=DESIGN_SIGNATURE),
                ck,compress=FALSE)
      }
    }
  }

  finite_count <- sum(apply(bs,1,function(v) all(is.finite(v))))
  cat("      finite replicate audit: ",finite_count,"/",R,"\n",sep="")
  if(finite_count!=R){
    stop("Bootstrap finite-replicate audit failed for ",tag,": ",
         finite_count,"/",R," complete finite replicates.",call.=FALSE)
  }

  rows <- list(); n <- 0L
  for(cp in comparators){
    for(me in metrics){
      v <- bs[,cp,me]
      ci <- quantile(v,c(.025,.975),na.rm=TRUE,type=7)

      src <- if(me %in% c('P90_AE','P95_AE','P99_AE')) MONTH else CELL
      pt <- pair_metric_table(src,me,cp)
      pt <- pt[
        pt$EVALUATION_LAYER==layer &
          (if(is.na(B)) is.na(pt$N_BOOT) else pt$N_BOOT==B),
        ,drop=FALSE
      ]
      n <- n+1L
      rows[[n]] <- data.frame(
        EVALUATION_LAYER=layer,N_BOOT=B,COMPARATOR=cp,METRIC=me,
        ESTIMATE=est[cp,me],
        CI_LOW=ci[1],CI_HIGH=ci[2],
        PROP_ACT_BETTER=mean(pt$DIFF_ACT_MINUS_COMPARATOR<0,na.rm=TRUE),
        N_MATCHED_CELLS=nrow(pt),
        BOOT_R=R,BLOCK_YEARS=BLOCK_YEARS,
        stringsAsFactors=FALSE
      )
    }
  }
  rbind_fill(rows)
}

inf_file <- file.path(OUT,'HIERARCHICAL_MOVING_YEAR_BLOCK_BOOTSTRAP_CI.csv')
INF <- rbind_fill(lapply(seq_len(nrow(layers)),function(i){
  bootstrap_layer(layers$EVALUATION_LAYER[i],layers$N_BOOT[i],N_BOOT_CI)
}))
INF$SUPERIOR_95CI <- INF$CI_HIGH < 0
write.csv(INF,inf_file,row.names=FALSE)
if (!nrow(INF) || any(!is.finite(INF$CI_LOW)) || any(!is.finite(INF$CI_HIGH))) {
  stop("Hierarchical-bootstrap CI table failed integrity checks.", call.=FALSE)
}
cat("\n[4/6] Hierarchical bootstrap .......... COMPLETE\n")
cat("      CI rows written: ", nrow(INF), "\n", sep="")

# ---- ACT vs QDM non-inferiority at P95/P99 -------------------------------------
NI0 <- INF[INF$COMPARATOR=='QDM' & INF$METRIC%in%c('P95_AE','P99_AE'),]
NI <- rbind_fill(lapply(NI_MARGINS_MM,function(delta){
  z <- NI0
  z$MARGIN_MM <- delta
  z$NONINFERIOR_95CI <- z$CI_HIGH < delta
  z
}))
write.csv(NI,file.path(OUT,'ACT_VS_QDM_NONINFERIORITY_SENSITIVITY.csv'),
          row.names=FALSE)
cat("[5/6] ACT vs QDM non-inferiority ..... COMPLETE\n")

# ---- Probabilistic calibration --------------------------------------------------
# Only four Stage-3 stochastic methods have predictive intervals. RAW/QDM do not.
calf<-file.path(FULL,'MASTER_STAGE3_PROBABILISTIC_CALIBRATION.csv')
if(file.exists(calf)){
  CAL<-read.csv(calf,check.names=FALSE)
  write.csv(CAL,file.path(OUT,'STAGE3_PROBABILISTIC_CALIBRATION_4METHODS_ALL_B.csv'),row.names=FALSE)
  cat("[6/6] Probabilistic calibration ..... COMPLETE\n")
} else {
  warning("MASTER_STAGE3_PROBABILISTIC_CALIBRATION.csv not found.")
  cat("[6/6] Probabilistic calibration ..... NOT FOUND\n")
}

# ---- Audit manifest -------------------------------------------------------------
aud<-data.frame(
  CHECK=c('Six methods present','RAW invariant across B','QDM invariant across B',
          'Stage3 stochastic methods','B values','Hierarchical CI replicates',
          'Year-block length','Temporal bootstrap date handling',
          'Descriptive daily association metrics','Inference implementation'),
  VALUE=c(paste(expected,collapse=';'),'PASS','PASS',
          'QMB_LIN;HARD90;FIXED94;ACT_QMB',paste(B_VALUES,collapse=';'),
          N_BOOT_CI,BLOCK_YEARS,'UNIQUE_PSEUDO_YEAR_SLOTS',
          'Pearson correlation;KGE','SHARED_LAYER_BOOTSTRAP_WITH_CHECKPOINT_RESUME_V1.1.6_CANONICAL_PAIRING')
)
write.csv(aud,file.path(OUT,'AUDIT_MANIFEST.csv'),row.names=FALSE)

# ---- Metric provenance ----------------------------------------------------------
PROV<-data.frame(
  INFERENCE_VERSION='v1.2.1',
  COMMON_METRIC_CORE=COMMON_METRIC_CORE_VERSION,
  EVALUATION_YEARS=paste(CM_TARGET_YEARS,collapse='-'),
  WET_DAY_DEFINITION=paste0('rainfall > ',CM_WET_THRESHOLD,' mm/day'),
  QUANTILE_TYPE=CM_QUANTILE_TYPE,
  EXTREME_PROBS=paste(CM_EXTREME_PROBS,collapse=';'),
  EXTREME_POOLING='pooled calendar month across 2015-2025; all days including zero/dry days',
  B_VALUES=paste(B_VALUES,collapse=';'),
  stringsAsFactors=FALSE
)
write.csv(PROV,file.path(OUT,'METRIC_PROVENANCE.csv'),row.names=FALSE)

cat('\nRUN_INFERENCE_v1.2.1 COMPLETE\nOutput: ',OUT,'\n',sep='')
