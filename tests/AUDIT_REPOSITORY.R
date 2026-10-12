cat("Adaptive-QMB repository preflight audit\n")
cat("Working directory:", normalizePath(getwd(), winslash="/", mustWork=FALSE), "\n\n")

required_files <- c(
  "code/QMB_KOKIC_CORE.R",
  "code/QMB_ADAPTIVE_QDMT_CORE_v1.2.1.R",
  "code/COMMON_METRIC_CORE_v1.1.0.R",
  "code/CHRIPS_GCM_CORE_v1.1.0.R",
  "code/RUN_FULL_ANALYSIS_v1.2.4.R",
  "code/RUN_INFERENCE_v1.2.1.R",
  "code/PLOT_CODE_v1.0.1.R",
  "README.md",
  "CITATION.cff",
  "LICENSE_PENDING.md"
)

failed <- FALSE

cat("1. Required files\n")
for (f in required_files) {
  ok <- file.exists(f)
  cat(if (ok) "PASS " else "FAIL ", f, "\n", sep="")
  if (!ok) failed <- TRUE
}

cat("\n2. R syntax\n")
r_files <- c(
  required_files[grepl("\\.R$", required_files)],
  "tests/RUN_ANALYSIS_PILOT_v1.0.1.R"
)
for (f in r_files[file.exists(r_files)]) {
  ok <- tryCatch({parse(file=f); TRUE}, error=function(e) {
    cat("FAIL ", f, ": ", conditionMessage(e), "\n", sep="")
    FALSE
  })
  if (ok) cat("PASS ", f, "\n", sep="") else failed <- TRUE
}

cat("\n3. Citation placeholders\n")
if (file.exists("CITATION.cff")) {
  cff <- readLines("CITATION.cff", warn=FALSE)
  placeholders <- grep("REPLACE|ADD_GITHUB|COMPLETE_BEFORE|ISI_", cff, value=TRUE, ignore.case=TRUE)
  if (length(placeholders)) {
    cat("FAIL unresolved CITATION.cff placeholder(s):\n")
    cat(paste0("  ", placeholders, "\n"))
    failed <- TRUE
  } else {
    cat("PASS no known placeholders\n")
  }
}

cat("\n4. Input CSV count\n")
dat <- list.files("data/adaptive", pattern="\\.csv$", recursive=TRUE, full.names=TRUE)
if (length(dat) == 0L) {
  cat("SKIP no input CSVs bundled (expected for public code repository)\n")
} else if (length(dat) == 48L) {
  cat("PASS 48 paired input CSV files found\n")
} else {
  cat("FAIL found ", length(dat), " CSV files; expected 48 for a reproduction run\n", sep="")
  failed <- TRUE
}

cat("\n5. Repository hygiene\n")

# Check only files actually tracked by Git.
# Local ignored metadata such as .Rproj.user may legitimately exist
# on a developer machine but must never be committed to the repository.
git_available <- nzchar(Sys.which("git"))

if (git_available) {
  
  tracked <- tryCatch(
    system2(
      "git",
      c("ls-files"),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character(0)
  )
  
  forbidden_patterns <- c(
    "(^|/)\\.Rproj\\.user(/|$)",
    "(^|/)\\.DS_Store$",
    "(^|/)__MACOSX(/|$)",
    "(^|/)\\._"
  )
  
  found_bad <- tracked[
    Reduce(
      `|`,
      lapply(
        forbidden_patterns,
        function(p) grepl(p, tracked)
      )
    )
  ]
  
  if (length(found_bad)) {
    cat("FAIL unwanted metadata is tracked by Git:\n")
    cat(paste0("  ", found_bad, "\n"))
    failed <- TRUE
  } else {
    cat("PASS no unwanted RStudio/macOS metadata tracked by Git\n")
  }
  
} else {
  
  cat("SKIP Git executable not available; repository hygiene not checked\n")
}

cat("\n6. Optional plotting dependency\n")
cat("ggplot2:", if(requireNamespace("ggplot2",quietly=TRUE)) "AVAILABLE" else "NOT INSTALLED (needed only for plotting)", "\n")

cat("\nPreflight audit ", if (failed) "FAILED" else "PASSED", ".\n", sep="")
if (failed) stop("Repository preflight audit failed.", call.=FALSE)
