#!/usr/bin/env Rscript
# Side-by-side benchmark: Original vs Optimized Guitar
# Run: micromamba run -n guitar Rscript tests/test_compare.R
#
# Requires Guitar_orig/ directory next to Guitar/

library(ggplot2)

pkg_orig <- normalizePath("../Guitar_orig")
pkg_opt  <- normalizePath(".")
outdir   <- "tests/figures"
dir.create(outdir, showWarnings = FALSE)

cat("============================================================\n")
cat("  Guitar: Original vs Optimized Comparison\n")
cat("============================================================\n\n")

run_bench <- function(label) {
  library(Guitar, quietly = TRUE)
  library(GenomicFeatures, quietly = TRUE)

  txdb_file <- system.file("extdata", "mm10_toy.sqlite", package = "Guitar")
  bed1      <- system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed", package = "Guitar")
  txdb <- loadDb(txdb_file)

  # Warm-up
  invisible(GuitarPlot(txTxdb = txdb, stBedFiles = list(bed1),
                       enableCI = FALSE, pltTxType = c("tx")))

  res <- list()

  # 1: tx no CI (deterministic)
  t0 <- proc.time()
  res$p_tx_noci <- GuitarPlot(txTxdb = txdb, stBedFiles = list(bed1),
                               enableCI = FALSE, pltTxType = c("tx"))
  res$time_tx_noci <- (proc.time() - t0)["elapsed"]

  # 2: tx with CI (200 resamples for fair comparison)
  t0 <- proc.time()
  res$p_tx_ci <- GuitarPlot(txTxdb = txdb, stBedFiles = list(bed1),
                              enableCI = TRUE, pltTxType = c("tx"),
                              CI_ResamplingTime = 200)
  res$time_tx_ci200 <- (proc.time() - t0)["elapsed"]

  # 3: mrna no CI
  t0 <- proc.time()
  res$p_mrna_noci <- GuitarPlot(txTxdb = txdb, stBedFiles = list(bed1),
                                 enableCI = FALSE, pltTxType = c("mrna"))
  res$time_mrna_noci <- (proc.time() - t0)["elapsed"]

  # 4: all types CI
  t0 <- proc.time()
  GuitarPlot(txTxdb = txdb, stBedFiles = list(bed1),
             enableCI = TRUE, pltTxType = c("tx", "mrna", "ncrna"),
             CI_ResamplingTime = 200, miscOutFilePrefix = tempfile())
  res$time_all_ci200 <- (proc.time() - t0)["elapsed"]

  # Extract density data from no-CI plot for numerical comparison
  res$density <- ggplot_build(res$p_tx_noci)$data[[1]]

  detach("package:Guitar", unload = TRUE)
  return(res)
}

# ---- Phase 1: Install & benchmark ORIGINAL ----
cat("--- Installing ORIGINAL package ---\n")
invisible(system2("R", c("CMD", "INSTALL", "--no-multiarch", pkg_orig),
                  stdout = FALSE, stderr = FALSE))
cat("--- Benchmarking ORIGINAL ---\n")
orig <- run_bench("original")

# Save original plots
ggsave(file.path(outdir, "cmp_original_tx_noci.pdf"),  orig$p_tx_noci,   width = 10, height = 6)
ggsave(file.path(outdir, "cmp_original_tx_ci.pdf"),    orig$p_tx_ci,     width = 10, height = 6)
ggsave(file.path(outdir, "cmp_original_mrna_noci.pdf"),orig$p_mrna_noci, width = 10, height = 6)

# ---- Phase 2: Install & benchmark OPTIMIZED ----
cat("\n--- Installing OPTIMIZED package ---\n")
invisible(system2("R", c("CMD", "INSTALL", "--no-multiarch", pkg_opt),
                  stdout = FALSE, stderr = FALSE))
cat("--- Benchmarking OPTIMIZED ---\n")
opt <- run_bench("optimized")

# Save optimized plots
ggsave(file.path(outdir, "cmp_optimized_tx_noci.pdf"),  opt$p_tx_noci,   width = 10, height = 6)
ggsave(file.path(outdir, "cmp_optimized_tx_ci.pdf"),    opt$p_tx_ci,     width = 10, height = 6)
ggsave(file.path(outdir, "cmp_optimized_mrna_noci.pdf"),opt$p_mrna_noci, width = 10, height = 6)

# ---- Results ----
cat("\n============================================================\n")
cat("  SPEED COMPARISON\n")
cat("============================================================\n\n")
cat(sprintf("  %-30s %10s %10s %10s\n", "Scenario", "Original", "Optimized", "Speedup"))
cat(sprintf("  %-30s %10s %10s %10s\n", "--------", "--------", "---------", "-------"))

scenarios <- list(
  list("tx, no CI",           "time_tx_noci"),
  list("tx, CI (200 resamp)", "time_tx_ci200"),
  list("mrna, no CI",         "time_mrna_noci"),
  list("all types, CI (200)", "time_all_ci200")
)
for (s in scenarios) {
  o <- orig[[s[[2]]]]
  n <- opt[[s[[2]]]]
  cat(sprintf("  %-30s %9.2fs %9.2fs %9.1fx\n", s[[1]], o, n, o / n))
}

cat("\n============================================================\n")
cat("  NUMERICAL COMPARISON (tx, no CI - deterministic)\n")
cat("============================================================\n\n")
d_orig <- orig$density
d_opt  <- opt$density
cat(sprintf("  Original:  %d points, y range [%.4f, %.4f]\n",
            nrow(d_orig), min(d_orig$y), max(d_orig$y)))
cat(sprintf("  Optimized: %d points, y range [%.4f, %.4f]\n",
            nrow(d_opt), min(d_opt$y), max(d_opt$y)))
if (nrow(d_orig) == nrow(d_opt)) {
  max_diff <- max(abs(d_orig$y - d_opt$y))
  cat(sprintf("  Max |y_orig - y_opt|: %.2e\n", max_diff))
  if (max_diff < 1e-10) {
    cat("  Result: IDENTICAL density curves\n")
  } else if (max_diff < 0.01) {
    cat("  Result: Negligible difference (< 0.01)\n")
  } else {
    cat("  Result: WARNING - curves differ\n")
  }
}

cat(sprintf("\n  Figures saved to %s/cmp_*.pdf\n", outdir))
cat("  Compare original vs optimized PDFs side by side.\n\n")
