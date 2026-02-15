#!/usr/bin/env Rscript
# Benchmark test - measures execution time of key operations
# Run: micromamba run -n guitar Rscript tests/test_benchmark.R

library(Guitar)
library(GenomicFeatures)

txdb <- loadDb(system.file("extdata", "mm10_toy.sqlite", package = "Guitar"))
bed1 <- system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed", package = "Guitar")

bench <- function(label, expr) {
  t1 <- proc.time()
  result <- eval(expr)
  elapsed <- (proc.time() - t1)["elapsed"]
  cat(sprintf("  %-35s %6.2fs\n", label, elapsed))
  invisible(result)
}

cat("== Benchmarks ==\n")
bench("makeGuitarTxdb",
  makeGuitarTxdb(txdb = txdb, txPrimaryOnly = FALSE))

bench("GuitarPlot tx, no CI",
  GuitarPlot(txTxdb = txdb, stBedFiles = list(bed1),
             enableCI = FALSE, pltTxType = c("tx")))

bench("GuitarPlot tx, CI (200 resamp)",
  GuitarPlot(txTxdb = txdb, stBedFiles = list(bed1),
             enableCI = TRUE, pltTxType = c("tx")))

bench("GuitarPlot all types, no CI",
  GuitarPlot(txTxdb = txdb, stBedFiles = list(bed1),
             enableCI = FALSE, pltTxType = c("tx","mrna","ncrna"),
             miscOutFilePrefix = tempfile()))

bench("GuitarPlot all types, CI",
  GuitarPlot(txTxdb = txdb, stBedFiles = list(bed1),
             enableCI = TRUE, pltTxType = c("tx","mrna","ncrna"),
             miscOutFilePrefix = tempfile()))

cat("\nDone.\n")
