#!/usr/bin/env Rscript
# Basic functionality tests for Guitar
# Run: micromamba run -n guitar Rscript tests/test_basic.R

library(GuitarFast)
library(GenomicFeatures)

passed <- 0
failed <- 0

test <- function(name, expr) {
  tryCatch({
    result <- eval(expr)
    if (isTRUE(result)) {
      cat(sprintf("  PASS: %s\n", name))
      passed <<- passed + 1
    } else {
      cat(sprintf("  FAIL: %s (returned %s)\n", name, deparse(result)))
      failed <<- failed + 1
    }
  }, error = function(e) {
    cat(sprintf("  FAIL: %s (%s)\n", name, e$message))
    failed <<- failed + 1
  })
}

# Setup
txdb <- loadDb(system.file("extdata", "mm10_toy.sqlite", package = "GuitarFast"))
bed1 <- system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed", package = "GuitarFast")

cat("== makeGuitarTxdb ==\n")
gt <- makeGuitarTxdb(txdb = txdb, txPrimaryOnly = FALSE)
test("returns a list", is.list(gt))
test("has txTypes", !is.null(gt$txTypes))
test("contains tx", "tx" %in% gt$txTypes)
test("contains mrna", "mrna" %in% gt$txTypes)
test("contains ncrna", "ncrna" %in% gt$txTypes)
test("tx has componentWidth matrix", is.matrix(gt$tx$componentWidth))
test("tx has startPoint matrix", is.matrix(gt$tx$startPoint))
test("mrna componentWidthAverage sums > 0", sum(gt$mrna$componentWidthAverage) > 0)

cat("\n== GuitarPlot no CI ==\n")
p1 <- GuitarPlot(txTxdb = txdb, stBedFiles = list(bed1),
                 enableCI = FALSE, pltTxType = c("tx"))
test("returns ggplot", inherits(p1, "gg"))

cat("\n== GuitarPlot with CI ==\n")
p2 <- GuitarPlot(txTxdb = txdb, stBedFiles = list(bed1),
                 enableCI = TRUE, pltTxType = c("tx"), CI_ResamplingTime = 50)
test("returns ggplot", inherits(p2, "gg"))

cat("\n== GuitarPlot mrna only ==\n")
p3 <- GuitarPlot(txTxdb = txdb, stBedFiles = list(bed1),
                 enableCI = FALSE, pltTxType = c("mrna"))
test("returns ggplot", inherits(p3, "gg"))

cat("\n== GuitarPlot ncrna only ==\n")
p4 <- GuitarPlot(txTxdb = txdb, stBedFiles = list(bed1),
                 enableCI = FALSE, pltTxType = c("ncrna"))
test("returns ggplot", inherits(p4, "gg"))

cat(sprintf("\n== Results: %d passed, %d failed ==\n", passed, failed))
if (failed > 0) quit(status = 1)
