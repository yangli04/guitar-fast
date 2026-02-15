#!/usr/bin/env Rscript
# Multi-group and GRangesList input tests
# Run: micromamba run -n guitar Rscript tests/test_multigroup.R

library(GuitarFast)
library(GenomicFeatures)
library(rtracklayer)

passed <- 0
failed <- 0

test <- function(name, expr) {
  tryCatch({
    result <- eval(expr)
    if (isTRUE(result)) {
      cat(sprintf("  PASS: %s\n", name))
      passed <<- passed + 1
    } else {
      cat(sprintf("  FAIL: %s\n", name))
      failed <<- failed + 1
    }
  }, error = function(e) {
    cat(sprintf("  FAIL: %s (%s)\n", name, e$message))
    failed <<- failed + 1
  })
}

txdb <- loadDb(system.file("extdata", "mm10_toy.sqlite", package = "GuitarFast"))
bed1 <- system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed", package = "GuitarFast")
bed2 <- system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed6.bed", package = "GuitarFast")

cat("== Two BED files, custom names ==\n")
p1 <- GuitarPlot(txTxdb = txdb,
                 stBedFiles = list(bed1, bed2),
                 stGroupName = c("BED12", "BED6"),
                 enableCI = FALSE, pltTxType = c("mrna"))
test("returns ggplot", inherits(p1, "gg"))

cat("\n== GRangesList input ==\n")
gr1 <- blocks(import(bed1))
gr2 <- import(bed2)
p2 <- GuitarPlot(txTxdb = txdb,
                 stGRangeLists = list(gr1, gr2),
                 stGroupName = c("Peaks1", "Peaks2"),
                 enableCI = FALSE, pltTxType = c("tx"))
test("returns ggplot", inherits(p2, "gg"))

cat("\n== Single GRangesList with CI ==\n")
p3 <- GuitarPlot(txTxdb = txdb,
                 stGRangeLists = list(gr1),
                 enableCI = TRUE, pltTxType = c("tx"),
                 CI_ResamplingTime = 50)
test("returns ggplot", inherits(p3, "gg"))

cat("\n== mapFilterTranscript = TRUE ==\n")
p4 <- GuitarPlot(txTxdb = txdb,
                 stBedFiles = list(bed1, bed2),
                 stGroupName = c("A", "B"),
                 mapFilterTranscript = TRUE,
                 enableCI = FALSE, pltTxType = c("mrna"))
test("returns ggplot", inherits(p4, "gg"))

cat(sprintf("\n== Results: %d passed, %d failed ==\n", passed, failed))
if (failed > 0) quit(status = 1)
