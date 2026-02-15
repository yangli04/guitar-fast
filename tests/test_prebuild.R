#!/usr/bin/env Rscript
# Test GuitarTxdb save/load prebuild functionality
# Run: micromamba run -n guitar Rscript tests/test_prebuild.R

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
save_file <- tempfile(fileext = ".RData")

cat("== Save GuitarTxdb ==\n")
# Build and save via txGuitarTxdbSaveFile
# The save appends prefix, so we use a simpler approach: build + save manually
gt <- makeGuitarTxdb(txdb = txdb, txPrimaryOnly = FALSE)
save(guitarTxdb, list = "guitarTxdb", envir = list2env(list(guitarTxdb = gt)), file = save_file)
test("save file exists", file.exists(save_file))
test("save file > 0 bytes", file.info(save_file)$size > 0)

cat("\n== Load and plot from saved GuitarTxdb ==\n")
p <- GuitarPlot(txGuitarTxdb = save_file,
                stBedFiles = list(bed1),
                enableCI = FALSE, pltTxType = c("tx"))
test("returns ggplot from cached txdb", inherits(p, "gg"))

cat("\n== Load with CI ==\n")
p2 <- GuitarPlot(txGuitarTxdb = save_file,
                 stBedFiles = list(bed1),
                 enableCI = TRUE, pltTxType = c("tx"),
                 CI_ResamplingTime = 50)
test("returns ggplot with CI from cached txdb", inherits(p2, "gg"))

unlink(save_file)

cat(sprintf("\n== Results: %d passed, %d failed ==\n", passed, failed))
if (failed > 0) quit(status = 1)
