# GuitarFast

**GuitarFast** is a performance-optimized fork of the [Guitar](https://bioconductor.org/packages/release/bioc/html/Guitar.html) Bioconductor package (v2.26.0). It visualizes the distribution of RNA-related genomic features across transcript landmarks -- TSS, start codon, stop codon, and TES -- with **3-6x faster** execution than the original.

## What It Does

GuitarFast takes genomic feature coordinates (e.g., m6A peaks, RNA binding protein sites, RNA modification sites) and plots their density distribution relative to RNA transcript structure. This reveals whether features cluster near specific landmarks like the stop codon, which is a hallmark of m6A methylation.

Three transcript views are supported:

| View | Components | Use Case |
|:-----|:-----------|:---------|
| **tx** | Promoter - RNA - Tail | All transcripts, general overview |
| **mrna** | Promoter - 5'UTR - CDS - 3'UTR - Tail | Protein-coding transcripts |
| **ncrna** | Promoter - ncRNA - Tail | Long non-coding RNAs |

## Installation

### Prerequisites

GuitarFast depends on Bioconductor packages. Install them first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicFeatures", "rtracklayer",
                       "AnnotationDbi", "GenomicRanges"))
```

### From GitHub (recommended)

```r
# Using remotes (lightweight)
install.packages("remotes")
remotes::install_github("yliuchicago/GuitarFast")

# Or using devtools
# install.packages("devtools")
# devtools::install_github("yliuchicago/GuitarFast")
```

To install with vignettes:

```r
remotes::install_github("yliuchicago/GuitarFast", build_vignettes = TRUE)
```

### From Source

```bash
git clone https://github.com/yliuchicago/GuitarFast.git
R CMD INSTALL GuitarFast
```

### Optional

For faster confidence interval computation:

```r
install.packages("matrixStats")
```

## Vignettes

| Vignette | Description | Links |
|:---------|:------------|:------|
| Quick Start Guide | Installation, basic usage, all input formats | [HTML](https://yangli04.github.io/guitar-fast/Guitar-Quick-Start.html) / [source](vignettes/Guitar-Quick-Start.Rmd) |
| Benchmark | Performance comparison with original Guitar | [HTML](https://yangli04.github.io/guitar-fast/GuitarFast-Benchmark.html) / [source](vignettes/GuitarFast-Benchmark.Rmd) |

## Quick Start

```r
library(GuitarFast)
library(GenomicFeatures)

# Load transcriptome annotation
txdb <- loadDb(system.file("extdata", "mm10_toy.sqlite", package = "GuitarFast"))

# BED file with m6A peaks
bed_file <- system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed",
                         package = "GuitarFast")

# Plot distribution on transcripts
GuitarPlot(txTxdb = txdb,
           stBedFiles = list(bed_file),
           enableCI = FALSE,
           pltTxType = c("tx"))
```

## Usage Examples

### Plot with Confidence Intervals

Bootstrap confidence intervals show statistical reliability of the density curve.

```r
GuitarPlot(txTxdb = txdb,
           stBedFiles = list(bed_file),
           enableCI = TRUE,
           pltTxType = c("tx"))
```

### mRNA-Specific Distribution

View distribution across 5'UTR, CDS, and 3'UTR:

```r
GuitarPlot(txTxdb = txdb,
           stBedFiles = list(bed_file),
           enableCI = FALSE,
           pltTxType = c("mrna"))
```

### Compare Multiple Groups

```r
bed12 <- system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed",
                      package = "GuitarFast")
bed6  <- system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed6.bed",
                      package = "GuitarFast")

GuitarPlot(txTxdb = txdb,
           stBedFiles = list(bed12, bed6),
           stGroupName = c("BED12", "BED6"),
           enableCI = FALSE,
           pltTxType = c("mrna"))
```

### Use GRanges Input Directly

```r
library(rtracklayer)

gr1 <- blocks(import(bed12))
gr2 <- import(bed6)

GuitarPlot(txTxdb = txdb,
           stGRangeLists = list(gr1, gr2),
           stGroupName = c("Peaks-BED12", "Peaks-BED6"),
           enableCI = FALSE,
           pltTxType = c("tx"))
```

### Save Pre-built Transcriptome for Reuse

Building the transcriptome coordinate system takes ~1-2 seconds. Build once and reuse:

```r
# Build and save
guitarTxdb <- makeGuitarTxdb(txdb = txdb, txPrimaryOnly = FALSE)
save(guitarTxdb, file = "my_guitarTxdb.RData")

# Load in later sessions
GuitarPlot(txGuitarTxdb = "my_guitarTxdb.RData",
           stBedFiles = list(bed_file),
           enableCI = TRUE,
           pltTxType = c("tx"))
```

### Use GTF/GFF Annotation

```r
# From GTF file
GuitarPlot(txGTF = "annotation.gtf",
           stBedFiles = list(bed_file),
           enableCI = FALSE)

# From genome version (downloads TxDb automatically)
GuitarPlot(txGenomeVer = "hg38",
           stBedFiles = list(bed_file))
```

### Save PDFs for All Transcript Types

```r
GuitarPlot(txTxdb = txdb,
           stBedFiles = list(bed_file),
           enableCI = TRUE,
           pltTxType = c("tx", "mrna", "ncrna"),
           miscOutFilePrefix = "my_output")
# Saves: my_output_tx_test.pdf, my_output_mrna_test.pdf, my_output_ncrna_test.pdf
```

## Benchmark: GuitarFast vs Guitar

GuitarFast achieves 3-6x speedup over the original Guitar package on all workloads. The results below were measured on the bundled toy dataset (mm10, 1000 m6A peaks):

| Scenario | Guitar (original) | GuitarFast | Speedup |
|:---------|-------------------:|-----------:|--------:|
| tx, no CI | 8.40s | **1.52s** | **5.5x** |
| tx, CI (200 resamples) | 9.37s | **1.56s** | **6.0x** |
| mrna, no CI | 4.22s | **1.40s** | **3.0x** |
| All types, CI (200 resamples) | 18.77s | **4.79s** | **3.9x** |

The density curves produced by GuitarFast are **numerically identical** to the original Guitar for the no-CI path (max |y_orig - y_opt| = 0.00).

### What Was Optimized

1. **`vapply(txRange, NROW, ...)` replaced with `elementNROWS()`** -- single largest bottleneck (6.7s to 0.008s). The original called R's `NROW` on each of ~2700 GRanges elements individually; `elementNROWS()` does it in one C-level call.

2. **Default `CI_ResamplingTime` reduced from 1000 to 200** -- benchmarked as producing statistically equivalent confidence intervals (mean CI range: 0.259 vs 0.264). Users can still set `CI_ResamplingTime = 1000` for the original behavior.

3. **Pre-allocated CI resampling matrix** -- replaced `replicate()` per-iteration allocation with a single pre-allocated 256 x N matrix and batch `sample.int()`.

4. **C-optimized quantile computation** -- uses `matrixStats::rowQuantiles()` when available (falls back gracefully).

5. **Vectorized `normalize()`** -- replaced `apply()` per-row calls with `max.col()` and `cbind` matrix indexing.

6. **Fixed growing GRanges** -- replaced O(n^2) `c(GRanges, GRanges)` in loop with list accumulation + `do.call(c, list)`.

7. **Removed duplicate rectangle annotations** -- visual bug fix where each transcript component rectangle was drawn twice.

8. **Fixed GuitarTxdb save/load** -- original used `read.table()` on binary RData files; now correctly uses `load()`.

## Parameter Reference

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `txTxdb` | NULL | TxDb object for gene annotation |
| `txGTF` / `txGFF` | NULL | Path to GTF/GFF annotation file |
| `txGenomeVer` | NULL | Genome version (hg19, hg38, mm9, mm10) |
| `txGuitarTxdb` | NULL | Path to saved GuitarTxdb .RData file |
| `stBedFiles` | NULL | List of BED file paths |
| `stGRangeLists` | NULL | List of GRanges/GRangesList objects |
| `stGroupName` | NULL | Custom group names for legend |
| `enableCI` | TRUE | Show bootstrap confidence intervals |
| `CI_ResamplingTime` | 200 | Number of bootstrap resamples |
| `pltTxType` | c("tx","mrna","ncrna") | Transcript types to plot |
| `headOrtail` | TRUE | Include 1kb flanking regions |
| `stSampleNum` | 10 | Points sampled per genomic feature |
| `miscOutFilePrefix` | NA | Set to save PDFs (one per txType) |
| `mapFilterTranscript` | TRUE | Filter transcript mappings by length |
| `txPrimaryOnly` | FALSE | Use only primary (longest) transcripts |
| `adjust` | 1 | Bandwidth adjustment for density estimation |

## Input Formats

### Genomic Features (what you want to plot)

- **BED files** (BED6 or BED12) via `stBedFiles`
- **GRanges / GRangesList objects** via `stGRangeLists`

### Transcriptome Annotation (the reference)

- **TxDb object** via `txTxdb` (e.g., from `GenomicFeatures::loadDb()`)
- **GTF file** via `txGTF`
- **GFF3 file** via `txGFF`
- **Genome version string** via `txGenomeVer` (e.g., "hg38", requires corresponding TxDb package)
- **Pre-built GuitarTxdb** via `txGuitarTxdb` (path to .RData file)

## License

GPL-2, same as the original Guitar package. Based on Guitar v2.26.0 by Cui et al. (2016).
