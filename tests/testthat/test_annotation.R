library(srnadiff)
library(testthat)

context("Checking annotation import")

dir     <- system.file("extdata", package="srnadiff", mustWork = TRUE)
gtfFile <- file.path(dir, "Homo_sapiens.GRCh38.76.gtf.gz")
gffFile <- file.path(dir, "mirbase21_GRCh38.gff3")

test_that("Importing miRNAs from GTF file", {
    annotation <- readWholeGenomeAnnotation(gtfFile)
    expect_equal(length(annotation), 164)
})

test_that("Importing precursor miRNAs from mirBase GFF3 file", {
    annotation <- readMiRBasePreAnnotation(gffFile)
    expect_equal(length(annotation), 98)
})

test_that("Importing mature miRNAs from mirBase GFF3 file", {
    gtfFile    <- file.path(dir, "mirbase21_GRCh38.gff3")
    annotation <- readMiRBaseMatureAnnotation(gffFile)
    expect_equal(length(annotation), 161)
})

test_that("Importing miRNAs from GTF file manually", {
    gtfFile    <- file.path(dir, "Homo_sapiens.GRCh38.76.gtf.gz")
    annotation <- readAnnotation(gtfFile, source="miRNA", feature="gene",
                                    name="gene_name")
    expect_equal(length(annotation), 164)
})
