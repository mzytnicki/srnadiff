library(srnadiff)
library(testthat)

context("Checking class constructors")

basedir    <- system.file("extdata", package="srnadiff", mustWork = TRUE)
sampleInfo <- read.csv(file.path(basedir, "dataInfo.csv"))
gtfFile    <- file.path(basedir, "Homo_sapiens.GRCh38.76.gtf.gz")
annotReg   <- readAnnotation(gtfFile, feature="gene", source="miRNA")
bamFiles   <- paste(file.path(basedir, sampleInfo$FileName), "bam", sep = ".")

test_that("Running default constructor", {
    object      <- new("sRNADiffExp",
                         bamFiles  =bamFiles,
                         sampleInfo=sampleInfo,
                         annotReg  =annotReg)
    expect_is(object, "sRNADiff")
})

test_that("Running implemented constructor", {
    srnaExp <- srnadiffExp(bamFiles, sampleInfo, annotReg)
    expect_is(exp, "sRNADiffExp")
})

test_that("Running implemented constructor with unbalanced number of replicates", {
    expect_error(srnadiffExp(bamFiles, sampleInfo[1:3], annotReg))
})

test_that("Running example constructor", {
    exp <- sRNADiffExample()
    expect_is(exp, "sRNADiffExp")
})
