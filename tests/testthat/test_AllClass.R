library(srnadiff)
library(testthat)

context("Checking class constructors")

basedir    <- system.file("extdata", package="srnadiff", mustWork = TRUE)
sampleInfo <- read.csv(file.path(basedir, "dataInfo.csv"))
gtfFile    <- file.path(basedir, "Homo_sapiens.GRCh38.76.gtf.gz")
annotReg   <- readAnnotation(gtfFile, feature="gene", source="miRNA")
bamFiles   <- file.path(basedir, sampleInfo$FileName)

test_that("Running default constructor", {
    object      <- new("srnadiffExp",
                         bamFiles  =bamFiles,
                         sampleInfo=sampleInfo,
                         annotReg  =annotReg)
    expect_is(object, "srnadiffExp")
})

test_that("Running implemented constructor", {
    object <- srnadiffExp(bamFiles, sampleInfo, annotReg)
    expect_is(object, "srnadiffExp")
})

test_that("Running implemented constructor with one condition replicates", {
    expect_error(srnadiffExp(bamFiles, sampleInfo[1:3, ], annotReg))
})

test_that("Running implemented constructor with missing sample information", {
    expect_error(srnadiffExp(bamFiles, sampleInfo[1:2], annotReg))
})

test_that("Running example constructor", {
    object <- srnadiffExample()
    expect_is(object, "srnadiffExp")
})
