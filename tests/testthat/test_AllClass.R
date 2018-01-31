library(srnadiff)
library(testthat)

context("Checking class constructors")

dir        <- system.file("extdata", package="srnadiff", mustWork = TRUE)
data       <- read.csv(file.path(dir, "data.csv"))
gtfFile    <- file.path(dir, "Homo_sapiens.GRCh38.76.gtf.gz")
annotation <- readWholeGenomeAnnotation(gtfFile)
bamFiles   <- file.path(dir, data$FileName)
replicates <- data$SampleName
conditions <- factor(data$Condition)

test_that("Running default constructor", {
    object      <- new("sRNADiff",
                         annotation       =annotation,
                         bamFileNames     =bamFiles,
                         replicates       =replicates,
                         conditions       =conditions,
                         design           =DataFrame(condition=conditions)
    )
    expect_is(object, "sRNADiff")
})

test_that("Running implemented constructor", {
    exp <- sRNADiffExp(annotation, bamFiles, replicates, conditions)
    expect_is(exp, "sRNADiff")
})

test_that("Running example constructor", {
    exp <- sRNADiffExample()
    expect_is(exp, "sRNADiff")
})
