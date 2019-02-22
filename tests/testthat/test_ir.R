library(srnadiff)
library(testthat)

context("Checking IR strategy")

exp    <- sRNADiffExample()
exp    <- computeNormalizationFactors(exp)
exp    <- computeLogFoldChange(exp)
ranges <- runAllIR(exp)

test_that("Running IR method", {
    expect_equal(length(ranges), 7)
})
