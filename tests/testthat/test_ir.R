library(srnadiff)
library(testthat)

context("Checking IR strategy")

exp    <- sRNADiffExample()
ranges <- runAllIR(exp)

test_that("Running IR method", {
    expect_equal(length(ranges), 0)
})
