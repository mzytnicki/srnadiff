library(srnadiff)
library(testthat)

context("Checking slicing strategy")

exp    <- sRNADiffExample()
ranges <- runSlice(exp)

test_that("Running slicing method", {
    expect_equal(length(ranges), 31)
})
