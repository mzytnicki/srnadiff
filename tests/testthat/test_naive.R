library(srnadiff)
library(testthat)

context("Checking naive strategy")

exp    <- sRNADiffExample()
exp    <- setStrategies(exp, TRUE, TRUE, TRUE, TRUE)
ranges <- runAllNaive(exp)

test_that("Running naive method", {
    expect_equal(length(ranges), 32)
})
