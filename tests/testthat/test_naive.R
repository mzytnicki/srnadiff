library(srnadiff)
library(testthat)

context("Checking naive strategy")

exp     <- sRNADiffExample()
exp     <- setStrategies(exp, FALSE, TRUE, FALSE, FALSE)
exp     <- runAll(exp)
regions <- regions(exp)

test_that("Running naive method", {
    expect_equal(length(regions), 32)
})
