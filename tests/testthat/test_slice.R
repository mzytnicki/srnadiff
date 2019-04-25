library(srnadiff)
library(testthat)

context("Checking slice strategy")

exp    <- srnadiffExample()
parameters(exp) <- srnadiffDefaultParameters
ranges <- runSlice(exp)

test_that("Running IR method", {
    expect_equal(length(ranges), 8)
})
