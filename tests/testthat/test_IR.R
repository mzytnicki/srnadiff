library(srnadiff)
library(testthat)

context("Checking IR strategy")

exp    <- srnadiffExample()
parameters(exp) <- srnadiffDefaultParameters
ranges <- runIR(exp)

test_that("Running IR method", {
    expect_equal(length(ranges), 8)
})
