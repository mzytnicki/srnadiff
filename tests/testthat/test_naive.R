library(srnadiff)
library(testthat)

context("Checking naive strategy")

exp    <- srnadiffExample()
parameters(exp) <- srnadiffDefaultParameters
ranges <- runNaive(exp)

test_that("Running naive method", {
    expect_equal(length(ranges), 24)
})
