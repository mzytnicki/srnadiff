library(srnadiff)
library(testthat)

context("Checking HMM strategy")

object    <- srnadiffExample()
parameters(object) <- srnadiffDefaultParameters
counts    <- buildDataHmm(object)
pvalues   <- computePvalues(object, counts, 1)
intervals <- hmm(object, counts, pvalues)

test_that("Building data", {
    expect_equal(dim(counts), c(271, 6))
})

test_that("Computing p-values", {
    expect_equal(length(pvalues), 271)
})

test_that("Running HMM", {
    expect_equal(length(intervals), 16)
})
