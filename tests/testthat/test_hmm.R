library(srnadiff)
library(testthat)

context("Checking HMM strategy")

exp       <- sRNADiffExample()
counts    <- buildDataHmm(exp)
pvalues   <- computePvalues(exp, counts)
intervals <- runHmm(exp, counts, pvalues)

test_that("Building data", {
    expect_equal(dim(counts), c(447, 6))
})

test_that("Computing p-values", {
    expect_equal(length(pvalues), 447)
})

test_that("Running HMM", {
    expect_equal(length(intervals), 4)
})
