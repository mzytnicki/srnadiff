library(srnadiff)
library(testthat)

context("Checking the p-value methods")

exp             <- srnadiffExample()
parameters(exp) <- srnadiffDefaultParameters

test_that("Running DESeq2 method", {
    exp <- srnadiff(exp, diffMethod = "DESeq2")
    expect_equal(length(regions(exp)), 18)
})

test_that("Running edgeR method", {
    exp <- srnadiff(exp, diffMethod = "edgeR")
    expect_equal(length(regions(exp)), 8)
})

test_that("Running baySeq method", {
    exp <- srnadiff(exp, diffMethod = "baySeq")
    expect_equal(length(regions(exp)), 8)
})

test_that("Running unknown method", {
    expect_error(srnadiff(exp, diffMethod = "unkown"))
})
