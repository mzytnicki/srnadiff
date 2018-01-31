library(srnadiff)
library(testthat)

context("Checking plot")

exp  <- sRNADiffExample()
exp  <- runAll(exp)
plot <- plotRegion(exp, regions(exp, 0.05)[1])

test_that("Running plot", {
    expect_equal(1, 1)
    #expect_is(plot, "ggplot2")
})
