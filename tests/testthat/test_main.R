library(srnadiff)
library(testthat)

context("Checking main functions")

exp <- sRNADiffExample()
exp <- setStrategies(exp, TRUE, TRUE, TRUE, TRUE)
exp <- runAll(exp)

test_that("Testing regions method", {
    expect_equal(length(regions(exp)), 38)
})

test_that("Running with different strategies", {
    exp2 <- sRNADiffExample()
    exp2 <- setStrategies(exp2, TRUE, FALSE, TRUE, TRUE)
    exp2 <- runAll(exp2)
    expect_equal(length(regions(exp2)), 31)
})

test_that("Running with different sizes", {
    exp2 <- sRNADiffExample()
    exp2 <- setStrategies(exp2, TRUE, TRUE, TRUE, TRUE)
    exp2 <- setSizes(exp2, 15, 30)
    exp2 <- runAll(exp2)
    expect_equal(length(regions(exp2)), 65)
})

test_that("Running with different minimum depth", {
    exp2 <- sRNADiffExample()
    exp2 <- setStrategies(exp2, TRUE, TRUE, TRUE, TRUE)
    exp2 <- setMinDepth(exp2, 1)
    exp2 <- runAll(exp2)
    expect_equal(length(regions(exp2)), 38)
})

test_that("Running with different merge distance", {
    exp2 <- sRNADiffExample()
    exp2 <- setStrategies(exp2, FALSE, TRUE, FALSE, FALSE)
    exp2 <- setMergeDistance(exp2, 1)
    exp2 <- runAll(exp2)
    expect_equal(length(regions(exp2)), 33)
})

test_that("Running with different minimum differences", {
    exp2 <- sRNADiffExample()
    exp2 <- setStrategies(exp2, TRUE, TRUE, TRUE, TRUE)
    exp2 <- setMinDifferences(exp2, 100)
    exp2 <- runAll(exp2)
    expect_equal(length(regions(exp2)), 38)
})

test_that("Running with different transition probabilities", {
    exp2 <- sRNADiffExample()
    exp2 <- setStrategies(exp2, FALSE, FALSE, FALSE, TRUE)
    exp2 <- setTransitionProbabilities(exp2, 0.5, 0.5)
    exp2 <- runAll(exp2)
    expect_equal(length(ranges), 1)
})

test_that("Running with different emission probabilities", {
    exp2 <- sRNADiffExample()
    exp2 <- setStrategies(exp2, FALSE, FALSE, FALSE, TRUE)
    exp2 <- setEmissionProbabilities(exp2, 0.75)
    exp2 <- runAll(exp2)
    expect_equal(length(regions(exp2)), 27)
})

test_that("Running with different emission threshold", {
    exp2 <- sRNADiffExample()
    exp2 <- setStrategies(exp2, FALSE, FALSE, FALSE, TRUE)
    exp2 <- setEmissionThreshold(exp2, 0.5)
    exp2 <- runAll(exp2)
    expect_equal(length(regions(exp2)), 27)
})

test_that("Running with different number of overlapping base pairs", {
    exp2 <- sRNADiffExample()
    exp2 <- setStrategies(exp2, TRUE, TRUE, TRUE, TRUE)
    exp2 <- setMinOverlap(exp2, 20)
    exp2 <- runAll(exp2)
    expect_equal(length(regions(exp2)), 17)
})

test_that("Running several threads", {
    exp2 <- sRNADiffExample()
    exp2 <- setStrategies(exp2, TRUE, TRUE, TRUE, TRUE)
    exp2 <- setNThreads(exp2, 4)
    exp2 <- runAll(exp2)
    expect_equal(regions(exp), regions(exp2))
})

test_that("Running main function", {
    expect_equal(length(regions(exp)), 38)
})

test_that("Running redundant removal", {
    regions  <- GRanges(seqnames = "14", strand = "+",
                        ranges = IRanges(start = c(60000000, 60000100),
                        width = 10))
    padj     <- c(0.01,  0.001)
    regions2 <- removeRedundant(regions, padj)
    expect_equal(length(regions2), 2)
})

