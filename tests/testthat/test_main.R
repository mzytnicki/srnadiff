library(srnadiff)
library(testthat)

context("Checking main functions")

exp <- sRNADiffExample()
exp <- runAll(exp)

test_that("Testing regions method", {
    expect_equal(length(regions(exp)), 43)
})

test_that("Running with different strategies", {
    exp2 <- sRNADiffExample()
    exp2 <- setStrategies(exp2, TRUE, FALSE, TRUE, TRUE)
    exp2 <- runAll(exp2)
    expect_equal(length(exp2@regions), 35)
})

test_that("Running with different sizes", {
    exp2 <- sRNADiffExample()
    exp2 <- setSizes(exp2, 15, 30)
    exp2 <- runAll(exp2)
    expect_equal(length(exp2@regions), 76)
})

test_that("Running with different minimum depth", {
    exp2 <- sRNADiffExample()
    exp2 <- setMinDepth(exp2, 1)
    exp2 <- runAll(exp2)
    expect_equal(length(exp2@regions), 38)
})

test_that("Running with different merge distance", {
    exp2   <- sRNADiffExample()
    exp2   <- setMergeDistance(exp2, 1)
    ranges <- runAllNaive(exp2)
    expect_equal(length(ranges), 33)
})

test_that("Running with different minimum differences", {
    exp2   <- sRNADiffExample()
    exp2   <- setMinDifferences(exp2, 100)
    ranges <- runAllSlice(exp2)
    expect_equal(length(ranges), 35)
})

test_that("Running with different transition probabilities", {
    exp2   <- sRNADiffExample()
    exp2   <- setTransitionProbabilities(exp2, 0.5, 0.5)
    ranges <- runAllHmm(exp2)
    expect_equal(length(ranges), 4)
})

test_that("Running with different emission probabilities", {
    exp2   <- sRNADiffExample()
    exp2   <- setEmissionProbabilities(exp2, 0.75)
    ranges <- runAllHmm(exp2)
    expect_equal(length(ranges), 3)
})

test_that("Running with different emission threshold", {
    exp2   <- sRNADiffExample()
    exp2   <- setEmissionThreshold(exp2, 0.5)
    ranges <- runAllHmm(exp2)
    expect_equal(length(ranges), 13)
})

test_that("Running with different number of overlapping base pairs", {
    exp2 <- sRNADiffExample()
    exp2 <- setMinOverlap(exp2, 20)
    exp2 <- runAll(exp2)
    expect_equal(length(exp2@regions), 20)
})

test_that("Running several threads", {
    exp2 <- sRNADiffExample()
    exp2 <- setNThreads(exp2, 4)
    exp2 <- runAll(exp2)
    expect_equal(exp@regions, exp2@regions)
})

test_that("Running main function", {
    expect_equal(length(exp@regions), 43)
})

