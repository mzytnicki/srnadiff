library(srnadiff)
library(testthat)

context("Checking main functions")

exp <- srnadiffExample()
exp <- srnadiff(exp)

test_that("Testing regions method", {
    expect_equal(length(regions(exp)), 26)
})

test_that("Running with different strategies", {
    exp2 <- srnadiffExample()
    exp2 <- srnadiff(exp2, segMethod=c("annotation", "hmm"))
    expect_equal(length(regions(exp2)), 176)
})

test_that("Running with different sizes", {
    exp2 <- srnadiffExample()
    exp2 <- srnadiff(exp2, useParameters=list(minSize=15, maxSize=30))
    expect_equal(length(regions(exp2)), 43)
})

test_that("Running with different minimum depth", {
    exp2 <- srnadiffExample()
    exp2 <- srnadiff(exp2, useParameters=list(minDepth=1))
    expect_equal(length(regions(exp2)), 703)
})

test_that("Running with different transition probabilities", {
    exp2 <- srnadiffExample()
    exp2 <- srnadiff(exp2, segMethod="hmm",
                     useParameters=list(noDiffToDiff=0.5,diffToNoDiff=0.5))
    expect_equal(length(ranges), 1)
})

test_that("Running with different emission probabilities", {
    exp2 <- srnadiffExample()
    exp2 <- srnadiff(exp2, segMethod="hmm", useParameters=list(emission=0.75))
    expect_equal(length(regions(exp2)), 9)
})

test_that("Running with different emission threshold", {
    exp2 <- srnadiffExample()
    exp2 <- srnadiff(exp2, segMethod="hmm",
                     useParameters=list(emissionThreshold=0.5))
    expect_equal(length(regions(exp2)), 28)
})

test_that("Running with different number of overlapping base pairs", {
    exp2 <- srnadiffExample()
    exp2 <- srnadiff(exp2, segMethod="hmm", useParameters=list(minOverlap=15))
    expect_equal(length(regions(exp2)), 16)
})

test_that("Running several threads", {
    exp2 <- srnadiffExample()
    exp2 <- srnadiff(exp2, nThreads=2)
    expect_equal(regions(exp), regions(exp2))
})

test_that("Running main function", {
    expect_equal(length(regions(exp)), 26)
})

test_that("Running redundant regions removal", {
    regions  <- GRanges(seqnames = "14", strand = "+",
                        ranges = IRanges(start = c(60000000, 60000100),
                        width = 10))
    regions2 <- removeRedundant(regions)
    expect_equal(length(regions2), 2)
})

