context("Checking main functions")

exp <- srnadiffExample()

test_that("Testing regions method", {
    exp <- srnadiff(exp)
    expect_equal(length(regions(exp)), 26)
})

test_that("Running with different strategies", {
    exp <- srnadiff(exp, segMethod=c("annotation", "hmm"))
    expect_equal(length(regions(exp)), 178)
})

test_that("Running with different sizes", {
    exp <- srnadiff(exp, useParameters=list(minSize=15, maxSize=30))
    expect_equal(length(regions(exp)), 26)
})

test_that("Running with different minimum depth", {
    exp <- srnadiff(exp, useParameters=list(minDepth=1))
    expect_equal(length(regions(exp)), 26)
})

test_that("Running with different transition probabilities", {
    exp <- srnadiff(exp, segMethod="hmm",
                     useParameters=list(noDiffToDiff=0.5,diffToNoDiff=0.5))
    expect_equal(length(ranges), 1)
})

test_that("Running with different emission probabilities", {
    exp <- srnadiff(exp, segMethod="hmm", useParameters=list(emission=0.75))
    expect_equal(length(regions(exp)), 17)
})

test_that("Running with different emission threshold", {
    exp <- srnadiff(exp, segMethod="hmm",
                     useParameters=list(emissionThreshold=0.5))
    expect_equal(length(regions(exp)), 17)
})

test_that("Running with different number of overlapping base pairs", {
    exp <- srnadiff(exp, segMethod="hmm", useParameters=list(minOverlap=15))
    expect_equal(length(regions(exp)), 17)
})

test_that("Running several threads", {
    exp <- srnadiff(exp)
    exp2 <- srnadiffExample()
    exp2 <- srnadiff(exp2, nThreads=2)
    expect_equal(regions(exp), regions(exp2))
})

test_that("Running main function", {
    exp <- srnadiff(exp)
    expect_equal(length(regions(exp)), 26)
})

test_that("Running redundant regions removal", {
    regions  <- GRanges(seqnames = "14", strand = "+",
                        ranges = IRanges(start = c(60000000, 60000100),
                        width = 10))
    regions2 <- removeRedundant(regions)
    expect_equal(length(regions2), 2)
})

