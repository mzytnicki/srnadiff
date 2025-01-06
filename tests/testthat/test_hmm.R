context("Checking HMM strategy")

object    <- srnadiffExample()
parameters(object) <- srnadiffDefaultParameters
counts    <- buildDataHmm(object)
pvalues   <- computePvalues(object, counts, 1)$padj
intervals <- hmm(object, counts, pvalues)

test_that("Building data", {
    expect_equal(dim(counts), c(447, 6))
})

test_that("Computing p-values", {
    expect_equal(length(pvalues), 447)
})

test_that("Running HMM", {
    expect_equal(length(intervals), 17)
})
