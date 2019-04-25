library(srnadiff)
library(testthat)

context("Checking plot")

exp <- srnadiffExample()
exp <- srnadiff(exp)
plot <- plotRegions(exp, regions(exp, 0.05)[1])

test_that("Running plot", {
    expect_true(is.list(plot))
    expect_equal(length(plot), 3)
    expect_is(plot[[1]], "GenomeAxisTrack")
    expect_is(plot[[2]], "AnnotationTrack")
    expect_is(plot[[3]], "DataTrack")
})
