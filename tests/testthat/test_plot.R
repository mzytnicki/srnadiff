library(srnadiff)
library(testthat)

context("Checking plot")

exp <- srnadiffExample()
exp <- srnadiff(exp)
plot <- plotRegions(exp, regions(exp, 0.05)[1])

test_that("Running plot", {
    expect_true(is.list(plot))
    expect_equal(length(plot), 4)
    expect_is(plot[[1]], "AnnotationTrack")
    expect_is(plot[[2]], "GenomeAxisTrack")
    expect_is(plot[[3]], "AnnotationTrack")
    expect_is(plot[[4]], "DataTrack")
})
