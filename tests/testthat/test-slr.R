library(amalgam)

set.seed(1)
A <- randAcomp(6)

test_that("as.slr works as expected", {

  expect_equal(
    as.slr(A)[,1],
    log(A[,1] / A[,2])
  )
})

test_that("asSLR argument passed correctly", {

  expect_error(
    amalgam(A, n.amalgams = 3, asSLR = TRUE, maxiter = 5)
  )

  res <- amalgam(A, n.amalgams = 6, asSLR = FALSE, maxiter = 5)
  expect_equal(
    res$SLR,
    NULL
  )

  res <- amalgam(A, n.amalgams = 6, asSLR = TRUE, maxiter = 5)
  expect_equal(
    ncol(res$SLR),
    3
  )
})
