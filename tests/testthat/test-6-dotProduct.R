library(amalgam)
library(compositions)
`%*%` <- compositions::`%*%` # needed for testthat?

set.seed(1)
x <- randAcomp(20)[1,]
w <- randAcomp(20)[1,]

X <- randAcomp(4, 5)
W <- t(randAcomp(3, 5))

test_that("adot and amat agree with compositions scalar product", {

  expect_equal(
    adot(x, w),
    acomp(x) %*% acomp(w)
  )

  expect_equal(
    as.vector(amat(x, w)),
    acomp(x) %*% acomp(w)
  )

  expect_equal(
    amat(X, W)[3,2],
    acomp(X[3,]) %*% acomp(W[,2])
  )
})
