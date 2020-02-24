# A is a matrix where columns 6 and 7 are identical
library(amalgam)
A <- randAcomp(5)
A <- cbind(A, cbind(rep(.2, 5), rep(.2, 5)))
A <- sweep(A, 1, rowSums(A), "/")

# B is a matrix where columns 6 and 7 are merged
B <- A
B[,6] <- B[,6] + B[,7]
B <- B[,-7]

test_that("SKL divergence has distributional equivalence", {

  expect_equal(
    SKL(A),
    SKL(B)
  )
})
