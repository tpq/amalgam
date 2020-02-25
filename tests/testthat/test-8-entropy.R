library(amalgam)
library(entropy)

test_that("pre-KL shrinkage matches KL with shrinkage", {

  set.seed(1)
  y1 = c(4, 2, 3, 1, 10, 4)
  y2 = c(2, 3, 7, 1, 4, 3)
  y1c <- freqs(y1, method = "shrink")
  y2c <- freqs(y2, method = "shrink")

  expect_equal(
    KL.shrink(y1, y2)[1],
    KL.empirical(y1c, y2c)
  )
})

test_that("our KLD matches entropy package", {

  set.seed(1)
  A <- randAcomp(5)

  for(i in 1:nrow(A)){
    expect_equal(
      amalgam:::KLD(A[1,], A[i,]),
      KL.empirical(A[1,], A[i,])
    )
  }
})

test_that("fast SKL matches slow SKL", {

  set.seed(1)
  A <- randAcomp(5)

  expect_equal(
    round(
      SKL(A),
      8
    ),
    round(
      amalgam:::SKL_slow(A),
      8
    )
  )
})

test_that("SKL divergence has distributional equivalence", {

  # A is a matrix where columns 6 and 7 are identical
  A <- randAcomp(5)
  A <- cbind(A, cbind(rep(.2, 5), rep(.2, 5)))
  A <- sweep(A, 1, rowSums(A), "/")

  # B is a matrix where columns 6 and 7 are merged
  B <- A
  B[,6] <- B[,6] + B[,7]
  B <- B[,-7]

  expect_equal(
    SKL(A),
    SKL(B)
  )
})
