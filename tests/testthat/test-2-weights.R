library(amalgam)
bin <- c(0,0,0,0,0,0,0,1,1,1,1,1)

test_that("Nto1 weights function works as expected", {

  Nto1 <- weight.Nto1(bin, n.amalgams = 3) # nbits = 2
  expect_equal(
    colSums(Nto1),
    c(1,0,2)
  )

  Nto1 <- weight.Nto1(bin, n.amalgams = 7) # nbits = 3
  expect_equal(
    colSums(Nto1),
    c(0,0,1,0,0,0,1)
  )
})

test_that("NtoN weights function works as expected", {

  NtoN <- weight.NtoN(bin, n.amalgams = 1)
  expect_equal(
    nrow(NtoN),
    12
  )
  expect_equal(
    rowSums(NtoN)[1],
    0
  )

  NtoN <- weight.NtoN(bin, n.amalgams = 2)
  expect_equal(
    nrow(NtoN),
    6
  )
  expect_equal(
    rowSums(NtoN)[1],
    0
  )

  NtoN <- weight.NtoN(bin, n.amalgams = 3)
  expect_equal(
    nrow(NtoN),
    4
  )
  expect_equal(
    rowSums(NtoN)[1],
    0
  )

  NtoN <- weight.NtoN(bin, n.amalgams = 4)
  expect_equal(
    nrow(NtoN),
    3
  )
  expect_equal(
    rowSums(NtoN)[1],
    0
  )
})
