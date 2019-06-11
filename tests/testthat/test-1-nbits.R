library(amalgam)
bin <- c(0,0,0,0,0,0,0,1,1,1,1,1)

test_that("nbits function deparses binary input correctly", {

  expect_equal(
    int.from2bit(bin),
    c(0,0,0,1,3,3)
  )

  expect_equal(
    int.from2bit(bin),
    int.fromNbit(bin, nbits = 2)
  )

  expect_equal(
    int.fromNbit(bin, nbits = 4),
    c(0,1,15)
  )
})
