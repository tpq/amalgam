library(amalgam)
library(compositions)
library(GA)

set.seed(1)
A <- randAcomp(6)
ARGS <- prepareArgs(A, z = c("A", "A", "A", "B", "B", "B"))

res <- amalgam(A, maxiter = 5, z = c("A", "A", "A", "B", "B", "B"),
               objective = objective.keepDist)
test_that("objective.keepDist makes sense", {

  expect_equal(
    objective.keepDist(res$solution, ARGS),
    res$fitness
  )

  AMALG <- as.matrix(dist(clr(res$amalgams)))
  expect_equal(
    cor(ARGS$TARGET[lower.tri(ARGS$TARGET)],
        AMALG[lower.tri(AMALG)]),
    res$fitness
  )

  expect_equal(
    res$original %*% res$weights,
    res$amalgams
  )

  expect_equal(
    length(res$solution),
    res$totalBits
  )
})

res <- amalgam(A, maxiter = 5, z = c("A", "A", "A", "B", "B", "B"),
               objective = objective.maxRDA)
test_that("objective.maxRDA makes sense", {

  expect_equal(
    objective.maxRDA(res$solution, ARGS),
    res$fitness
  )

  v <- vegan::rda(ilr(res$amalgams), data.frame(ARGS$z))
  expect_equal(
    sum(v$CCA$eig / v$tot.chi),
    res$fitness
  )
})
