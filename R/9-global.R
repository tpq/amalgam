#' Create Random Poisson Data
#'
#' This function makes a random matrix by sampling rows
#'  from a Poisson distribution.
#'
#' @param nrows An integer. The number of rows in the matrix.
#' @param ncols An integer. The number of columns in the matrix.
#' @param lambda Argument passed to \code{rpois}.
#' @return A matrix.
#' @export
randPois <- function (nrows, ncols = nrows, lambda = 100){
  counts <- stats::rpois(nrows * ncols, lambda = lambda)
  matrix(counts, nrows, ncols)
}

#' Create Random Acomp Data
#'
#' This function makes a random matrix by sampling rows
#'  from a Poisson distribution, then closing them.
#'
#' @inheritParams randPois
#' @return A matrix.
#' @export
randAcomp <- function(nrows, ncols = nrows, lambda = 100){
  A <- compositions::acomp(randPois(nrows, ncols, lambda))
  class(A) <- "matrix"
  A
}

#' Package Check
#'
#' Checks whether the user has the required package installed.
#'  For back-end use only.
#'
#' @param package A character string. An R package.
packageCheck <- function(package){

  if(!requireNamespace(package, quietly = TRUE)){
    stop("Uh oh! This amalgam method depends on ", package, ".")
  }
}
