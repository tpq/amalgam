#' Aitchison Scalar Product for Vectors
#'
#' This function computes the Aitchison scalar product of
#'  two compositions, defined as the sum of the product
#'  of all pairwise log-ratios of \code{x} and of \code{w},
#'  scaled by the number of components.
#'
#' @param x A composition of node values.
#' @param w A composition of node weights.
#' @return A scalar product.
#'
#' @examples
#' x <- randAcomp(20)[1,]
#' w <- randAcomp(20)[1,]
#' adot(x, w)
#' @export
adot <- function(x, w){

  if(!identical(length(x), length(w))){
    stop("adot fails")
  }

  sum <- 0
  for(i in 1:length(x)){
    for(j in 1:length(w)){
      sum <- sum + log(x[i]/x[j])*log(w[i]/w[j])
    }
  }
  mean <- sum/(2*length(x))

  return(mean)
}

#' Aitchison Scalar Product for Matrices
#'
#' This function computes the Aitchison scalar product of
#'  two matrices, The first matrix should describe N sample
#'  compositions (as rows) with D components (as columns),
#'  while the second matrix should describe D weight
#'  compositions (as rows) for M outputs (as columns).
#'
#' @param X A matrix of sample compositions.
#' @param W A matrix of weight compositions.
#' @return A matrix of scalar products.
#'
#' @examples
#' X <- randAcomp(4, 5)
#' W <- t(randAcomp(3, 5))
#' amat(X, W)
#' @export
amat <- function(X, W){

  # If vector, make X a matrix with 1 row
  if(is.null(dim(X))){
    X <- t(X)
  }

  # If vector, make W a matrix with 1 column
  if(is.null(dim(W))){
    W <- as.matrix(W)
  }

  if(!identical(ncol(X), nrow(W))){
    stop("non-conformable arguments")
  }

  mat <- matrix(0, nrow(X), ncol(W))
  for(i in 1:nrow(X)){
    for(j in 1:ncol(W)){
      mat[i,j] <- adot(X[i,], W[,j])
    }
  }

  return(mat)
}
