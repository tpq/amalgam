#' Convert Bits to Weights Matrix
#'
#' This function builds a weights matrix from a binary vector.
#'  The resultant weights matrix with have \code{n.amalgams} columns
#'  and a number of rows based on the size of the binary vector.
#'  The resultant weights matrix will have at most 1 value per row.
#'  It is possible for row sums or column sums to equal 0.
#'
#' The interpretation of this weights matrix is that each component
#'  will only contribute to one amalgam node.
#'
#' @param codon A binary vector.
#' @param n.amalgams The number of columns in the weights matrix.
#' @examples
#' # For 3 amalgams, read as "node 1", "node 3", "node 2", "node 0"
#' weight.Nto1(c(0,1,1,1,1,0,0,0), 3)
#' @export
weight.Nto1 <- function(codon, n.amalgams){

  int <- int.fromNbit(codon, nbits = ceiling(log2(n.amalgams+1)))

  # Surplus nbits may produce false amalgams
  int[int > n.amalgams] <- 0

  nrow <- length(int)
  W <- matrix(0, nrow, n.amalgams)
  for(i in 1:nrow(W)){
    amalgTo <- int[i]
    W[i,amalgTo] <- 1
  }
  W
}

#' Convert Bits to Weights Matrix
#'
#' This function builds a weights matrix from a binary vector.
#'  The resultant weights matrix with have \code{n.amalgams} columns
#'  and a number of rows based on the size of the binary vector.
#'  The resultant weights matrix will have rows that sum to 1.
#'  It is possible for row sums or column sums to equal 0.
#'
#' The interpretation of this weights matrix is that each component
#'  can contribute to multiple amalgam nodes, but must contribute
#'  an equal amount to each node.
#'
#' @param codon A binary vector.
#' @param n.amalgams The number of columns in the weights matrix.
#' @examples
#' # For 2 amalgams, read as "node 2", "node 1+2", "node 1", "node 0"
#' weight.NtoN(c(0,1,1,1,1,0,0,0), 2)
#' @export
weight.NtoN <- function(codon, n.amalgams){

  W <- t(matrix(codon, ncol = length(codon)/n.amalgams, nrow = n.amalgams))
  W <- sweep(W, 1, FUN = "/", rowSums(W))
  W[is.na(W)] <- 0
  W
}
