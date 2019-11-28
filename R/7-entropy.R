#' Calculate Shannon's Index
#' @param p A composition.
shannon <- function(p){

  non0 <- p[p>0]
  -1* sum(non0*log(non0))
}

#' Calculate a Normalized Shannon's Index
#' @param p A composition.
uniformity <- function(p){

  H <- shannon(p)
  log(length(p)) - H
}

#' Calculate Fitness of Binary Vector
#'
#' This objective function seeks to preserve the relative entropy
#'  of the amalgamation for all samples.
#'
#' @param codon A binary vector.
#' @param ARGS Handled by \code{\link{prepareArgs}}.
#' @export
objective.keepEntropy <- function(codon, ARGS){

  W <- do.call(ARGS$weights, list(codon, n.amalgams = ARGS$n.amalgams))
  A <- ARGS$x %*% W

  # Don't allow zeros!
  if(any(A == 0)){
    return(-1000000)
  }

  if(ARGS$asSLR){
    stop("This objective does not support summed log-ratios.")
  }

  ref <- apply(ARGS$x, 1, uniformity)
  new <- apply(A, 1, uniformity)
  -1 * sum(abs(ref - new))
}

#' Calculate Fitness of Binary Vector
#'
#' This objective function seeks to maximize the (average)
#'  between-group difference in the entropy of the amalgamation.
#'  This function expects exactly two groups.
#'
#' @param codon A binary vector.
#' @param ARGS Handled by \code{\link{prepareArgs}}.
#' @export
objective.diffEntropy <- function(codon, ARGS){

  W <- do.call(ARGS$weights, list(codon, n.amalgams = ARGS$n.amalgams))
  A <- ARGS$x %*% W

  # Don't allow zeros!
  if(any(A == 0)){
    return(-1000000)
  }

  if(ARGS$asSLR){
    stop("This objective does not support summed log-ratios.")
  }

  new <- apply(A, 1, shannon)
  groups <- ARGS$z[,1]
  if(length(unique(groups)) != 2) stop("This objective requires exactly two groups.")
  grp1 <- groups == unique(groups)[1]
  grp2 <- groups == unique(groups)[2]
  mean(new[grp1]) - mean(new[grp2])
}
