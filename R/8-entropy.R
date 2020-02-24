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
  A <- sweep(A, 1, rowSums(A), "/")

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
  A <- sweep(A, 1, rowSums(A), "/")

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

#' Calculate Symmetric Kullback-Leibler divergence
#'
#' This function computes the (symmetric) Kullback-Leibler divergence.
#'
#' @param A A matrix of compositional data.
#' @return A matrix of distances.
#' @export
SKL <- function(A){

  pairwise <- function(x, y){

    Dxy <- entropy::KL.empirical(x, y)#, verbose = FALSE)
    Dyx <- entropy::KL.empirical(y, x)#, verbose = FALSE)
    re <- (Dxy + Dyx) / 2
    return(re[1])
  }

  res <- matrix(0, nrow(A), nrow(A))
  for(i in 1:nrow(A)){
    for(j in 1:nrow(A)){
      if(i < j){
        res[i, j] <- pairwise(A[i,], A[j,])
        res[j, i] <- res[i, j]
      }
    }
  }
  return(res)
}

#' Calculate Fitness of Binary Vector
#'
#' This objective function seeks to maximize the correlation
#'  between the (symmetric) Kullback-Leibler divergence of the
#'  complete composition and the (symmetric) Kullback-Leibler
#'  divergence of the amalgamation.
#'
#' @param codon A binary vector.
#' @param ARGS Handled by \code{\link{prepareArgs}}.
#' @export
objective.keepSKL <- function(codon, ARGS){

  W <- do.call(ARGS$weights, list(codon, n.amalgams = ARGS$n.amalgams))
  A <- ARGS$x %*% W

  # Don't allow zeros!
  if(any(A == 0)){
    return(-1000000)
  }

  # Calculate distance for SLR or amalgams
  if(ARGS$asSLR){
    stop("This objective function does not support summed log-ratios.")
  }else{
    NEWDIST <- SKL(A)
  }

  # Maximize correlation between AMALG dist and TARGET dist
  DIST <- ARGS$TARGET
  Od <- DIST[lower.tri(DIST)]
  Ad <- NEWDIST[lower.tri(NEWDIST)]
  stats::cor(Od, Ad)
}
