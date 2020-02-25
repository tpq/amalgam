#' Calculate Shannon's Index
#' @param p A composition.
shannon <- function(p){

  # p <- p[p>0]
  -1 * sum(p*log(p))
}

#' Calculate a Normalized Shannon's Index
#' @param p A composition.
uniformity <- function(p){

  H <- shannon(p)
  log(length(p)) - H
}

#' Calculate Kullback-Leibler Divergence
#' @param p,q A composition.
KLD <- function(p, q){

  # no0 <- p>0 & q>0
  # p <- p[no0]
  # q <- q[no0]
  sum(p * log(p / q))
}

#' Calculate Symmetric Kullback-Leibler divergence
#'
#' This function computes the Symmetric Kullback-Leibler divergence.
#'
#' @param A A matrix of compositional data.
#' @return A matrix of distances.
SKL_slow <- function(A){

  pairwise <- function(x, y){

    Dxy <- KLD(x, y)
    Dyx <- KLD(y, x)
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

#' Calculate Symmetric Kullback-Leibler divergence
#'
#' This function computes the Symmetric Kullback-Leibler divergence.
#'
#' @param A A matrix of compositional data.
#' @return A matrix of distances.
#' @export
SKL <- function(A){

  # sum(p * log(p/q))
  #   = sum(p * log(p) - p * log(q))
  #   = sum(p * log(p)) - sum(p * log(q))
  term2 <- A %*% t(log(A))
  term1 <- diag(term2) %*% t(rep(1, nrow(A)))
  KLD <- term1 - term2
  SKL <- (KLD + t(KLD))/2
  return(SKL)
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
objective.maxDEntropy <- function(codon, ARGS){

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
  A <- sweep(A, 1, rowSums(A), "/")

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

#' Calculate Fitness of Binary Vector
#'
#' This function is adapted from \code{\link{objective.maxRDA}}.
#'  It uses the PCoA of the SKL divergences to maximize the percent of the
#'  total variance in the PCoA that can be explained by the matrix z.
#'  In formula notation, RDA(PCoA ~ z).
#'
#' @param codon A binary vector.
#' @param ARGS Handled by \code{\link{prepareArgs}}.
#' @export
objective.maxDSKL <- function(codon, ARGS){

  W <- do.call(ARGS$weights, list(codon, n.amalgams = ARGS$n.amalgams))
  A <- ARGS$x %*% W
  A <- sweep(A, 1, rowSums(A), "/")

  # Don't allow zeros!
  if(any(A == 0)){
    return(-1000000)
  }

  if(ARGS$asSLR){
    stop("This objective function does not support summed log-ratios.")
  }

  # Maximize variance in Z explained by AMALG
  tryCatch({

    divergences <- SKL(A)
    df <- min(nrow(A)-1, ncol(A)-1)
    PCoA <- stats::cmdscale(divergences, k = df)
    v <- vegan::rda(PCoA, ARGS$z)

  }, error = function(e){

    v <- numeric() # handle error when SVD fails to converge
  })

  if(length(v) == 0) return(-1000000)
  varExplained <- sum(v$CCA$eig / v$tot.chi)
  if(length(varExplained) == 0) return(-1000000) # handle when no CCA inertia
  return(varExplained)
}
