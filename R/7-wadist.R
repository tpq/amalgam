#' Calculate Weighted Aitchison Distance
#'
#' This function computes the weighted aitchison distance proposed by
#'  Greenacre & Lewi (2009) in \code{doi:10.1007/s00357-009-9027-y}.
#'  Unlike the unweighted Aitchison distance, the weighted Aitchison
#'  distance has distributional equivalence.
#'
#' @param A A matrix of compositional data.
#' @return A matrix of distances.
#' @export
wadist <- function(A){

  # close the column weights (per Greenacre & Lewi "for notational convenience")
  col_weight <- colSums(A)
  col_weight <- col_weight / sum(col_weight)

  # logM is the (log of) the weighted geometric mean
  # -- defined as sum(col_weight * log(x)) / sum(col_weight)
  # -- but sum(col_weight) = 1 so we can ignore it
  logM <- apply(A, 1, function(row) sum(col_weight * log(row)))

  # like a CLR, we subtract the log geometric mean from the log data
  centered <- sweep(log(A), 1, logM, "-")

  # now, we use the centered data to calculate a weighted distance
  newdf <- sweep(centered, 2, col_weight, function(x, y) x * sqrt(y))
  as.matrix(stats::dist(newdf))
}

#' Calculate Fitness of Binary Vector
#'
#' This objective function seeks to maximize the correlation
#'  between the weighted Aitchison distance of the complete composition
#'  and the weighted Aitchison distance of the amalgamation.
#'
#' @param codon A binary vector.
#' @param ARGS Handled by \code{\link{prepareArgs}}.
#' @export
objective.keepWADIST <- function(codon, ARGS){

  W <- do.call(ARGS$weights, list(codon, n.amalgams = ARGS$n.amalgams))
  A <- ARGS$x %*% W

  # Don't allow zeros!
  if(any(A == 0)){
    return(-1000000)
  }

  # Calculate distance for SLR or amalgams
  if(ARGS$asSLR){
    slr <- as.slr(A)
    NEWDIST <- as.matrix(stats::dist(slr))
  }else{
    NEWDIST <- wadist(A)
  }

  # Maximize correlation between AMALG dist and TARGET dist
  DIST <- ARGS$TARGET
  Od <- DIST[lower.tri(DIST)]
  Ad <- NEWDIST[lower.tri(NEWDIST)]
  stats::cor(Od, Ad)
}
