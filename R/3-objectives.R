#' Calculate Fitness of Binary Vector
#'
#' This objective function seeks to maximize the correlation
#'  between the Aitchison distance of the complete composition
#'  and the Aitchison distance of the amalgamation.
#'
#' @param codon A binary vector.
#' @param ARGS Handled by \code{\link{prepareArgs}}.
#' @export
objective.keepDist <- function(codon, ARGS){

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
    clr <- compositions::clr(A)
    NEWDIST <- as.matrix(stats::dist(clr))
  }

  # Maximize correlation between AMALG dist and TARGET dist
  DIST <- ARGS$TARGET
  Od <- DIST[lower.tri(DIST)]
  Ad <- NEWDIST[lower.tri(NEWDIST)]
  stats::cor(Od, Ad)
}

#' Calculate Fitness of Binary Vector
#'
#' This objective function seeks to maximize the percent of the
#'  total variance present in the matrix z that can be explained
#'  by the amalgams.
#'
#' @param codon A binary vector.
#' @param ARGS Handled by \code{\link{prepareArgs}}.
#' @export
objective.maxRDA <- function(codon, ARGS){

  W <- do.call(ARGS$weights, list(codon, n.amalgams = ARGS$n.amalgams))
  A <- ARGS$x %*% W

  # Don't allow zeros!
  if(any(A == 0)){
    return(-1000000)
  }

  # Maximize variance in Z explained by AMALG
  tryCatch({

    if(ARGS$asSLR){

      slr <- as.slr(A)
      v <- vegan::rda(ARGS$z, slr)

    }else{

      ilr <- compositions::ilr(A)
      v <- vegan::rda(ARGS$z, ilr)
    }
  }, error = function(e){

    v <- numeric() # handle error when SVD fails to converge
  })

  if(length(v) == 0) return(-1000000)
  varExplained <- sum(v$CCA$eig / v$tot.chi)
  if(length(varExplained) == 0) return(-1000000) # handle when no CCA inertia
  return(varExplained)
}

#' Calculate Fitness of Binary Vector
#'
#' This objective function seeks to maximize the percent of the
#'  total variance present in the amalgamation that can be explained
#'  by the matrix z.
#'
#' @param codon A binary vector.
#' @param ARGS Handled by \code{\link{prepareArgs}}.
#' @export
objective.maxRDA2 <- function(codon, ARGS){

  W <- do.call(ARGS$weights, list(codon, n.amalgams = ARGS$n.amalgams))
  A <- ARGS$x %*% W

  # Don't allow zeros!
  if(any(A == 0)){
    return(-1000000)
  }

  # Maximize variance in Z explained by AMALG
  tryCatch({

    if(ARGS$asSLR){

      slr <- as.slr(A)
      v <- vegan::rda(slr, ARGS$z)

    }else{

      ilr <- compositions::ilr(A)
      v <- vegan::rda(ilr, ARGS$z)
    }
  }, error = function(e){

    v <- numeric() # handle error when SVD fails to converge
  })

  if(length(v) == 0) return(-1000000)
  varExplained <- sum(v$CCA$eig / v$tot.chi)
  if(length(varExplained) == 0) return(-1000000) # handle when no CCA inertia
  return(varExplained)
}

# objective.maxInvDiag <- function(codon, ARGS){
#
# }
