#' The \code{amalgam} Package
#'
#' @description
#' Welcome to the \code{amalgam} package!
#'
#' The \code{amalgam} function finds a useful data-driven amalgamation
#'  for the provided compositional data set (where rows are samples and
#'  columns are components). It works like this:
#'
#' 1. The genetic algorithm suggests random binary vector solutions.
#'
#' 2. Each binary vector is turned into a weights matrix that defines a
#'  rule for amalgamating the composition. See \code{\link{weight.Nto1}}
#'  and \code{\link{weight.NtoN}} for the weights matrix types.
#'
#' 3. An objective function determines the goodness of the weights matrix.
#'  See \code{\link{objective.keepDist}} and \code{\link{objective.maxRDA}}
#'  for the objective function types.
#'
#' 4. The genetic algorithm "breeds" and "mutates" the best solutions.
#'  Over time, the optimal solution emerges.
#'
#' @author Thom Quinn
#'
#' @param x A matrix. The input data. Rows are samples and columns are components.
#' @param n.amalgams An intger. How many components the amalgamation should have.
#' @param maxiter An integer. How long the genetic algorithm should run.
#' @param z A matrix. The constraining matrix. Optional.
#' @param objective A function. The objective function. See above.
#' @param weights A function. The weights function. See above.
#' @param asSLR A boolean. Toggles whether to turn the amalgams into a set
#'  of summed log-ratios (SLRs). See \code{\link{as.slr}}.
#' @param ... Arguments passed to \code{GA::ga} function.
#' @return An \code{amalgam} S3 object.
#' @examples
#' simData <- randAcomp(5, 10)
#' result <- amalgam(simData, n.amalgams = 3, objective = objective.keepDist)
#' print(result)
#' plot(result)
#' @export
amalgam <- function(x, n.amalgams = 3, maxiter = ncol(x)*10, z = NULL,
                    objective = objective.keepDist,
                    weights = weight.Nto1,
                    asSLR = FALSE, ...){

  ARGS <- prepareArgs(x = x, n.amalgams = n.amalgams, maxiter = maxiter, z = z,
                      objective = objective, weights = weights, asSLR = asSLR, ...)

  # Pass ... to GA::ga without breaking objective function
  call <- list(type = "binary", fitness = ARGS$objective,
               nBits = ARGS$totalBits, maxiter = ARGS$maxiter,
               ARGS = ARGS) # args given to objective function
  call <- c(call, ARGS$forGA) # args given to ga function

  # Use the provided objective function
  res <- do.call(GA::ga, call)

  # Calculate fitness of the best solution
  solution <- res@solution[1,]
  FITNESS <- do.call(ARGS$objective, list(codon = solution, ARGS = ARGS))

  # Get the best amalgams
  W <- do.call(ARGS$weights, list(codon = solution, n.amalgams = ARGS$n.amalgams))
  A <- ARGS$x %*% W
  rownames(W) <- colnames(ARGS$x)
  rownames(A) <- rownames(ARGS$x)

  # Get SLRs (optional)
  if(ARGS$asSLR){
    slr <- as.slr(A)
  }else{
    slr <- NULL
  }

  # Prepare S3 results object
  out <- list(
    "ga" = res,
    "totalBits" = ARGS$totalBits,
    "solution" = solution,
    "fitness" = FITNESS,
    "constraints" = ARGS$z,
    "amalgams" = A,
    "original" = ARGS$x,
    "original.no0" = ARGS$x.no0,
    "weights" = W,
    "SLR" = slr)
  class(out) <- "amalgam"
  return(out)
}

#' Argument Handler Function
#'
#' This function handles arguments for \code{\link{amalgam}}.
#'
#' @inheritParams amalgam
#' @return A list of arguments.
#' @export
prepareArgs <- function(x, n.amalgams = 3, maxiter = ncol(x)*10, z = NULL,
                        objective = objective.keepDist,
                        weights = weight.Nto1,
                        asSLR = FALSE, ...){

  # Coerce x as.matrix (needed for data.frame and acomp input)
  rowWeights <- rowSums(x)
  x <- as.matrix(x)
  class(x) <- "matrix"
  x <- sweep(x, 1, rowWeights, "/")
  message("Alert: Compositional data closed to sum to 1.")

  # Coerce z as.data.frame
  z <- as.data.frame(z)

  # Collect arguments as a list
  ARGS <- list(x = x, n.amalgams = n.amalgams, maxiter = maxiter, z = z,
               objective = objective, weights = weights, asSLR = asSLR, ...)
  ARGS$forGA <- as.list(substitute(list(...)))[-1]

  # Replace zeros if needed...
  if(any(x == 0)){
    message("Alert: Replacing zeros with zCompositions for TARGET calculation (if applicable).")
    packageCheck("zCompositions")
    ARGS$x.no0 <- zCompositions::cmultRepl(x, method = "CZM")
  }else{
    ARGS$x.no0 <- x
  }

  # Calculate "TARGET" from complete data (if applicable)
  if(identical(objective, objective.keepDist)){

    # Find Aitchison distance for non-zero data
    ARGS$TARGET <- as.matrix(stats::dist(compositions::clr(ARGS$x.no0)))
    message("Alert: Aitchison distance TARGET calculation complete.")
  }

  # Calculate "TARGET" from complete data (if applicable)
  if(identical(objective, objective.keepWADIST)){

    # Find Aitchison distance for non-zero data
    ARGS$TARGET <- wadist(ARGS$x.no0)
    message("Alert: Weighted Aitchison distance TARGET calculation complete.")
  }

  # Calculate "TARGET" from complete data (if applicable)
  if(identical(objective, objective.keepSKL)){

    # if(all(rowWeights <= 1)){
    #   stop("The SKL objective requires count-based input.")
    # }
    #
    # # NOTE: SKL requires count-based inputs!
    # ARGS$x <- sweep(ARGS$x, 1, rowWeights, "*")
    # ARGS$x.no0 <- sweep(ARGS$x.no0, 1, rowWeights, "*")
    # message("Alert: Compositional data re-scaled to counts (needed for the SKL objective).")

    # Find SKL divergence for non-zero data
    ARGS$TARGET <- SKL(ARGS$x.no0)
    message("Alert: Kullback-Leibler divergence TARGET calculation complete.")
  }

  if(identical(objective, objective.maxRDA) |
     identical(objective, objective.maxRDA2) |
     identical(objective, objective.diffEntropy)){

    if(nrow(ARGS$z) == 0){
      stop("Please provide a valid constraining matrix to argument 'z'.")
    }
  }

  # Calculate totalBits needed based on parameters
  if(identical(ARGS$weights, weight.Nto1)){

    ARGS$totalBits <- ncol(ARGS$x) * ceiling(log2(ARGS$n.amalgams+1))
    message("Alert: ", ARGS$totalBits, " bits needed for genetic algorithm.")

  }else if(identical(ARGS$weights, weight.NtoN)){

    ARGS$totalBits <- ncol(ARGS$x) * ARGS$n.amalgams
    message("Alert: ", ARGS$totalBits, " bits needed for genetic algorithm.")

  }else{

    stop("Provided 'weights' method not supported.")
  }

  return(ARGS)
}

#' Turn Components into Log-Ratios
#'
#' The \code{as.slr} function turns a composition
#'  (or \code{amalgam} object) into a set of log-ratios by taking
#'  the log of the first column over the second, the log of the
#'  third column over the fourth, and so on.
#'
#' When the input is an \code{amalgam} object, the result
#'  is a set of summed log-ratios (SLRs).
#'
#' @param amalgams A \code{matrix} or \code{amalgam} object.
#' @return A \code{matrix}.
#' @export
as.slr <- function(amalgams){

  if(class(amalgams) == "amalgam"){
    return(as.slr(amalgams$amalgams))
  }

  if(ncol(amalgams) %% 2 == 1){
    stop("Must have an even number of amalgams to turn them into summed log-ratios (SLRs).")
  }

  odds <- 1:ncol(amalgams) %% 2 == 1
  evens <- 1:ncol(amalgams) %% 2 == 0

  slr <- log(amalgams[,odds] / amalgams[,evens])

  return(slr)
}
