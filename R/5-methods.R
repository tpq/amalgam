#' Show \code{amalgam} Object
#'
#' This function shows an \code{amalgam} object.
#'
#' @param x An \code{amalgam} object.
#' @param ... Not used.
#' @export
print.amalgam <- function(x, ...){

  cat("AMALGAM OBJECT\n")
  cat("\n")

  cat("Dimension of original data:\n")
  print(dim(x$original))
  cat("\n")

  cat("Dimension of amalgamated data:\n")
  print(dim(x$amalgams))
  cat("\n")

  cat("Size of amalgams:\n")
  print(colSums(x$weights))
  cat("\n")

  cat("Size of contributions:\n")
  print(table(rowSums(x$weights)))
  cat("\n")

  cat("Fitness:\n")
  print(x$fitness)
  cat("\n")

  cat("Use plot() method!")
}

#' Plot \code{amalgam} Object
#'
#' This function plots an \code{amalgam} object.
#'
#' @param x An \code{amalgam} object.
#' @param col A character vector. Used to color the figure.
#' @param center A logical. Toggles whether to center the ternary plot.
#' @param scale A logical. Toggles whether to scale the ternary plot.
#' @param a1,a2,a3 The amalgams to show in the ternary plot.
#' @param ... Not used.
#' @importFrom compositions var.acomp
#' @export
plot.amalgam <- function(x, col = rep(1, nrow(x$amalgams)),
                         center = FALSE, scale = FALSE,
                         a1 = 1, a2 = 2, a3 = 3,
                         ...){

  cols <- viridis::viridis(length(unique(col)))[factor(col)]

  graphics::par(mfrow=c(2,2))

  pc <- stats::prcomp(compositions::ilr(x$original))$x
  vars <- apply(pc,2,stats::var)
  vars <- round(vars / sum(vars) * 100, 1)
  graphics::plot(pc[,1], pc[,2], col = cols,
                 xlab = paste0("PC1 of full ilr (",
                               vars[1], "%)"),
                 ylab = paste0("PC2 of full ilr (",
                               vars[2], "%)"),
                 main = "ilr-PCA of full data", pch = 16)

  graphics::plot(x$ga, legend = FALSE)

  pc <- stats::prcomp(compositions::ilr(x$amalgams))$x
  vars <- apply(pc,2,stats::var)
  vars <- round(vars / sum(vars) * 100, 1)
  graphics::plot(pc[,1], pc[,2], col = cols,
                 xlab = paste0("PC1 of amalgam ilr (",
                               vars[1], "%)"),
                 ylab = paste0("PC2 of amalgam ilr (",
                               vars[2], "%)"),
                 main = "ilr-PCA of amalgams", pch = 16)

  # compositions::acomp uses "robust" option for center/scale
  #  which gives error without library(compositions)
  if(is.null(getOption("robust"))){
    options("robust" = FALSE)
  }

  graphics::plot(compositions::acomp(x$amalgams[,c(a1,a2,a3)]),
                 col = cols, center = center, scale = scale,
                 main = "3-part amalgam", pch = 16)

  graphics::par(mfrow=c(1,1))
}
