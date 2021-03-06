% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7-wadist.R
\name{objective.keepWADIST}
\alias{objective.keepWADIST}
\title{Calculate Fitness of Binary Vector}
\usage{
objective.keepWADIST(codon, ARGS)
}
\arguments{
\item{codon}{A binary vector.}

\item{ARGS}{Handled by \code{\link{prepareArgs}}.}
}
\description{
This objective function seeks to maximize the correlation
 between the weighted Aitchison distance of the complete composition
 and the weighted Aitchison distance of the amalgamation.
}
