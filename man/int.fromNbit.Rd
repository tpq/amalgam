% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1-nbits.R
\name{int.fromNbit}
\alias{int.fromNbit}
\title{Convert Bits to Integers}
\usage{
int.fromNbit(codon, nbits = 4)
}
\arguments{
\item{codon}{A binary vector.}

\item{nbits}{The number of bits used per integer.}
}
\description{
This function converts bits into integers.
}
\examples{
int.fromNbit(c(0,0,1,1,0,0,1,1), nbits = 4)
}
