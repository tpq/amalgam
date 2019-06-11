#' Convert Bits to Integers
#'
#' This function converts bits into integers.
#'
#' @param codon A binary vector.
#' @param nbits The number of bits used per integer.
#' @examples
#' int.fromNbit(c(0,0,1,1,0,0,1,1), nbits = 4)
#' @export
int.fromNbit <- function(codon, nbits = 4){

  bitMatrix <- sapply(1:nbits, function(index) codon[1:nbits == index])
  if(length(codon) == nbits){
    bitString <- paste0(bitMatrix, collapse = "")
  }else{
    bitString <- apply(bitMatrix, 1, paste0, collapse = "")
  }
  strtoi(bitString, base = 2)
}

#' Convert Bits to Integers
#'
#' This function converts bits into integers.
#'
#' Use \code{\link{int.fromNbit}} instead.
#'
#' @param codon A binary vector.
#' @examples
#' int.from2bit(c(0,0,1,1,0,0,1,1))
#' @export
int.from2bit <- function(codon){

  bits <- paste0(codon[c(TRUE,FALSE)], codon[c(FALSE,TRUE)])
  strtoi(bits, base = 2)
}
