#' Average bendability coefficients of all possible k-mers
#'
#' Produces all possible k-mers and calculates their average bendability coefficients.
#' Upper limit for k is determined by the size of available RAM.
#'
#' @param scale One from "con", "conrigid", "dnase", "dnaserigid", "nuc", "nucrigid".
#' @param k Number of consecutive trinucleotides for which to calculate the average
#'    bendability coefficient.
#' @param sequence.out Whether to output k-mer sequences (TRUE) or indexes of
#'    alphabetically sorted permutations (FALSE). Defaults to FALSE.
#'
#' @return A four-column data.table:
#' \itemize{
#'  \item{sequence} : {sequence of length k+2}
#'  \item{index} : {alternatively, index of k-mer permutation}
#'  \item{Lbend} : {average bendability of k-mer}
#'  \item{pref} : {first two nucleotides of sequence/permutation}
#'  \item{suff} : {last two nucleotides of sequence/permutation}
#' }
#' @details Parameter \strong{k} is defined as number of consecutive trinucleotides rather
#' than number of nucleotides (window size) because a bendability coefficient is
#' defined for a trinucleotide. Size of window is actually k+2 nucleotides.
#'
#' Due to memory restrictions, setting \strong{sequence.out} to FALSE is recommended for
#' bigger values of k. Sequences corresponding to permutation indexes can easily
#' be retrieved using function \code{\link[arrangements]{permutations}} from
#' package \code{\link[arrangements]{arrangements-package}}.
#'
#' Dinucleotide prefixes and suffixes are represented with numbers 1-16 (corresponding
#' to their position when ordered alphabetically) to facilitate later matching and
#' reduce required memory.
#' @export
#'
#' @examples
#' LookupTable("con", 1, sequence.out=TRUE)
#' LookupTable("dnase", 3, sequence.out=TRUE)
#' LookupTable("con", 8, sequence.out=FALSE)
#'
#' # retrieving the 1st, 100th, 1000th and 10000th sequence from the last call:
#' v <- c(1, 100, 1000, 10000)
#' b <- c("A","C","G","T")
#' apply(arrangements::permutations(b, 8+2, replace=TRUE, index=v), 1, paste, collapse="")
LookupTable <- function(scale, k, sequence.out=F) {
  # calculates average bendability coefficient for all possible (k+2)mers
  bendtable <- TrinucleotideScales(scale)

  # due to periodic nature of all permutations of k+2 nucleotides, sum of bendability coefficients per permutation
  # can be calculated as follows:
  LookupVector <- function(bendtable, k) {
    if(k==1) v <- bendtable$refbend
    else v <- rep(LookupVector(bendtable, k-1), each=4) + bendtable$refbend
    return(v)
  }

  # sequences are saved immediately if k=1. for bigger ks, indexes of all permutations are stored in memory and
  # translated to sequences only if requested. averages calculated if needed.
  if(k==1) {
    lookup <- bendtable[, `:=`(sequence=ref, Lbend=refbend, pref=rep(1:16, each=4), suff=rep(1:16, times=4))][, c(3:6)]
  } else {
    lookup <- data.table(LookupVector(bendtable, k))
    lookup <- cbind(seq(nrow(lookup)), lookup/k, rep(1:16, each=(4^(k+2))/16), rep(1:16, times=(4^(k+2))/16))
    setnames(lookup, c("index", "Lbend", "pref", "suff"))
  }

  # "translate" indexes to sequences if requested:
  if((k!=1) && (sequence.out==T)) {
    sequences <- apply(arrangements::permutations(c("A", "C", "G", "T"), k+2, replace=T, index=lookup$index), 1, paste, collapse="")
    lookup <- lookup[, index:=sequences]
    setnames(lookup, "index", "sequence")
  }

  # lookup table contains sequence or permutation index, bendability coefficient, and suffix and prefix
  return(lookup)
}
utils::globalVariables(c("ref", "refbend", "Lbend"))
