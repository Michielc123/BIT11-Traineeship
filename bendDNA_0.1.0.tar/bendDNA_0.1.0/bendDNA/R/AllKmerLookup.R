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
#' @return A two-column data.table: sequence or index of k-mer, and its average
#'    bendability coefficient.
#' @details Parameter \strong{k} is defined as number of consecutive trinucleotides rather
#' than number of nucleotides (window size) because a bendability coefficient is
#' defined for a trinucleotide. Size of window is actually k+2 nucleotides.
#'
#' Due to memory restrictions, setting \strong{sequence.out} to FALSE is recommended for
#' bigger values of k. Sequences corresponding to permutation indexes can easily
#' be retrieved using function \code{\link[arrangements]{permutations}} from
#' package \code{\link[arrangements]{arrangements-package}}.
#' @export
#'
#' @examples
#' AllKmerLookup("con", 1, sequence.out=TRUE)
#' AllKmerLookup("dnase", 3, sequence.out=TRUE)
#' AllKmerLookup("con", 8, sequence.out=FALSE)
#'
#' # retrieving the 1st, 100th, 1000th and 10000th sequence from the last call:
#' v <- c(1, 100, 1000, 10000)
#' apply(arrangements::permutations(c("A","C","G","T"), 8+2, replace=TRUE, index=v), 1, paste, collapse="")
AllKmerLookup <- function(scale, k, sequence.out=F) {
  # calculates average bendability coefficient for all possible (k+2)mers
  bendtable <- TrinucleotideScales(scale)

  # due to periodic nature of all permutations of k+2 nucleotides, sum of bendability coefficients per permutation
  # can be calculated as follows:
  LookupVector <- function(bendtable, k) {
    if(k==1) v <- bendtable$refbend
    else v <- rep(LookupVector(bendtable, k-1), each=4) + bendtable$refbend
    return(v)
  }

  # sequences are saved immediately for k=1. for bigger k's, indexes of all permutations are stored in memory and
  # translated to sequences only if requested. averages calculated if needed.
  if(k==1) {
    lookup <- bendtable[, `:=`(sequence=ref, avg.bend=refbend)][, c(3:4)]
  } else {
    lookup <- data.table(LookupVector(bendtable, k))
    lookup <- cbind(seq(nrow(lookup)), lookup/k)
    setnames(lookup, c("index", "avg.bend"))
  }

  # "translate" indexes to sequences if requested:
  if((k!=1) && (sequence.out==T)) {
    sequences <- apply(arrangements::permutations(c("A", "C", "G", "T"), k+2, replace=T, index=lookup$index), 1, paste, collapse="")
    lookup <- lookup[, index:=sequences]
    setnames(lookup, "index", "sequence")
  }

  return(lookup)
}
