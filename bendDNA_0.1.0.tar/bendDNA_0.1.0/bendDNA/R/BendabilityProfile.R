#' Bendability profile of DNA sequence
#'
#' Splits sequence into k-mers with two-nucleotide overlap and calculates average
#' bendability per k-mer according to the chosen scale.
#'
#' @param query DNAString or a character string.
#' @param scale One from "con", "conrigid", "dnase", "dnaserigid", "nuc", "nucrigid".
#' @param k Number of consecutive trinucleotides for which to calculate the average
#'    bendability coefficient.
#'
#' @return A three-column data.table:
#' \itemize{
#'  \item{kmer} : {sequence of length k+2}
#'  \item{position} : {relative position of k-mer in query sequence}
#'  \item{avg.bend} : {average bendability of k-mer}
#' }
#' @details Parameter \strong{k} is defined as number of consecutive trinucleotides rather
#' than number of nucleotides (window size) because a bendability coefficient is
#' defined for a trinucleotide. Size of window is actually k+2 nucleotides.
#' @export
#'
#' @examples
#' BendabilityProfile("TGATTCCTAAAGTCA", "con", 1)
#' BendabilityProfile("TGATTCCTAAAGTCA", "con", 3)
#' BendabilityProfile("TGATTCCTAAAGTCA", "nuc", 3)
BendabilityProfile <- function(query, scale, k) {

  bendtable <- TrinucleotideScales(scale)

  # split sequence into trinucleotides and match corresponding bendability coefficients from bendtable:
  trinucleotides <- data.table(stringr::str_sub(query, seq(1, stringr::str_length(query)-2, 1), seq(3, stringr::str_length(query), 1)))
  bends <- bendtable[trinucleotides, on=.(ref=V1), nomatch=0]$refbend

  # split sequence into kmers (subsequences of length k+2, with 2 nt overlap):
  querykmers <- data.table(stringr::str_sub(query, seq(1, stringr::str_length(query)-k-1, k), seq(k+2, stringr::str_length(query), k)))

  # calculate average bendability coefficient for each k-mer (i.e. k trinucleotides):
  if(k==1) means <- bends
  else means <- caTools::runmean(bends, k, endrule = "NA", alg="fast", align="right")[seq(k, stringr::str_length(query)-2, k)]

  # preserve original order of subsequences in query and output the table:
  queryprofile <- querykmers[, ':=' (pos=seq(nrow(querykmers)), avgbend=means)]
  setnames(queryprofile, c("kmer", "position", "avg.bend"))
  return(queryprofile)
}
