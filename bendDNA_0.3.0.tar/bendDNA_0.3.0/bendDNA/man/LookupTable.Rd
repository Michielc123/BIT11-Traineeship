% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LookupTable.R
\name{LookupTable}
\alias{LookupTable}
\title{Average bendability coefficients of all possible k-mers}
\usage{
LookupTable(scale, k, sequence.out = F)
}
\arguments{
\item{scale}{One from "con", "conrigid", "dnase", "dnaserigid", "nuc", "nucrigid".}

\item{k}{Number of consecutive trinucleotides for which to calculate the average
bendability coefficient.}

\item{sequence.out}{Whether to output k-mer sequences (TRUE) or indexes of
alphabetically sorted permutations (FALSE). Defaults to FALSE.}
}
\value{
A four-column data.table:
\itemize{
 \item{sequence} : {sequence of length k+2}
 \item{index} : {alternatively, index of k-mer permutation}
 \item{Lbend} : {average bendability of k-mer}
 \item{pref} : {first two nucleotides of sequence/permutation}
 \item{suff} : {last two nucleotides of sequence/permutation}
}
}
\description{
Produces all possible k-mers and calculates their average bendability coefficients.
Upper limit for k is determined by the size of available RAM.
}
\details{
Parameter \strong{k} is defined as number of consecutive trinucleotides rather
than number of nucleotides (window size) because a bendability coefficient is
defined for a trinucleotide. Size of window is actually k+2 nucleotides.

Due to memory restrictions, setting \strong{sequence.out} to FALSE is recommended for
bigger values of k. Sequences corresponding to permutation indexes can easily
be retrieved using function \code{\link[arrangements]{permutations}} from
package \code{\link[arrangements]{arrangements-package}}.

Dinucleotide prefixes and suffixes are represented with numbers 1-16 (corresponding
to their position when ordered alphabetically) to facilitate later matching and
reduce required memory.
}
\examples{
LookupTable("con", 1, sequence.out=TRUE)
LookupTable("dnase", 3, sequence.out=TRUE)
LookupTable("con", 8, sequence.out=FALSE)

# retrieving the 1st, 100th, 1000th and 10000th sequence from the last call:
v <- c(1, 100, 1000, 10000)
b <- c("A","C","G","T")
apply(arrangements::permutations(b, 8+2, replace=TRUE, index=v), 1, paste, collapse="")
}
