#' Clustering and consensus of similar sequences
#'
#' Intended to present output of function MatchBendability in a more read-friendly
#' format. The output can contain very large numbers of highly similar sequences,
#' so clustering and calculating a consensus can help provide a clearer picture.
#' DISCLAIMER: due to the nature of bendability problem and limitations of IUPAC
#' extended genetic alphabet, output of this function is not intended for further machine
#' use. Not all sequences which can be inferred from consensus output will have satisfied
#' the initial bendability requirements (for example AA/TT will be marked as WW, same as
#' AT/TA).
#'
#' @param sequences A DNAStringSet object or a character vector of same-length DNA
#'    sequences.
#' @param cutoff Numeric, between 0 and 1. Distance required for separation of clusters.
#' @param sample.seq If it is NULL (default), all sequences from input will be used in
#'    clustering and consensus calculation. If a single integer x is provided, x
#'    input sequences will be chosen at random. Can also be a vector of indexes of
#'    desired sequences (their row positions in the input data.table).
#'
#' @return Variable number of consensus sequences written in extended genetic alphabet
#'    (IUPAC_CODE_MAP).
#' @details Similar sequences are clustered together using inexact method
#' (sequence-only method, doesn't use a distance matrix). Number of clusters depends
#' on sequence similarity and cutoff value (lower values will produce more clusters).
#' @export
#'
#' @examples
#' dt <- MatchBendability("TGATTCCTAAAGTCA", "con", k=3, tolerance=0.1)
#' ClusterAndConsensus(dt$sequence, cutoff=0.01)
#' ClusterAndConsensus(dt$sequence, cutoff=0.5)
ClusterAndConsensus <- function(sequences, cutoff, sample.seq=NULL){
  # Intended to present Bendability and ApplyBendability output in a more human-readable format. Outputs of those
  # functions usually contain very large numbers of highly similar sequences, so clustering and making consensus
  # sequences can help provide a clearer picture for a human user.
  # Cutoff value should be in range 0-1. Lower values will produce more clusters.
  # DISCLAIMER: due to the nature of bendability problem and limitations of IUPAC extended genetic alphabet, output
  # of this function is not intended for further machine use. Not all sequences which can be inferred from consensus
  # output will have satisfied the initial bendability requirements (for example AA/TT will be marked as WW, same as
  # AT/TA).

  sq <- Biostrings::DNAStringSet(sequences)

  # if sample.seq is not provided, all sequences are used for clustering;
  # if sample.seq is a single integer, that many sequences are chosen at random (recommended for bigger outputs);
  # alternatively, sample.seq can be a vector of indexes of sequences from output we wish to take into account:
  if(is.null(sample.seq)) { sq <- sq
  } else if(length(sample.seq)==1) { sq <- sq[sample(length(sq), sample.seq)]
  } else { sq <- sq[sample.seq, ] }

  # similar sequences are clustered together using inexact method (sequence-only method, doesn't use a distance
  # matrix). number of clusters depends on sequence similarity and cutoff value.
  cluster <- DECIPHER::IdClusters(myDistMatrix = NULL, method = "inexact", myXStringSet = sq, cutoff = cutoff)
  t <- data.table(cbind(sq, cluster))
  t <- split(t, t$clust)

  # a consensus sequence is calculated for each cluster:
  consensus <- sapply(t, function(x){Biostrings::consensusString(Biostrings::DNAStringSet(x[["sq"]]))})
  consensus <- data.table(cbind(consensus), unique(cluster))
  return(consensus)
}
