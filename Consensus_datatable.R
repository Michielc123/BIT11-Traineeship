Consensus_list_datatable <- function(list, cutoff){

  Consensus <- function(mf, cutoff){

    b <- DNAStringSet(mf)
    c <- IdClusters(myDistMatrix = NULL, method = "inexact", myXStringSet = b, cutoff = cutoff)
    b <- as.data.table(b)
    c <- as.data.table(c)
    d <- cbind(b,c)
    d <- split(d, d$cluster)
    sequence <- sapply(d, function(x){consensusString(DNAStringSet(x[["x"]]))})
    sequence <- as.data.table(sequence)
  }

  return(lapply(list, function(x){Consensus(x[["sequence"]], cutoff)}))

}
