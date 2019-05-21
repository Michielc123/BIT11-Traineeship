Consensus_list_vector <- function(list, cutoff){

  Consensus <- function(mf, cutoff){

    b <- DNAStringSet(mf)
    c <- IdClusters(myDistMatrix = NULL, method = "inexact", myXStringSet = b, cutoff = cutoff)
    b <- as.data.table(b)
    c <- as.data.table(c)
    d <- cbind(b,c)
    d <- split(d, d$cluster)
    e <- as.vector(unlist(sapply(d, function(x){consensusString(DNAStringSet(x[["x"]]))})))
  }

  return(lapply(list, function(x){Consensus(x[["sequence"]], cutoff)}))

}
