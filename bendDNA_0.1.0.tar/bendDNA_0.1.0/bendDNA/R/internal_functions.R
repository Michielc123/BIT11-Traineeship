
TrinucleotideScales <- function(scale) {
  # choose one from con, conrigid, dnase, dnaserigid, nuc, nucrigid
  bendtable <- scales[, c("V1", scale), with=F]
  setnames(bendtable, c("ref", "refbend"))
  return(bendtable)
}


LookupTable <- function(scale, k) {
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
  # later translated to sequences if a match is found. averages calculated if needed.
  if(k==1) {
    lookup <- bendtable[, `:=`(sequence=ref, Lbend=refbend, pref=rep(1:16, each=4), suff=rep(1:16, times=4), bend=refbend)][, c(3:7)]
  } else {
    lookup <- data.table(LookupVector(bendtable, k))
    lookup <- cbind(seq(nrow(lookup)), lookup/k, rep(1:16, each=(4^(k+2))/16), rep(1:16, times=(4^(k+2))/16))
    setnames(lookup, c("index", "Lbend", "pref", "suff"))
    lookup <- lookup[, bend:=Lbend]
  }

  # lookup table contains sequence/permutation index, suffix and prefix for that sequence, and two copies of bendability
  # coefficients (the extra one is needed for data.table merging later - hope to fix that issue).
  return(lookup)
}


QueryProfile <- function(query, bendtable, k, tolerance) {
  # depends on finished bendability table rather than scale to avoid potential conflicts in function MatchBendability

  # split query into trinucleotides and match corresponding bendability coefficients from bendtable:
  trinucleotides <- data.table(stringr::str_sub(query, seq(1, stringr::str_length(query)-2, 1), seq(3, stringr::str_length(query), 1)))
  bends <- bendtable[trinucleotides, on=.(ref=V1), nomatch=0]$refbend

  # split query into kmers (subsequences of length k+2, with 2 nt overlap):
  querykmers <- data.table(stringr::str_sub(query, seq(1, stringr::str_length(query)-k-1, k), seq(k+2, stringr::str_length(query), k)))

  # calculate average bendability coefficient for each kmer (i.e. k trinucleotides):
  if(k==1) means <- bends
  else means <- caTools::runmean(bends, k, endrule = "NA", alg="fast", align="right")[seq(k, stringr::str_length(query)-2, k)]

  # preserve original order of subsequences in query; add lower and upper bounds of interval for matching:
  queryprofile <- querykmers[, ':=' (pos=seq(nrow(querykmers)), Qbend=means, LB=means-tolerance, UB=means+tolerance)]
  return(queryprofile)
}


Bendability <- function(query, k, tolerance, bendtable, lookup) {
  # scale is already determined by bendtable and lookup (k as well, so be careful it remains the same)
  queryprofile <- QueryProfile(query, bendtable, k, tolerance)

  # query kmers are matched with tolerance (matched on interval) with all possible elements from lookup table:
  setNumericRounding(2)
  x <- lookup[queryprofile, on=.(bend>=LB, bend<=UB), allow.cartesian=T, nomatch=0
              ][, deltabend:=abs(Qbend-Lbend)
                ][, c(8,1,10,3,4)]

  # for k>1, actual sequences are added instead of indexes:
  if(k!=1) {
    sequences <- apply(arrangements::permutations(c("A", "C", "G", "T"), k+2, replace=T, index=x$index), 1, paste, collapse="")
    x <- x[, index:=sequences]
    setnames(x, "index", "sequence")
  }

  # consecutive kmers which overlap (suffix of one matches the prefix of next) are merged into full sequences:
  setorder(x, pos)
  x <- split(x, by="pos")
  out <- x[[1]]
  for(i in 2:nrow(queryprofile)) {
    out <- out[x[[i]], on=.(suff=pref), allow.cartesian=T, nomatch=0
               ][, ':='(sequence=paste0(sequence, stringr::str_sub(i.sequence,3,k+2)), deltabend=(deltabend+i.deltabend), pref=pref, suff=i.suff)
                 ][, .(sequence, deltabend, pref, suff)]
  }

  # calculate approximate string distance between match and query:
  out <- cbind(out, adist(out$sequence, query))[, c(1,2,5,3,4)]
  setnames(out, "V1", "strdist")
  setorder(out, deltabend, strdist)

  # output contains additional information: deltabend (cummulative difference between bend. coefficients of query
  # and its match) and strdist (generalized Levenshtein distance)
  return(out)
}

