
TrinucleotideScales <- function(scale) {
  # choose one from con, conrigid, dnase, dnaserigid, nuc, nucrigid
  bendtable <- scales[, c("V1", scale), with=F]
  setnames(bendtable, c("ref", "refbend"))
  return(bendtable)
}
utils::globalVariables("scales")


QueryProfile <- function(query, bendtable, k, tolerance) {
  # depends on finished bendability table rather than scale to avoid potential conflicts in function MatchBendability

  # split query into trinucleotides and match corresponding bendability coefficients from bendtable:
  trinucleotides <- data.table(stringr::str_sub(query, seq(1, stringr::str_length(query)-2, 1), seq(3, stringr::str_length(query), 1)))
  bends <- bendtable[trinucleotides, on=list(ref=V1), nomatch=0]$refbend

  # split query into kmers (subsequences of length k+2, with 2 nt overlap):
  querykmers <- data.table(stringr::str_sub(query, seq(1, stringr::str_length(query)-k-1, k), seq(k+2, stringr::str_length(query), k)))

  # calculate average bendability coefficient for each kmer (i.e. k trinucleotides):
  if(k==1) means <- bends
  else means <- frollmean(bends, k, fill=NA, algo="fast", align="right")[seq(k, stringr::str_length(query)-2, k)]

  # preserve original order of subsequences in query; add lower and upper bounds of interval for matching:
  queryprofile <- querykmers[, ':=' (pos=seq(nrow(querykmers)), Qbend=means, LB=means-tolerance, UB=means+tolerance)]
  return(queryprofile)
}
utils::globalVariables("V1")


Bendability <- function(query, k, tolerance, bendtable, lookup) {
  # scale is already determined by bendtable and lookup (k as well, so be careful it remains the same)
  queryprofile <- QueryProfile(query, bendtable, k, tolerance)

  # query kmers are matched with tolerance (matched on interval) with all possible elements from lookup table:
  setNumericRounding(2)
  x <- lookup[queryprofile, on=list(bend>=LB, bend<=UB), allow.cartesian=T, nomatch=0
              ][, deltabend:=abs(Qbend-Lbend)
                ][, c(8,1,10,3,4)]

  # if k-mers are still represented by indexes, those are translated to sequences:
  if(names(x)[2] == "index") {
    sequences <- apply(arrangements::permutations(c("A", "C", "G", "T"), k+2, replace=T, index=x$index), 1, paste, collapse="")
    x <- x[, index:=sequences]
    setnames(x, "index", "sequence")
  }

  setorder(x, pos)

  # consecutive kmers which overlap (suffix of one matches the prefix of next) are merged into full sequences.
  # it's done iteratively by merging pairs of neighbouring sequences, until the whole thing is reconstructed:
  while(max(x$pos>1)) {

    # first, split matches according to odd or even position:
    xo <- x[pos%%2==1]
    xe <- x[pos%%2==0][, pos:=(pos-1)]

    # (what to do if there's an odd number of groups):
    if(max(xo$pos)>max(xe$pos)) {
      lasteven <- xe[pos==max(pos), ][xo[pos==max(pos), ], on=list(suff=pref), allow.cartesian=T, nomatch=0
                                      ][, ':='(sequence=paste0(sequence, stringr::str_sub(i.sequence,start=3)), deltabend=(deltabend+i.deltabend), pref=pref, suff=i.suff)
                                        ][, list(pos, sequence, deltabend, pref, suff)]
      xo <- xo[!(pos==max(pos)), ]
      xe <- rbind(xe[!(pos==max(pos)), ], lasteven)
    }

    # now merge pairs (first group with second, third with fourth etc.):
    x <- xo[xe, on=list(pos=pos, suff=pref), allow.cartesian=T, nomatch=0
            ][, ':='(sequence=paste0(sequence, stringr::str_sub(i.sequence,start=3)), deltabend=(deltabend+i.deltabend), pref=pref, suff=i.suff)
              ][, pos:=.GRP, by="pos"
                ][, list(pos, sequence, deltabend, pref, suff)]
  }

  out <- x[, list(sequence, deltabend, pref, suff)]

  # calculate approximate string distance between match and query:
  out <- cbind(out, utils::adist(out$sequence, query))
  names(out)[5] <- "strdist"
  out <- out[, c(1,2,5,3,4)]
  setorder(out, deltabend, strdist)

  # output contains additional information: deltabend (cummulative difference between bend. coefficients of query
  # and its match) and strdist (generalized Levenshtein distance)
  return(out)
}
utils::globalVariables(c("bend", "LB", "UB", "deltabend", "Qbend", "Lbend", "index", "pos", "pref", "i.sequence", "i.deltabend", "i.suff", "suff", "strdist"))

