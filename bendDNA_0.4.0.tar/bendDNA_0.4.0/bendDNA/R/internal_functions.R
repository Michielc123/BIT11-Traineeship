
Bendability <- function(query, scale, k, tolerance, lookup) {

  # get bendability profile of query, together with tolerances:
  bendprofile <- BendabilityProfile(query, scale, k, plot.it=F)
  querykmers <- data.table(stringr::str_sub(query, seq(1, stringr::str_length(query)-k-1, k), seq(k+2, stringr::str_length(query), k)))
  queryprofile <- cbind(querykmers, data.table(pos = names(bendprofile), transpose(bendprofile))[!1, ])
  names(queryprofile) <- c("sequence", "pos", "Qbend")
  queryprofile <- queryprofile[, ':=' (pos=as.integer(pos), Qbend=as.numeric(Qbend))][, ':=' (LB=Qbend-tolerance, UB=Qbend+tolerance)]

  # query kmers are matched with tolerance (matched on interval) with all possible elements from lookup table:
  setNumericRounding(2)
  x <- lookup[queryprofile, on=list(bend>=LB, bend<=UB), allow.cartesian=T, nomatch=0
              ][, deltabend:=abs(Qbend-Lbend)
                ][, c(8,1,10,3,4)]

  # if k-mers in lookup are still represented by indexes, those are translated to sequences:
  if(names(x)[2] == "index") {
    sequences <- apply(arrangements::permutations(c("A", "C", "G", "T"), k+2, replace=T, index=x$index), 1, paste, collapse="")
    x <- x[, index:=sequences]
    setnames(x, "index", "sequence")
  }

  setorder(x, pos)
  iter <- floor(log2(max(x$pos)))

  # consecutive kmers which overlap (suffix of one matches the prefix of next) are merged into full sequences.
  # it's done iteratively by merging pairs of neighbouring sequences, until the whole thing is reconstructed:
  for(i in 1:iter) {

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
utils::globalVariables(c("bend", "LB", "UB", "Qbend", "index", "pos"))
