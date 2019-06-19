#' Find sequences with similar bendability profiles
#'
#' Calculates bendability profile of query sequence and finds all other sequences
#' which have similar (or same) profile.
#'
#' @param query DNAString or a character string.
#' @param scale One from "con", "conrigid", "dnase", "dnaserigid", "nuc", "nucrigid".
#' @param k Number of consecutive trinucleotides for which to calculate the average
#'    bendability coefficient.
#' @param tolerance All k-mers with average bendability coefficient which falls into
#'    interval <query k-mer bendability +/- tolerance> are matched with query.
#' @param wsize Window size for calculations. Should be (n*k)+2. See Details.
#' @param output.list Applicable if wsize is set. Whether to output a list of
#'    data.tables which contain matches for one window, or a single data.table
#'    with full-length sequences. Defaults to FALSE.
#' @param random.out Whether to output all results (FALSE, default) or a random
#'    selection (TRUE).
#' @param size.cutoff Applicable when random.out=TRUE. If a result has more than
#'    size.cutoff matches, those matches are grouped by prefix-suffix combination and
#'    a fraction of each group is chosen at random. Defaults to 20.
#' @param frac Applicable when random.out=TRUE (see above). For a group with x members,
#'    size of fraction kept is defined as max(1, round(x*frac)). Defaults to 0.2.
#' @param lookup Optional. Output of function AllKmerLookup. If not supplied, the
#'    function makes lookup table automatically for each run.
#'
#' @return Single data.table, or a list of data.tables:
#' \itemize{
#'  \item{sequence} : {sequence which matches bendability profile of query}
#'  \item{deltabend} : {cummulative difference between bendability coefficients of match and query}
#'  \item{strdist} : {generalized Levenshtein distance between match and query}
#' }
#' @details Higher tolerance values are not recommended in combination with bigger k.
#'
#' If parameter \strong{wsize} is set, query sequence will be split to chunks for processing.
#' That can speed up the execution and enable exporting results as a list, which is
#' useful in case a run with default parameters failed or took too long to finish.
#' Chunks are of length wsize, with a two-nucleotide overlap beetween consecutive
#' ones. Chosen value must satisfy the conditions \emph{wsize=(n*k)+2} and
#' \emph{length(query)>=2*wsize}. The downside is that last <wsize nucleotides of
#' query sequence will not be processed.
#'
#' Setting \strong{output.list} to FALSE increases the runtime and memory consumption
#' dramatically, especially for longer sequences and higher tolerance values
#' respectively . If the execution failed with default parameters, try setting
#' output.list to TRUE, or processing your sequence in chunks.
#'
#' Setting \strong{random.out} to TRUE will increase speed. Recommended if a random subset
#' of matches is enough for your application or a run with random.out=F already failed.
#'
#' Lookup table depends only on parameters scale and k. If you plan to apply the
#' function on multiple sequences without changing those two parameters, providing
#' \strong{lookup} speeds up the process because the table doesn't have to be
#' created again for each run.
#' @export
#'
#' @examples
#' MatchBendability("TGATTCCTAAAGTCA", "con", k=1, tolerance=0.5)
#' MatchBendability("TGATTCCTAAAGTCA", "con", k=3, tolerance=0.1)
MatchBendability <- function(query, scale, k, tolerance, wsize=NULL, output.list=F, random.out=F, size.cutoff=20, frac=0.2, lookup=NULL) {

  # get bendability scale table:
  bendtable <- TrinucleotideScales(scale)

  # get lookup table, if it isn't already provided:
  if(is.null(lookup)) lookup <- LookupTable(scale, k, sequence.out=F)

  # additional column with bendabilities (a sad artefact needed for data.table merging, which I will try to fix):
  lookup <- lookup[, bend:=Lbend]

  # if a window has more than size.cutoff matches, group them by prefix-suffix combination and keep a fraction
  # of each group at random (size of fraction is determined by parameter frac). recommended if a random subset
  # of matches is enough for your application, and/or a very large number of matches is expected (happens with
  # long sequences, big k, big tolerance):
  PickRandomly <- function(dt, size.cutoff, frac) {
    if(nrow(dt)>size.cutoff) dt <- dt[, .SD[sample(.N, max(1, round(.N*frac)))], by=list(pref, suff)][, c(3:5, 1:2)]
    return(dt)
  }

  if(is.null(wsize)){

    # if wsize is not defined, query sequence is processed as a whole:
    out <- Bendability(query, k, tolerance, bendtable, lookup)
    if(random.out==T) out <- PickRandomly(out, size.cutoff, frac)
    out <- out[, c(1,2,3)]

  } else {

    # if wsize is defined, query is split into windows of size wsize. function Bendability is applied to each window:
    querywindows <- as.data.table(stringr::str_sub(query, seq(1, stringr::str_length(query)-wsize+1, wsize-2), seq(wsize, stringr::str_length(query), wsize-2)))
    l <- apply(querywindows, 1, Bendability, k, tolerance, bendtable, lookup)

    # resulting list is filtered to leave only entries which have a predecessor and a successor in the preceding and
    # next windows, respectively:
    filtfw <- sapply(seq(1, length(l)-1), function(x) which(l[[x]]$suff %in% l[[x+1]]$pref), simplify=F)
    filtfw[[length(l)]] <- seq(nrow(l[[length(l)]]))
    filtrev <- sapply(seq(length(l), 2), function(x) which(l[[x]]$pref %in% l[[x-1]]$suff), simplify=F)
    filtrev[[length(l)]] <- seq(nrow(l[[1]]))
    filtrev <- rev(filtrev)
    filtall <- mapply(intersect, filtfw, filtrev)
    filtered <- sapply(seq(length(l)), function(x) l[[x]][c(filtall[[x]]),], simplify = F)

    if(random.out==T) filtered <- lapply(filtered, PickRandomly, size.cutoff, frac)

    # output can be a list of data.tables, each data.table representing all matches for one window. if full sequences
    # are required, the thing has to be iteratively merged, which will increase the runtime:
    if(output.list==T) { out <- lapply(filtered, function(x) x[, c(1:3)])
    } else {

      x <- rbindlist(filtered, idcol=T)
      setnames(x, ".id", "pos")

      # consecutive windows which overlap (suffix of one matches the prefix of next) are merged into full sequences.
      # it's done iteratively by merging pairs of neighbouring sequences, until the whole thing is reconstructed:
      while(max(x$pos>1)) {

        # first, split matches according to odd or even position:
        xo <- x[pos%%2==1]
        xe <- x[pos%%2==0][, pos:=(pos-1)]

        # (what to do if there's an odd number of groups):
        if(max(xo$pos)>max(xe$pos)) {
          lasteven <- xe[pos==max(pos), ][xo[pos==max(pos), ], on=list(suff=pref), allow.cartesian=T, nomatch=0
                                          ][, ':='(sequence=paste0(sequence, stringr::str_sub(i.sequence,start=3)), deltabend=(deltabend+i.deltabend), strdist=(strdist+i.strdist), pref=pref, suff=i.suff)
                                            ][, list(pos, sequence, deltabend, strdist, pref, suff)]
          xo <- xo[!(pos==max(pos)), ]
          xe <- rbind(xe[!(pos==max(pos)), ], lasteven)
        }

        # now merge pairs (first group with second, third with fourth etc.):
        x <- xo[xe, on=list(pos=pos, suff=pref), allow.cartesian=T, nomatch=0
                ][, ':='(sequence=paste0(sequence, stringr::str_sub(i.sequence,start=3)), deltabend=(deltabend+i.deltabend), strdist=(strdist+i.strdist), pref=pref, suff=i.suff)
                  ][, pos:=.GRP, by="pos"
                    ][, list(pos, sequence, deltabend, strdist, pref, suff)]
      }

      out <- x[, list(sequence, deltabend, strdist)]
    }
  }
  return(out)
}

utils::globalVariables(c("pref", "i.sequence", "deltabend", "i.deltabend", "strdist", "i.strdist", "i.suff", "suff"))
