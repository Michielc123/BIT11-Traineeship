#' Bendability profile of DNA sequence
#'
#' Splits sequence(s) into (k+2)-mers with two-nucleotide overlap and calculates average
#' bendability per (k+2)-mer according to the chosen scale. Plots the result.
#'
#' @param sequences A DNAStringSet object or a character vector.
#' @param scale One from "con", "conrigid", "dnase", "dnaserigid", "nuc", "nucrigid".
#' @param k Number of consecutive trinucleotides for which to calculate the average
#'    bendability coefficient.
#' @param plot.it Whether to plot the bendability profile(s) or not. Defaults to TRUE.
#' @param smt Degree of line smoothing. Between 0 and 1, defaults to 0 (no smoothing).
#'
#' @return A data.table, or a plot and a data.table. Each input sequence is represented
#' by one row in data.table, which contains average bendability coefficients for (k+2)-mers.
#' @details Parameter \strong{k} is defined as number of consecutive trinucleotides rather
#' than number of nucleotides (or window size) because a bendability coefficient is
#' defined for a trinucleotide. Size of window is actually k+2 nucleotides.
#' @export
#'
#' @examples
#' BendabilityProfile("TGATTCCTAAAGTCA", "con", 1, plot.it=TRUE)
#'
#' s <- c("TGATTCCTAAAGTCA", "CCTGAAATGCTAGCGT", "AAACTAGCCTCGATG")
#' BendabilityProfile(s, "con", 1, plot.it=TRUE, smt=0.3)
#' BendabilityProfile(s, "con", 3, plot.it=FALSE)
BendabilityProfile <- function(sequences, scale, k, plot.it=TRUE, smt=0) {

  bendtable <- scales[, c("V1", scale), with=F]
  setnames(bendtable, c("ref", "refbend"))
  sq <- as.character(sequences)

  # funtion Yaxis returns (average) bendability coefficients given one sequence:
  Yaxis <- function(query, k) {
    trinucleotides <- data.table(stringr::str_sub(query, seq(1, stringr::str_length(query)-2, 1), seq(3, stringr::str_length(query), 1)))
    bends <- bendtable[trinucleotides, on=list(ref=V1), nomatch=0]$refbend
    if(k==1) { means <- bends
    } else { means <- frollmean(bends, k, fill=NA, algo="fast", align="right")[seq(k, stringr::str_length(query)-2, k)] }
    return(means)
  }

  # apply function Yaxis to all sequences:
  l <- lapply(sq, Yaxis, k)
  dt <- rbindlist(lapply(l, function(x) data.table(t(x))), fill = TRUE)
  names(dt) <- as.character(1:ncol(dt))

  if(!is.null(names(sequences))) { out <- cbind(as.data.table(names(sequences)), dt)
  } else { out <- cbind(as.data.table(rownames(dt)), dt) }

  if(plot.it==TRUE) {

    # adjust format for plotting:
    melted <- melt(t(dt))
    names(melted) <- c("position", "sequence", "value")
    melted <- melted[!(is.na(melted$value)), ]
    melted$sequence <- as.factor(melted$sequence)

    # plot it:
    plot <- ggplot(melted, aes(x = position, y = value, group = sequence, col = sequence)) +
      ggformula::geom_spline(spar = smt) +
      labs(x = "\nWindow", y = "Bendability\n", colour = NULL) +
      ggtitle(paste0("Average bendability on window size ", k+2, ", scale ", scale)) +
      theme(plot.title = element_text(hjust = 0.5, size = 35),
            legend.text = element_text(size = 20),
            legend.key.size = unit(0.8,"cm"),
            plot.margin = unit(c(1,1,1,1), "cm"),
            legend.margin = margin(0,0,0,100),
            axis.text.x = element_text(size = 30, margin = margin(c(20,0,0,0))),
            axis.text.y = element_text(size = 30, margin = margin(c(0,20,0,0))),
            axis.title = element_text(size = 30)) +
      guides(col = guide_legend(ncol = 1))

    return(list(out, plot))

  } else { return(out) }

}
utils::globalVariables(c("V1", "value", "position"))
