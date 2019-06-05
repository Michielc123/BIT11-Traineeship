#' Plot bendability of sequence(s)
#'
#' @param sequences A DNAStringSet object or a character vector of same-length DNA
#'    sequences.
#' @param scale One from "con", "conrigid", "dnase", "dnaserigid", "nuc", "nucrigid".
#' @param k Number of consecutive trinucleotides for which to calculate and plot the 
#'    average bendability coefficient.
#' @param sample.seq If it is NULL (default), all sequences from input will be plotted. 
#'    If a single integer is provided, that many input sequences will be chosen at 
#'    random. Can also be a vector of indexes of desired sequences (their row positions 
#'    in the input data.table).
#'
#' @return A plot.
#' @export
#'
#' @examples
PlotBendability <- function(sequences, scale, k, sample.seq=NULL, r) {
  # takes full-sequence output of functions Bendability and ApplyBendability and plots (average) bend.coefficients
  # for all or chosen sequences 
  
  bendtable <- BendabilityTable(scale)
  
  # if sample.seq is not provided, all sequences are plotted;
  # if sample.seq is a single integer, that many sequences are chosen at random;
  # alternatively, sample.seq can be a vector of indexes of sequences from output we wish to plot:
  if(is.null(sample.seq)) { index <- 1:length(sequences)
  } else if(length(sample.seq)==1) { index <- sample(1:length(sequences), sample.seq)
  } else { index <- sample.seq }
  sq <- as.character(sequences[index])
  
  # funtion Yaxis returns (average) bendability coefficients given one sequence:
  Yaxis <- function(query, k) {
    trinucleotides <- data.table(str_sub(query, seq(1, width(query)-2, 1), seq(3, width(query), 1)))
    bends <- bendtable[trinucleotides, on=.(ref=V1), nomatch=0]$refbend
    if(k==1) { means <- bends
    } else { means <- runMean(bends, k)[seq(k, width(query)-2, k)] }
    return(means)
  }
  
  # apply function Yaxis to all chosen sequences:
  mat <- sapply(sq, Yaxis, k)
  colnames(mat) <- index
  melted <- melt(mat, measure.vars = colnames(mat))
  print(class(melted))
  print(class(melted$Var1))
  print(class(melted$value))
  melted$Var2 <- as.factor(melted$Var2)
  
  ggplot(melted, aes(x = Var1, y = value, group = Var2, col = Var2)) + 
    geom_spline(spar = r) + 
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
    guides(col = guide_legend(ncol = 2))
  
}

