#' Bendability plotting for a list of sequences.
#'
#' `Plot()` takes a list of sequences and plots the average bendability scores of a certain amout of randomly selected sequences.
#'
#' @param df Dataframe with sequence column.
#' @param k Number of consecutive trinucleotides to calculate the average from.
#' @param y (optional) Number of randomly selected sequences from the dataset.
#'
#' @return A plot of the average bendability Scores.
#'
#' @import ggplot2
#' @import data.table
#' @import Biostrings
#' @import plyr
#' @import stringr
#' @import TTR
#' @import dplyr
#'
#' @export


Plot <- function(df, k, y=NULL){

  tabel <- fread("table_bendability_UPPERCASE.txt", sep = "\t", header = F)
  tabel <- as.data.table(tabel)
  colnames(tabel) <- c("Trinucl", "Score")
  tabel$Score <- as.numeric(gsub(",",".",tabel$Score))
  setkey(tabel, Trinucl)

  df <- as.data.frame(as.character(df$sequence))
  colnames(df) <- "sequence"

  if(!is.null(y)){
      df$id <- rownames(df)
      df <- sample_n(df, y)
  }

  Windowmaker <- function(l,k){

      st2 <- seq(1, nchar(as.character(l)) - 2)
      en2 <- seq(3, nchar(as.character(l)))

      TriNucl <- str_sub(l, st2, en2)
      TriNucl <- as.data.table(TriNucl, sorted = F)
      colnames(TriNucl) <- "Trinucl"
      a <- merge(TriNucl,tabel, by = "Trinucl", sort = F)
      if(k > 1){runMean(a$Score, n = k)[seq(k,length(a$Score),k)]}
      else{return(a$Score)}

  }

  b <- ldply(lapply(df$sequence, Windowmaker, k = k), rbind)
  rownames(b) <- df$id
  b <- as.data.table(t(b))
  colnames(b) <- gsub("V","",colnames(b))

  mlt <- melt(b, measure.vars = colnames(b))
  mlt$id <- rep(1:nrow(b), nrow(mlt)/nrow(b))

  ggplot(mlt, aes(x = as.factor(id), y = value, group = variable, col = variable)) +
    geom_line() +
    labs(x = "\nWindow", y = "Bendability\n", colour = NULL) +
    ggtitle(paste0("Average bendability with window size: ", k, "\n")) +
    theme(plot.title = element_text(hjust = 0.5, size = 35),
          legend.text = element_text(size = 17),
          legend.key.size = unit(0.8,"cm"),
          plot.margin = unit(c(1,3,3,3), "cm"),
          legend.margin = margin(0,0,0,100),
          axis.text.x = element_text(size = 30, margin = margin(c(20,0,0,0))),
          axis.text.y = element_text(size = 30, margin = margin(c(0,20,0,0))),
          axis.title = element_text(size = 30)) +
    guides(col = guide_legend(ncol = 2))

  ggsave("Plot.pdf", width = unit(27, "cm"), height = unit(19,"cm"))
}
