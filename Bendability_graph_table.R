BendScoresGraph <- function(mf, wdw) {

  library("ggplot2")
  library("data.table")
  library("Biostrings")
  library("plyr")
  library("stringr")

  ### Reading in the bendability table ###

  setwd("C:/Users/Michiel/Desktop/Sequences")
  tabel <- fread("table_bendability_UPPERCASE.txt", sep = "\t", header = F)
  tabel <- as.data.table(tabel)
  colnames(tabel) <- c("Trinucl", "Score")
  tabel$Score <- as.numeric(gsub(",",".",tabel$Score))
  setkey(tabel, Trinucl)

  ### Reading in the Multifasta file ###

  mfasta <- readDNAStringSet(mf, "fasta")

  ### Average calculation per window function ###

  AvgCalc <- function(seq){

    dna <- DNAString(seq)
    view <- successiveViews(dna, width=rep(wdw, length(dna)/wdw))

    Windowmaker <- function(x){

      strt <- seq(1, nchar(x) - 2)
      end <- seq(3, nchar(x))

      TriNucl <- Views(x, start = strt, end = end)
      TriNucl <- as.data.table(TriNucl)
      colnames(TriNucl) <- "Trinucl"
      setkey(TriNucl, Trinucl)

      a <- merge(tabel, TriNucl)
      avg <- mean(as.numeric(a$Score))
      return(avg)
    }

    sapply(view, Windowmaker)

  }

  ### Applying the function to Multifasta and creating the bendability table ###

  bend_per_gene <- lapply(mfasta, AvgCalc)

  bend_per_gene <- ldply(bend_per_gene, rbind)
  bend_per_gene <- t(bend_per_gene)
  colnames(bend_per_gene) <- bend_per_gene[1,]
  bend_per_gene <- bend_per_gene[-1,]
  rownames(bend_per_gene) <- c()

  ### Plotting the bendabilities ###

  mlt <- melt(bend_per_gene, measure.vars = colnames(bend_per_gene))

  plt <- ggplot(na.omit(mlt), aes(x = Var1 * wdw, y = as.numeric(as.character(value)), col = Var2, group = Var2)) +
    geom_smooth(se = F, method = 'loess', formula = y ~ x) +
    labs(x = "Position", y = "Bendability", title = "Bendability per window",
         colour = NULL) +
    theme(plot.title = element_text(hjust = 0.5))

  return(list("graph" = plt, "table" = as.data.frame(bend_per_gene)))

}

# Benchmarking of the function

replicate(10, system.time(BendScoresGraph("sequence.txt", 20)))
replicate(10, system.time(BendScoresGraph_stringr("sequence.txt", 20)))


l <- BendScoresGraph("sequence.txt", 20)
BendScoresGraph_stringr("sequence.txt", 20)
k$table
all.equal(l$table,k$table)

