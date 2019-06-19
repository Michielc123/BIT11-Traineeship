library("ggplot2")
library("data.table")
library("Biostrings")
library("dplyr")
library("plyr")

### Reading in the bendability table ###

setwd("C:/Users/Michiel/Desktop/Sequences")
tabel <- fread("table_bendability_UPPERCASE.txt", sep = "\t", header = F)
tabel <- as.data.table(tabel)
colnames(tabel) <- c("Trinucl", "Score")
tabel$Score <- as.numeric(gsub(",",".",tabel$Score))
setkey(tabel, Trinucl)

### Defining the function ###

AvgCalc <- function(seq, wdw){
  
  df <- data.frame(Window = "", Bendability = "", stringsAsFactors = FALSE)
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

### Reading in the Multifasta file

mfasta <- readDNAStringSet("multifasta_same.txt", "fasta")
bend_per_gene <- lapply(mfasta, AvgCalc, wdw = 15)

bend_per_gene <- ldply(bend_per_gene, rbind)
bend_per_gene <- t(bend_per_gene)
colnames(bend_per_gene) <- bend_per_gene[1,]
bend_per_gene <- bend_per_gene[-1,]
bend_per_gene <- as.data.frame(bend_per_gene)
rownames(bend_per_gene) <- c()

mlt <- melt(bend_per_gene, measure.vars = colnames(bend_per_gene))

mlt$seq <- rep(1:nrow(bend_per_gene),ncol(bend_per_gene))

ggplot(na.omit(mlt), aes(x = seq * 15, y = as.numeric(value), col = variable, group = variable)) + 
  geom_smooth(se = F, method = 'loess', formula = y ~ x) 



