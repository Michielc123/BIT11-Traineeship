
# Function used to merge a list of dataframes consisting of small sequences. Merging only the the ones who overlap.
# Output: List of all possible sequences with a 2 nucleotide overlap

library(ggplot2)
library(data.table)
library(Biostrings)
library(plyr)
library(stringr)
library(DECIPHER)

setwd("C:/Users/Michiel/Desktop/Sequences")
load("example_longerseq.RData")

l <- lapply(l, function(x){x$sequence})
l
k <- list(l[[3]], l[[4]])
k

func2 <- function(c){
  
  func <- function(x,y){
      a <- grep(str_sub(x,nchar(x)-4, nchar(x)), y, value = T)
      if(length(a) != 0){paste0(str_sub(x,1,nchar(x)-5),a)}
    }
      
  c[[1]] <- unlist(lapply(c[[1]],func,y = c[[2]]))
  c[[2]] <- NULL
  return(c)
}

while(length(l) > 1){
  l <- func2(l)
  length(l[[1]])
}



