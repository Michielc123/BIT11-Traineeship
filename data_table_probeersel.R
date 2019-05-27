library(ggplot2)
library(data.table)
library(Biostrings)
library(dplyr)
library(stringr)
library(DECIPHER)

setwd("C:/Users/Michiel/Desktop/Sequences")
load("example_longerseq.RData")

l <- lapply(l, function(x){as.data.table(x$sequence)})

func2 <- function(k){

  func <- function(a, b){
  r <- str_sub(a, start = nchar(a) - 4, end = nchar(a))
  d <- b[V1 %like% paste0("^",r)]
  s <- as.data.table(lapply(d, function(x){paste0(str_sub(a,1,nchar(a)-5), x)}))

  }

  k[[1]] <- rbindlist(lapply(k[[1]][,V1], func, b = k[[2]]))
  k[[1]] <- k[[1]][nchar(V1) == max(nchar(V1))]
  k[[2]] <- NULL
  return(k)
  
}


while(length(l > 1)){
  l <- func2(l)
  print(length(l))
}

