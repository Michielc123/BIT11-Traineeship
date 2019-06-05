library(ggplot2)
library(data.table)
library(Biostrings)
library(plyr)
library(stringr)
library(DECIPHER)

setwd("C:/Users/Michiel/Desktop/Sequences")
load("example_longerseq.RData")

l <- lapply(l, function(x){as.data.table(x$sequence)})
l
k <- list(l[[3]], l[[4]])
k

func <- function(x,y){
  a <- grep(str_sub(x,nchar(x)-4, nchar(x)), y, value = T)
  if(length(a) != 0){paste0(str_sub(x,1,nchar(x)-5),a)}
}

func(k[[1]], k[[2]])

func2 <- function(c){
  
  func <- function(x,y){
    a <- grep(str_sub(x,nchar(x)-4, nchar(x)), y, value = T)
    if(length(a) != 0){paste0(str_sub(x,1,nchar(x)-5),a)}
  }
  
  unlist(lapply(c[[1]],func,y = c[[2]]))
  
}

func2(k)

while(length(l) > 1){
  l <- func2(l)
  length(l[[1]])
}

