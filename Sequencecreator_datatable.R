# Function used to merge a list of dataframes consisting of small sequences. Merging only the the ones who overlap.
# Making use of a different approach

library(ggplot2)
library(data.table)
library(Biostrings)
library(dplyr)
library(stringr)
library(DECIPHER)

setwd("C:/Users/Michiel/Desktop/Sequences")
load("example_longerseq.RData")

l <- lapply(l, function(x){as.data.table(x$sequence)})

func4 <- function(x){
  
  func3 <- function(a, b){
    r <- str_sub(a, start = nchar(a) - 4, end = nchar(a))
    d <- empty(b[V1 %like% paste0("^",r)])
    return(d)
  }
  
  x[[1]][, log := unlist(lapply(x[[1]][,V1],func3, b = x[[2]]))]
  x[[1]] <- x[[1]][log == F]
  x[[1]]$log <- NULL
  
  return(x)
}

for (i in seq(1,length(l) - 1)){
  j = i + 1
  k <- l[i:j]
  print(k)
  k <- func4(k)
  l[[i]] <- k[[1]]
}

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