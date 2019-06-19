library("ggplot2")
library("data.table")
library("Biostrings")
library("ggthemes")

### Reading in the bendability table

setwd("C:/Users/Michiel/Desktop")
tabel <- fread("table_bendability_UPPERCASE.txt", sep = "\t", header = F)
tabel <- as.data.table(tabel)
colnames(tabel) <- c("Trinucl", "Score")
tabel$Score <- as.numeric(gsub(",",".",tabel$Score))
setkey(tabel, Trinucl)

### Average Bendability calculation

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
    data.frame(Window = as.character(x),Bendability=avg,stringsAsFactors = FALSE)
    
    }
  
  sapply(view, Windowmaker)

}

test <- AvgCalc("ATATTATTAATTAATATCGCGGCGCGATGTCGTATATTACGCGCTGATGCTAGCTAGCTAGCATGCTAGCTAGCTGTTTTTTTTTTTTTTTTTTTCGATCGATCGATGCATCGATCGATGCTAGCTAATTATATATATATCGCGCTGACTGACTAGCTGATCGATCG", 20)
test <- t(test)
test <- as.data.frame(test)
test$Window <- as.character(test$Window)
test$Bendability <- as.numeric(test$Bendability)

## Plotting the Bendability

ggplot(test, aes(as.factor(sort(as.numeric(rownames(test)))), Bendability, group = 1)) + 
  stat_smooth(aes(colour = "#FF2400"), se = F , method = 'loess') + 
  geom_point(shape = 17, color = "#0E4D92", size = 2) + 
  ylab("Average Bendability") +
  scale_x_discrete(labels = paste(substr(test$Window,1,5), "...")) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "None", 
        axis.title.x = element_blank())



c <- sapply(b, AvgCalc, 200)
colnames(c) <- c("Average Bendability")
