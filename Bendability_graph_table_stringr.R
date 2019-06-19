library(ggplot2)
library(data.table)
library(Biostrings)
library(plyr)
library(stringr)

############################################################
### Function for reading in Multifasta of long sequences ###
############################################################

  BendScoresGraph <- function(mf, wdw) {

    ### Reading in the bendability table ###

    tabel <- fread("table_bendability_UPPERCASE.txt", sep = "\t", header = F)
    tabel <- as.data.table(tabel)
    colnames(tabel) <- c("Trinucl", "Score")
    tabel$Score <- as.numeric(gsub(",",".",tabel$Score))
    setkey(tabel, Trinucl)

    ### Reading in the Multifasta file ###

    mfasta <- readDNAStringSet(mf, "fasta")

    ### Average calculation per window function ###

    AvgCalc <- function(seq){

      st1 <- seq(1,nchar(as.character(seq)) - wdw + 1,wdw)
      en1 <- seq(wdw,nchar(as.character(seq)),wdw)
      z <- str_sub(seq, st1, en1)

      Windowmaker <- function(x){

        st2 <- seq(1, nchar(as.character(x)) - 2)
        en2 <- seq(3, nchar(as.character(x)))

        TriNucl <- str_sub(x, st2, en2)
        TriNucl <- as.data.table(TriNucl)
        colnames(TriNucl) <- "Trinucl"
        setkey(TriNucl, Trinucl)

        a <- merge(tabel, TriNucl)
        avg <- mean(as.numeric(a$Score))
        return(avg)
      }

      as.vector(unlist(sapply(z, Windowmaker)))

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

    ggsave(paste0(gsub("\\..*","",mf),".pdf"))

    return(list("graph" = plt, "table" = as.data.frame(bend_per_gene)))

  }

#############################################################
### Function for reading in Multifasta of short sequences ###
#############################################################

  BendScoresGraph_small <- function(mf){

    ### Reading in the bendability table and ###

    tabel <- fread("table_bendability_UPPERCASE.txt", sep = "\t", header = F)
    tabel <- as.data.table(tabel)
    colnames(tabel) <- c("Trinucl", "Score")
    tabel$Score <- as.numeric(gsub(",",".",tabel$Score))
    setkey(tabel, Trinucl)

    mfasta <- readDNAStringSet(mf, "fasta")

    bend_per_gene <- lapply(mfasta, function(x){

      st2 <- seq(1, nchar(as.character(x)) - 2)
      en2 <- seq(3, nchar(as.character(x)))

      TriNucl <- str_sub(x, st2, en2)
      TriNucl <- as.data.frame(TriNucl)
      colnames(TriNucl) <- "Trinucl"
      TriNucl <- as.data.table(TriNucl)
      a <- merge(TriNucl, tabel, by = "Trinucl", sort = F)
      return(a$Score)

    })

    bend_per_gene <- ldply(bend_per_gene, rbind)
    bend_per_gene <- t(bend_per_gene)
    colnames(bend_per_gene) <- bend_per_gene[1,]
    bend_per_gene <- bend_per_gene[-1,]
    rownames(bend_per_gene) <- c()

    mlt <- melt(bend_per_gene, measure.vars = colnames(bend_per_gene))

    plt <- ggplot(na.omit(mlt), aes(x = Var1 + 2, y = as.numeric(as.character(value)), col = Var2, group = Var2)) +
      #geom_smooth(se = F, method = 'loess', formula = y ~ x) +
      geom_line() +
      labs(x = "Position", y = "Bendability", title = "Bendability per window",
           colour = NULL) +
      theme(plot.title = element_text(hjust = 0.5))

    ggsave(paste0(gsub(".*","",mf),".pdf"))

  return(list("graph" = plt, "table" = as.data.frame(bend_per_gene)))

  }
