IUPAC <- function(fasta){

  code <- c("y", "r", "w", "s", "k", "m", "b", "d", "h", "v", "n", "x")
  anzVar <- c(2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4) 
  
  ambinuc <- data.frame(cbind(code, anzVar, a = FALSE, c = FALSE, g = FALSE, t = FALSE), stringsAsFactors = FALSE)
  
  ambinuc$a[c(2, 3, 6, 8:12)] <- TRUE 
  ambinuc$c[c(1, 4, 6, 7, 9:12)] <- TRUE 
  ambinuc$g[c(2, 4, 5, 7, 8, 10:12)] <- TRUE 
  ambinuc$t[c(1, 3, 5, 7:9, 11, 12)] <- TRUE 
  
  ambinuc$anzVar <- as.numeric(ambinuc$anzVar)
  
  seqn <- read.fasta(fasta, as.string = TRUE)
  
  seqn <- data.frame(cbind(seqnam = paste(">", names(seqn), sep = ""), seqn = seqn), stringsAsFactors = FALSE)
  seqn$seqn <- gsub(" ", "", seqn$seqn)
  
  for(i in 1:nrow(seqn)) {
    
    temp <- seqn$seqn[i]
    
    if (sum(is.element(ambinuc$code, unlist(strsplit(temp, split = "")))) == 0) {
      
      if(exists("fertig") == FALSE) fertig <- seqn[i,] else fertig <- rbind(fertig, seqn[i,])
    
      } 
    
    else {
      
      ersatz <- ambinuc[is.element(ambinuc$code, unlist(strsplit(temp, split = ""))), ]
      ersetzen <- unlist(strsplit(gsub("[a, c, g, t]", "", temp), split = ""))
      ersetzen <- data.frame(cbind(ersetzen, anzVar = NA), stringsAsFactors = FALSE)
      
      for(j in 1:nrow(ersetzen)) ersetzen$anzVar[j] <- ambinuc$anzVar[ambinuc$code == ersetzen$ersetzen[j]]
      
      ersetzen$anzVar <- as.numeric(ersetzen$anzVar)
      
      kombi <- prod(ersetzen$anzVar)
      disambseq <- data.frame(matrix(nrow = kombi, ncol = 2))
      colnames(disambseq) <- colnames(seqn)
      disambseq$seqn[1:kombi] <- temp
      disambseq$seqnam[1:kombi] <- paste(seqn$seqnam[i], rownames(disambseq), sep = "_")
      wdh <- kombi
      
      for(j in 1:nrow(ersetzen)) {
        
        wdh <- wdh/as.numeric(ersetzen$anzVar[j])
        wdhwdh <- kombi/(wdh*as.numeric(ersetzen$anzVar[j]))
        unambnuc <- rep(rep(names(which(t(ersatz[ersatz$code == ersetzen$ersetzen[j], ])[,1] == TRUE)), each = wdh), wdhwdh)
        
        for(k in 1:length(unambnuc)) disambseq$seqn[k] <- sub(ersetzen$ersetzen[j], unambnuc[k], disambseq$seqn[k])}
      
      if(exists("fertig") == FALSE) fertig <- disambseq else fertig <- rbind(fertig, disambseq)
    
    }
  }
  
  newfasta <- paste(dirname(fasta), paste("/disambig_", basename(fasta), sep = ""), sep = "")
  writeLines(paste(fertig$seqnam, fertig$seqn, sep = "\n"), newfasta)
  print(fertig)
}