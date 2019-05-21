library(ggplot2)
library(data.table)
library(Biostrings)
library(plyr)
library(stringr)
library(DECIPHER)

PlotIt <- function(a){

  a$deltaBC <- NULL
  a$strdist <- NULL
  mlt <- melt(a, "sequence")

  
  plt <- ggplot(mlt, aes(x = variable, y = as.numeric(as.character(value)), col = sequence, group = sequence)) +
    geom_line() +
    labs(x = "\nPosition", y = "Bendability\n",
         colour = NULL) +
    theme(plot.title = element_text(hjust = 1),
          legend.text = element_text(size = 17),
          legend.key.size = unit(0.8,"cm"),
          plot.margin = unit(c(1,3,3,3), "cm"),
          legend.margin = margin(0,0,0,100),
          axis.text.x = element_text(size = 30, margin = margin(c(20,0,0,0))),
          axis.text.y = element_text(size = 30, margin = margin(c(0,20,0,0))),
          axis.title = element_text(size = 30)) +
    scale_x_discrete(labels = seq(1,ncol(a) - 1)) +
    guides(col = guide_legend(ncol = 2))
  
  ggsave("PlotIt.pdf", width = unit(27, "cm"), height = unit(19,"cm"))

}
