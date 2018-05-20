


control <- "1000genom_chasm.tab"

ovMut <- read.delim("tcgaov_muts.tab",
                              header=T,	
                              row.names=1,
                              stringsAsFactors=F,
                              check.names=F)

ovMut <- as.data.frame(ovMut)
View(ovMut)