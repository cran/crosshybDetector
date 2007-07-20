"computeTm" <-
function(seqs, plot=FALSE){
  # Tm is computed using the approximated formula:
  # A -> 2
  # C -> 4
  # G -> 4
  # T -> 2
  s <- strsplit(toupper(seqs), "") 
  tmp <- sapply(s, function(z){ sum(z %in% c("A", "T"))*2 + sum(z %in% c("C", "G"))*4 })
  if (plot) { hist(tmp, main=paste("Histogram of Tm -- (n ", length(tmp), ")", sep="")) }
  return(tmp)
}

