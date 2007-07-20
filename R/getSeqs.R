"getSeqs" <-
function(file,start=NULL,stop=NULL){
  seqs     <- read.delim(file, header=FALSE, as.is=TRUE)
  
  # Cut the sequences if needed
  if((!is.null(start))&(!is.null(stop))){
    seqs[,2] <- substr(seqs[,2], start=start, stop=stop)
  }
  return(seqs)
}

