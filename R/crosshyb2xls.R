"crosshyb2xls" <-
function(raw, array=NULL, parent, children, arrayName, channel=c("red", "green"), probeNameID=c("ProbeName", "Name")){
  
  # Define variables
  if (channel == "red"){
    chName <- "Red"
    maf    <- maRf(raw[,array])
    mab    <- maRb(raw[,array]) 
  }else{
    chName <- "Green"
    maf    <- maGf(raw[,array])
    mab    <- maGb(raw[,array])
  }
  
  preName  <- gsub("\\.gpr$|_FEATURES.txt$","",arrayName)
  colNames <- c("ProbeNumber", 
                "ProbeName", 
                paste(chName, "Foreground_raw", sep=""), 
                paste(chName, "Background_raw", sep=""))
  
  # Create parent dataframe
  parent.df <- data.frame(parent, 
                          I(maInfo(maGnames(raw))[[probeNameID]][parent]), 
                          maf[parent], 
                          mab[parent])
  colnames(parent.df) <- colNames
  write.table(parent.df, file=paste(preName, "_corruptors", chName, ".xls", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE )
  
  # Create children dataframe    
  child.df  <- data.frame(child, 
                          I(maInfo(maGnames(raw))[[probeNameID]][child]),
                          maf[child],
                          mab[child])
  colnames(child.df)  <- colNames
  write.table(child.df, file=paste(preName, "_corrupted", chName, ".xls", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE )
}

