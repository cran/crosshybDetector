"crosshyb2xls.putative" <-
function(input, arrayName){

  preName  <- gsub("\\.gpr$|_FEATURES.txt$","",arrayName)
  # Write putative corruptors RED channel
  write.table(input$dataR,
              file=paste(preName, "_putative_corruptorsRed.xls", sep=""),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  # Write putative corruptors GREEN channel
  write.table(input$dataG,
              file=paste(preName, "_putative_corruptorsGreen.xls", sep=""),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  # Write putative corrupted RED/GREEN channel
  list2ascii(input$childrenR, file=paste(preName, "_putative_corruptedRed.xls", sep=""))
  list2ascii(input$childrenG, file=paste(preName, "_putative_corruptedGreen.xls", sep=""))
}

