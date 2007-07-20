"extractBadProbes" <-
function(input, pVal=0.001){

   badProbes    <- list(corruptorsR=vector(), corruptorsG=vector(), corruptorsRG=vector(),
                        corruptedR=vector(), corruptedG=vector(), corruptedRG=vector())

   # RED channel
   if (sum(input$dataR[,3] < pVal) > 0){
      child  <- unique(unlist(input$childrenR[input$dataR[,3] < pVal]))
      parent <- input$dataR[(input$dataR[,3] < pVal),1]
      
      # Update list of parents and childrens
      badProbes$corruptorsR <- parent
      badProbes$corruptedR <- child
    }
   
   # GREEN channel
   if (sum(input$dataG[,3] < pVal) > 0){
      child  <- unique(unlist(input$childrenG[input$dataG[,3] < pVal]))
      parent <- input$dataG[(input$dataG[,3] < pVal),1]
      
      # Update list of parents and childrens
      badProbes$corruptorsG <- parent
      badProbes$corruptedG <- child
    }

   # Extract overlap
   badProbes$corruptorsRG <- intersect(badProbes$corruptorsR, badProbes$corruptorsG)
   badProbes$corruptedRG  <- intersect(badProbes$corruptedR, badProbes$corruptedG)
   
   return(badProbes)
}

