"crosshybMCplot" <-
function(input, pVal=0.001, arrayName=NULL, doPlot=FALSE){
  if (is.null(arrayName)){
    print("Please select the array name!")
    return(FALSE)
  }

  plotName <- paste(gsub("\\.gpr$|_FEATURES.txt$","",arrayName), "_montecarloPvalue.png", sep = "")
  if (doPlot){ png(file=plotName, bg="white", width=400, height=600) }

  par(mfrow=c(2,1))
  if (prod(dim(input$dataR)) != 0){
    plot(-log10(input$dataR[,3]), xlab="Rank by Red raw intensity", ylab="-log10(p-value)")
    abline(h=-log10(pVal), col="red")
    title(paste("Corruptors - n: ", sum(input$dataR[,3] < pVal), " (p < ", pVal, ")", sep="" ), font.main=4, line=2)
    title(paste("Red channel -", arrayName), font.main=2, line=1, cex.main=0.8)
  }
  if (prod(dim(input$dataG)) != 0){
    plot(-log10(input$dataG[,3]), xlab="Rank by Green raw intensity", ylab="-log10(p-value)")
    abline(h=-log10(pVal), col="red")
    title(paste("Corruptors - n: ", sum(input$dataG[,3] < pVal), " (p < ", pVal, ")", sep="" ), font.main=4, line=2)
    title(paste("Green channel -", arrayName), font.main=2, line=1, cex.main=0.8)
  }
  if(doPlot){ dev.off() }
}

