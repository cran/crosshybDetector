"crosshybMAplot" <-
function(m, a, subset, badProbes, arrayName, doPlot=FALSE){
  
  plotName <- paste(gsub("\\.gpr$|_FEATURES.txt$","",arrayName), "_normMAWithBadProbes.png", sep = "")
  if (doPlot) {png(file=plotName, bg="white", width=600, height=600)}
   
  par(mfrow=c(1,1))
  plot(a[subset], m[subset], pch=16, cex=0.2, xlab="log2(average intensity)", ylab="log2(ratio)")
  title("MA plot - post-normalization", font.main=4, line=2)
  title(arrayName, font.main=2, line=1, cex.main=0.8)
    
  # Add saturated spots
  fields   <- c("corruptorsR", "corruptorsG", "corruptorsRG", "corruptedR", "corruptedG", "corruptedRG")
  cols     <- c("red4", "green4", "orange3", "red", "green", "orange")
  size     <- rep(c(1.5, 0.5), c(3,3))
  fieldPos <- match(fields, names(badProbes)) 
  for (i in 1:length(fields)){
    points(
          a[badProbes[[ fieldPos[i] ]] ],
          m[badProbes[[ fieldPos[i] ]] ],
          pch=16,
          cex=size[i],
          col=cols[i]
          )
  }
  # Add legend
  legend("topright", legend=fields,
         cex=0.8, pch=16, col=cols, pt.cex=size)  
              
   abline(h=1)
   abline(h=-1)
   if(doPlot) {dev.off()}
}

