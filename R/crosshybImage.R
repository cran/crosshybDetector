"crosshybImage" <-
function(raw, plate, parent, children, arrayName, channel=c("red", "green"), doPlot=FALSE){
  
  if (channel == "red") { 
    postName <- "crosshybImageR"
    data2plot <- maRf(raw)[,plate] 
  }else{ 
    postName <- "crosshybImageG"
    data2plot <- maGf(raw)[,plate]
  }
  data2plot[-(children)] <- NA
  data2plot[parent]      <- - max(data2plot, na.rm=TRUE)
  
  plotName <- paste(gsub("\\.gpr$|_FEATURES.txt$","",arrayName), "_", postName, ".png", sep = "")
  if (doPlot){ png(file=plotName, bg="white", width=400, height=600) }
  maImage.func(data2plot, L = maLayout(raw),
               col= maPalette(low = "black", high = channel, mid="white", k = 30),
               zlim=range(data2plot, na.rm=TRUE),
               #main=paste("Spatial distribution of corruptors/corrupted --", channel, "channel")
               )
  legend("topright", legend=c(paste("Corruptors:", length(parent)), paste("Corrupted:", length(children))),
          cex=0.8, pch=16, col=c("black", channel))
  if(doPlot) { dev.off() }
}

