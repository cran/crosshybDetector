\name{crosshybImage}
\alias{crosshybImage}
\title{ Spatial distribution of corruptor and corrupted probes }
\description{
  Draws an image with the spatial distribution on the array of probes flagged as corruptors or corrupted.
}
\usage{
crosshybImage(raw, plate, parent, children, arrayName,
              channel = c("red", "green"), doPlot = FALSE)
}
\arguments{
  \item{raw}{ an object of class \code{marrayRaw} }
  \item{plate}{ the plate number }
  \item{parent}{ a vector of probe numbers corresponding to the corruptor probes }
  \item{children}{ a vector of probe numbers corresponding to the corrupted probes }
  \item{arrayName}{ the name of the array. Used to create the output file name }
  \item{channel}{ the channel to use for drawing image }
  \item{doPlot}{ logical, if TRUE writes the image to file }
}
\value{
  This function creates an image with the spatial distribution of the corruptor 
  and corrupted probes on the array for the selected channel. If \code{doPlot} is TRUE, the image is
  written to a file and nothing is displayed. The color intensity of corrupted probes is proportional
  to the intensity value
}
\author{ Paolo Uva }
\seealso{ \code{\link{extractBadProbes}} }
\examples{
data(raw)
data(probeSeq)

\dontrun{ 
# Run crosshyb...
# This function will take several minutes to finish
crosshyb.out  <- crosshyb(raw, probeSeq, plate=1, numPermut=10000,
                          name="Name", probes=c("probes", "spike"),
                          satValue = 65535, maxProbes=100)                  
}

# ... or load directly crosshyb output
data(crosshyb.out)

# Extract corrupted and corruptor probes  
badprobes <- extractBadProbes(crosshyb.out, pVal=0.01)

# Write results for RED channel
if(length(badprobes$corruptorsR)){
  parent <- badprobes$corruptorsR
  child  <- badprobes$corruptedR
  crosshybImage(raw, plate = 1, parent=parent, children=child,
                arrayName = "myArray", channel="red", doPlot=FALSE)
}

}
\keyword{ misc }
