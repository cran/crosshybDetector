\name{crosshybMAplot}
\alias{crosshybMAplot}
\title{ MA-plot }
\description{
  Creates an MA-plot with color coding for corruptor and corrupted probes.
}
\usage{
crosshybMAplot(m, a, subset, badProbes, arrayName, doPlot = FALSE)
}
\arguments{
  \item{m}{ a vector containing log2(ratios) }
  \item{a}{ a vector containing log2(average intensities) }
  \item{subset}{ the subset of probes to be plotted }
  \item{badProbes}{ the output of \code{extractBadProbes} containing corruptor
  and corrupted probes for red and green channels }
  \item{arrayName}{ the name of the array. Used to create the output file name }
  \item{doPlot}{ logical, if TRUE writes the plot to file }
}
\details{
  An MA-plot is a plot of log-intensity ratios (M-values) versus
  log-intensity averages (A-values)
}
\value{
  If \code{doPlot} is FALSE a plot is created on the current graphics device.
  If \code{doPlot} is TRUE the plot is only written to file.
}
\author{ Paolo Uva }
\seealso{ \code{\link{extractBadProbes}} }
\examples{
data(raw)
data(probeSeq)
\dontrun{
# Run crosshyb - This function will take several minutes to finish
crosshyb.out <- crosshyb(raw, probeSeq, plate=1, numPermut=10000, probeNameID="Name",
                         probes="probes", satValue = 65535, maxProbes=100)                  
    
}

# ...or load directy crosshyb output
data(crosshyb.out)  

# Extract bad probes  
badProbes <- extractBadProbes(crosshyb.out, pVal=0.01)
   
# Loess normalization using only probes
# Use maNorm function from marray package
norm <- maNorm(raw, norm="l", subset=maControls(raw) == "probes")

# MA plot showing corruptor and corrupted probes
crosshybMAplot(m = maM(norm),
               a = maA(norm), 
               subset=maControls(raw) \%in\% c("probes", "spike"),
               badProbes=badProbes,
               arrayName="myArray",
               doPlot=FALSE) 
}
\keyword{ misc }
