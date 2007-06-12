\name{crosshybDetector-package}
\alias{crosshybDetector-package}
\alias{crosshybDetector}
\docType{package}
\title{
Detection of cross-hybridization events in microarray experiments
}
\description{
\code{crosshybDetector} is a package which calculates for each probe on the array the
probability of cross-hybridization by using the probe intensity values. The software provides
the user with the list of probes potentially corrupted and the associated p-values computed
by Monte Carlo simulations. Plots are generated allowing a global overview of the cross-hybridization
events in a microarray experiments. The package contains code from the dismissed pairseqsim package
 (removed in BioConductor > 1.9 ), created by Witold Wolski <witek96@users.sourceforge.net> and
  released under GNU LGPL license
}
\details{
\tabular{ll}{
Package: \tab crosshybDetector\cr
Type: \tab Package\cr
Version: \tab 1.0.3\cr
Date: \tab 2008-06-12\cr
License: \tab LGPL\cr
}
}
\author{
Paolo Uva <paolo\_uva@merck.com>
}
\references{
Uva, P., and de Rinaldis, E. CrosshybDetector: detection of cross-hybridization
events in DNA microarray experiments. Submitted.
}
\keyword{ package }
\seealso{
\code{\link{crosshyb}}
}
\examples{
\dontrun{
# This workflow will produce several images and files
# containing the results of crosshybDetector

pVal <- 0.01        # Threshold for corruptors
data(probeSeq)      # Vector of probe sequences
data(raw)           # Object of class marrayRaw containing one array

# Array names
arrayNames <- c("myArray")
for (k in 1:length(arrayNames)){
      
  # Run crosshyb
  crosshyb.out <- crosshyb(raw, probeSeq, plate=k, numPermut=10000,
                           probeNameID="Name", probes="probes",
                           satValue = 65535, maxProbes=50)                  
    
  # Write probes analyzed by crosshyb to file
  # Not nice for reading but useful for tracking
  crosshyb2xls.putative(crosshyb.out, arrayName=arrayNames[k])
    
  # Plot p-values obtained with Monte Carlo simulations
  crosshybMCplot(crosshyb.out, pVal, arrayName=arrayNames[k], doPlot=TRUE)
    
  # Extract bad probes  
  badProbes <- extractBadProbes(crosshyb.out, pVal)
                      
  # Plots for RED channel
  if(length(badProbes$corruptorsR)){
    parent <- badProbes$corruptorsR
    child  <- badProbes$corruptedR
      
    # Draw image plot
    crosshybImage(raw, plate = k, parent=parent, children=child, 
                  arrayName = arrayNames[k], channel="red", doPlot=TRUE)
      
    # Write parents and childrens to tab-delimited files (xls extension)
    crosshyb2xls(raw, array=k, parent=parent, children=child,
                 arrayName=arrayNames[k], channel="red", probeNameID="Name")
  }
    
  # Plots for GREEN channel
  if(length(badProbes$corruptorsG)){
    parent <- badProbes$corruptorsG
    child  <- badProbes$corruptedG
      
    # Draw image plot
    crosshybImage(raw, plate = k, parent=parent, children=child,
                  arrayName = arrayNames[k], channel="green", doPlot=TRUE)
      
    # Write parents and childrens to tab-delimited files (xls extension)
    crosshyb2xls(raw, array=k, parent=parent, children=child,
                 arrayName=arrayNames[k], channel="green", probeNameID="Name")
  }
   
  # Loess normalization using only probes
  # Use maNorm function from marray package
  norm <- maNorm(raw, norm="l", subset=maControls(raw) == "probes")

  # MA plot showing corruptor and corrupted probes
  crosshybMAplot(m = maM(norm[,k]),
                 a = maA(norm[,k]), 
                 subset=maControls(raw) \%in\% c("probes", "spike"),
                 badProbes=badProbes,
                 arrayName=arrayNames[k],
                 doPlot=TRUE) 
  rm(badProbes)
}
}

}
