\name{crosshybMCplot}
\alias{crosshybMCplot}
\title{ Plot of Monte Carlo simulations }
\description{
  Plot the results of Monte Carlo simulations.
}
\usage{
crosshybMCplot(input, pVal = 0.001, arrayName = NULL, doPlot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{ the output of \code{crosshyb} containing the pvalues
  obtained from Monte Carlo simulations }
  \item{pVal}{ if the pvalue of a probe is < \code{pVal}, the probe is flagged as corruptor }
  \item{arrayName}{ the name of the array, used to create the output file name }
  \item{doPlot}{ logical, if TRUE writes the plot to file }
}
\details{
  This function draws two scatter plots, one for each channel (red and green), containing the 
  pvalue from Monte Carlo simulations ordered by raw intensity. A red line corresponding to
  the \code{pVal} parameter is added. All the probes above the red line (i.e. whose pvalue < \code{pVal})
  are considered as corruptors and their number is shown in the title. Such a number reflects the amount of
  cross-hybridization in the microarray experiment.\cr
  The degree of the cross-hybridization effect is proportional to the relative
  abundance of the non-target sequence versus the target sequence. As a consequence, 
  abundant target sequences can generate high signal intensities on their related chip spots
  as well as increasing the intensity values of spots carrying similar probes. For this reason most of the 
  probes identified as corruptors fall into the left part of the plot corresponding to those probes with the
  higher raw intensity signals.
}
\value{
  If \code{doPlot} is FALSE the plot is created on the current graphics device.
  If \code{doPlot} is TRUE the plot is only written to file.
}
\author{ Paolo Uva }
\seealso{ \code{\link{crosshyb}} }
\examples{
\dontrun{
# Run crosshyb algorithm...
# This function will take several minutes to finish
data(raw)
data(probeSeq)
crosshyb.out <- crosshyb(raw, probeSeq, plate=1, numPermut=10000, probeNameID="Name",
                         probes=c("probes", "spike"), satValue = 65535, maxProbes=100)                      
}

# ... or load directly the crosshyb output
data(crosshyb.out)

# Plot p-values obtained with Monte Carlo simulations
crosshybMCplot(crosshyb.out, pVal=0.01, arrayName="myArray", doPlot=FALSE)

}
\keyword{ misc }
