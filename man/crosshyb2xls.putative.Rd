\name{crosshyb2xls.putative}
\alias{crosshyb2xls.putative}
\title{ Writing list of analyzed probes to file }
\description{
  Writes the lists of probes processed by \code{crosshyb} to files. 
  Those files are only used for tracking purposes and contain the ASCII
  representation of \code{crosshyb} output.
}
\usage{
crosshyb2xls.putative(input, arrayName)
}
\arguments{
  \item{input}{ The output of \code{crosshyb} }
  \item{arrayName}{ The name of the array. Used to create the output file name }
}
\value{
  This function writes a tab-delimited files for each channel of a dual channel
  microarray experiment into the current directory
}
\author{ Paolo Uva }
\seealso{ \code{\link{crosshyb}}, \code{\link{crosshybMCplot}} }
\examples{
\dontrun{
# Run crosshyb algorithm...
# This function will take several minutes to finish
data(raw)
data(probeSeq)
crosshyb.out <- crosshyb(raw, probeSeq, plate=1, numPermut=10000,
                         probeNameID="Name", probes=c("probes", "spike"),
                         satValue = 65535, maxProbes=100)                      

# ... or load directly the crosshyb output
data(crosshyb.out)

# Write probes analyzed by crosshyb to files
crosshyb2xls.putative(crosshyb.out, arrayName="myArray")
}
}
\keyword{ misc }
