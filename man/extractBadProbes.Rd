\name{extractBadProbes}
\alias{extractBadProbes}
\title{ Extract corruptor and corrupted probes }
\description{
  This function extracts the corruptor and corrupted probes for each channel
  from the output of \code{crosshyb}
}
\usage{
extractBadProbes(input, pVal = 0.001)
}
\arguments{
  \item{input}{ the output of \code{crosshyb} }
  \item{pVal}{ if the pvalue of a probe is < \code{pVal}, the probe is flagged as corruptor }
}
\value{
  A list containing:
  \item{corruptorsR }{index of probes flagged as corruptors based on red channel intensity}
  \item{corruptedR }{index of probes flagged as corrupted based on red channel intensity}
  \item{corruptorsG }{index of probes flagged as corruptors based on green channel intensity}
  \item{corruptedG }{index of probes flagged as corrupted based on green channel intensity}
  \item{corruptorsRG }{index of probes flagged as corruptors based on red AND green channel intensities
  - intersection of \code{corruptorsR} with \code{corruptorsG} }
  \item{corruptedRG }{index of probes flagged as corrupted based on red AND green channel intensities
  - intersection of \code{corruptedR} with \code{corruptedG} }
}
\author{ Paolo Uva }
\seealso{ \code{\link{crosshyb}} }
\examples{
\dontrun{
# Run crosshyb algorithm...
# This function will take several minutes to finish
data(raw)
data(probeSeq)
crosshyb.out <- crosshyb(raw, probeSeq, plate=1, numPermut=10000,
                         probeNameID="Name", probes=c("probes", "spike"),
                         satValue = 65535, maxProbes=100)                      
}

# ... or load directly the crosshyb output
data(crosshyb.out)

# Extract bad probes
badProbes <- extractBadProbes(crosshyb.out, pVal=0.01)

}
\keyword{ misc }
