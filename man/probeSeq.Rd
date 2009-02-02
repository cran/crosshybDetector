\name{probeSeq}
\alias{probeSeq}
\docType{data}
\title{ Sequence data}
\description{
  A vector of sequence data.
}
\usage{data(raw)}
\format{
  A vector containing 25760 sequences corresponding to the probes in \code{raw} dataset.
  Probes whose sequence is NA will not be used by \code{crosshyb}
}
\references{ Uva, P., and de Rinaldis, E. CrossHybDetector: detection of cross-hybridization events in DNA microarray experiments. BMC Bioinformatics 2008, 9:485 }

\examples{
data(probeSeq)
str(probeSeq)
}
\keyword{datasets}
