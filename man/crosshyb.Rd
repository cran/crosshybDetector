\name{crosshyb}
\alias{crosshyb}
\title{ Identify probes affected by crosshybridization }
\description{
  This function identifies probes that are potentially affected by cross-hybridization.
  The probability of a probe of being a corruptor is computed by Monte Carlo simulations.
  The algorithm is briefly described into the \code{crosshyb_func} help page.
}
\usage{
crosshyb(raw, probeSeq, plate = 1, probeNameID = "ProbeName",
         numPermut = 10000, probes = c("probes", "spike"),
         satValue = 65535, maxProbes = 50, delta = 10)
}
\arguments{
  \item{raw}{an object of class \code{marrayRaw} }
  \item{probeSeq}{a vector with the probe sequences }
  \item{plate}{integer, the array to analyze }
  \item{probeNameID}{the slot name of \code{raw} containing the probe names.
  Probes names are extracted with \code{maInfo(maGnames(raw))[[probeNameID]]} }
  \item{numPermut}{integer, the number of permutations for Monte Carlo simulations }
  \item{probes}{the type of probes that are analyzed for potential corruption. 
  Probe types are defined in \code{maControls(raw)} }
  \item{satValue}{integer, the value corresponding to saturation (default 65535) }
  \item{maxProbes}{integer, the maximun number of probes to analyze for potential corruption }
  \item{delta}{integer, threshold for difference in melting temperature. Pairs of probes whose delta Tm
  is lower than \code{delta} are considered as similar. Delta Tm i,k is computed as Tm i,i - Tm i,k,
  where Tm i,i is the Tm of the perfect pairing probe i vs probe i and Tm i,k is the Tm of the 
  imperfect pairing probe i vs probe k
  }
}
\value{
  A list containing:
  \item{dataR}{Dataframe with three columns containing the results of the analysis for 
  the putative corruptor probes applied to the red channel. Contains: i) the probe number, ii) the probe name and
  iii) the pvalue based on Monte Carlo simulations of being a corruptor (pvalue is corrected for
  multiple testing)}
  \item{dataG}{Same as \code{dataR}, but applied to the green channel}
  \item{childrenR}{A list containing the corrupted probes for each putative corruptor from the red channel}
  \item{childrenG}{As \code{childrenR}, but applied to the green channel}
}
\author{ Paolo Uva }
\seealso{ \code{\link{crosshybMCplot}}, \code{\link{extractBadProbes}}, \code{\link{crosshyb_func}} }
\examples{
\dontrun{
data(raw)
data(probeSeq)
# This function will take several minutes to finish
crosshyb.out <- crosshyb(raw, probeSeq, plate=1, numPermut=10000, 
                         probeNameID="Name", probes=c("probes", "spike"),
                         satValue = 65535, maxProbes=100)
}
}
\keyword{ misc }
