\name{crosshyb2xls}
\alias{crosshyb2xls}
\title{ Write corruptor and corrupted probes to files }
\description{
  This function writes \code{corruptor} and \code{corrupted} probes to tab-delimited files.
  Each file contains the probe names and the corresponding raw (before normalization) foreground and background signals.
}
\usage{
crosshyb2xls(raw, array = NULL, parent, children, arrayName,
             channel = c("red", "green"), probeNameID = c("ProbeName", "Name"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw}{ An object of class \code{marrayRaw} }
  \item{array}{ Integer. The array to analyze }
  \item{parent}{ A vector with the index of the corruptor probes, as returned 
                 by \code{extractBadProbes} }
  \item{children}{ A vector with the index of the corrupted probes, as returned
                   by \code{extractBadProbes} }
  \item{arrayName}{ The name of the array to analyze. It is used to create the output file names }
  \item{channel}{ The channel to analyze }
  \item{probeNameID}{ The slot name of \code{raw} containing the probe names.
  Probes names are extracted with \code{maInfo(maGnames(raw)[[probeNameID]]} }
}

\value{
  Write a tab-delimited files for each channel of a dual channel microarray experiment
  into the current directory. Each file has four columns:
  \code{ProbeNumber}, \code{ProbeName}, \code{Foreground\_Raw} and \code{Background\_Raw}
  containing respectively the probe number, probe name and the foreground and
  background raw signals (before normalization) 
}
\author{ Paolo Uva }

\seealso{ \code{\link{extractBadProbes}}, \code{\link{crosshyb}}, \code{\link{crosshybImage}} }
\examples{
data(raw)
data(probeSeq)

\dontrun{ 
# Run crosshyb...
# This function will take several minutes to finish
crosshyb.out  <- crosshyb(raw, probeSeq, plate=1, numPermut=10000, name="Name",
                          probes=c("probes", "spike"), satValue = 65535,
                          maxProbes=100)                  
}

# ... or load directly crosshyb output
data(crosshyb.out)

# Extract corrupted and corruptor probes  
badprobes <- extractBadProbes(crosshyb.out, pVal=0.05)

\dontrun{                   
# Write results for RED channel
if(!length(badprobes$corruptorsR)){
  parent <- badprobes$corruptorsR
  child  <- badprobes$corruptedR
  crosshyb2xls(raw, array=1, parent=parent, children=child,
               arrayName="myArray", channel="red", probeNameID="Name")
}
}

}
\keyword{ misc }
