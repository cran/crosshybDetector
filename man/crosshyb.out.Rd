\name{crosshyb.out}
\alias{crosshyb.out}
\docType{data}
\title{Output of crosshyb}
\description{
  This list contains an example of the output of \code{crosshyb}.
}
\usage{data(crosshyb.out)}

\value{
  \item{dataR}{Dataframe with three columns containing the results of the analysis for 
  the putative corruptor probes applied to the red channel. Contains: i) the probe number, ii) the probe name and
  iii) the pvalue based on Monte Carlo simulations of being a corruptor (pvalue is corrected for
  multiple testing)}
  \item{dataG}{Same as \code{dataR}, but applied to the green channel}
  \item{childrenR}{A list containing the corrupted probes for each putative corruptor from the red channel}
  \item{childrenG}{As \code{childrenR}, but applied to the green channel}
}

\references{ Uva, P., and de Rinaldis, E. CrossHybDetector: detection of cross-hybridization events in DNA microarray experiments. BMC Bioinformatics 2008, 9:485 }
\examples{
data(crosshyb.out)
str(crosshyb.out)
}
\keyword{datasets}
