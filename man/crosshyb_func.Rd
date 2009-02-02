\name{crosshyb_func}
\alias{crosshyb_func}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Computing cross-hybrization pvalues }
\description{
  It is usually called from within \code{crosshyb}, but can be invoked directly.
  The algorithm is described into the paper in \code{References} section 
}
\usage{
crosshyb_func(input, seqList, mysub, numPermut,
              probes, satValue, maxProbes, delta = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{ A dataframe containing the i) raw intensity value, 
  ii) the probe type (e.g. spike, probes, negative,...) and
  iii) the probe name }
  \item{seqList}{ An object of class \code{AASequenceList} created with \code{readFasta} }
  \item{mysub}{ Substitution matrix to be used for the alignment }
  \item{numPermut}{ Number of permutations for Monte Carlo simulations }
  \item{probes}{ A vector with the probes to be analyzed (e.g. probes, spike, etc) }
  \item{satValue}{ The value corresponding to saturated probes }
  \item{maxProbes}{ The maximun number of probes to analyze }
  \item{delta}{ Threshold for difference in melting temperature. Pairs of probes whose delta Tm
  is lower than \code{delta} are considered as similar. Delta Tm i,k is computed as Tm i,i - Tm i,k,
  where Tm i,i is the Tm of the percect pairing probe i vs probe i and Tm i,k is the Tm of the 
  imperfect pairing probe i vs probe k
  }
}
\value{
  A list containing:
  \item{data}{Dataframe with three columns containing the results of the analysis for 
  the putative corruptor probes. Contains: i) the probe number, ii) the probe name and
  iii) the pvalue based on Monte Carlo simulations of being a corruptor (pvalue is corrected for
  multiple testing)}
  \item{children}{A list containing the corrupted probes for each putative corruptor}
}
\references{ Uva, P., and de Rinaldis, E. CrossHybDetector: detection of cross-hybridization events in DNA microarray experiments. BMC Bioinformatics 2008, 9:485 }
\author{ Paolo Uva }
\seealso{ \code{\link{crosshyb}} }
\keyword{ misc }
