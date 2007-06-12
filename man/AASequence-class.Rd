\name{AASequence-class}
\docType{class}
\alias{AASequence-class}
\alias{frequency,AASequence-method}
\alias{initialize,AASequence-method}
\alias{selfalign}
%\alias{readFasta,AASequence-method}
%\alias{salign,AASequenceList,AASequence-method}
%\alias{salign,AASequence,AASequenceList-method}
%\alias{salign,AASequence,AASequence-method}
\alias{selfalign,AASequence,Submatrix-method}
\alias{show,AASequence-method}
%\alias{testalign,AASequence,AASequence-method}

\title{"AASequence" Object representation of a protein Sequence}
\description{Object representation of a protein Sequence}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AASequence", sequence, info, alphabet)}.
or by call to \code{AASequence("PROTINSEQUENCEHERE",info="my prot seq")}
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"character"} the sequence}
    \item{\code{info}:}{Object of class \code{"character"} the info }
  }
}
\section{Extends}{
Class \code{"character"}, directly.
}
\section{Methods}{
  \describe{
    \item{frequency}{\code{signature(x = "AASequence")}: returns the
      amino acid frequencies }
    \item{initialize}{\code{signature(.Object = "AASequence")}: internal }
    \item{readFasta}{\code{signature(object = "AASequence")}: read FASTA
    file format}
    \item{salign}{\code{signature(obj1 = "AASequenceList", obj2 = "AASequence")}: see: \code{\link{salign-methods}}}
    \item{salign}{\code{signature(obj1 = "AASequence", obj2 = "AASequenceList")}: see: \code{\link{salign-methods}}}
    \item{salign}{\code{signature(obj1 = "AASequence", obj2 = "AASequence")}: see: \code{\link{salign-methods}}}
    \item{selfalign}{\code{signature(object = "AASequence", sub =
	"Submatrix")}: Computes Smith-Waterman score of the alignment of
    the sequence with itself.}
    \item{show}{\code{signature(object = "AASequence")}: }
    \item{testalign}{\code{signature(obj1 = "AASequence", obj2 =
	"AASequence")}: internal, see: \code{\link{testalign-methods}}}
  }
}
\references{\url{http://r4proteomics.sourceforge.net}}
\author{Witold E. Wolski \email{witek96@users.sourceforge.net}}
%\note{ ~~further notes~~ }

\seealso{
\code{\link{salign-methods}},\code{\link{readFasta-methods}}
}
\examples{
seq1<-new("AASequence","MEDQVGFGFRPNDEEL",info="seq1")
seq2<-AASequence("seq2","VAISEVNICSY")
}
\keyword{classes}
