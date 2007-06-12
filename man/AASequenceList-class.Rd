\name{AASequenceList-class}
\docType{class}
\alias{AASequenceList-class}
\alias{[,AASequenceList-method}
\alias{[[<-,AASequenceList-method}
\alias{[<-,AASequenceList-method}
%\alias{readFasta,AASequenceList-method}
%\alias{salign,AASequenceList,NULL-method}
%\alias{salign,AASequenceList,AASequence-method}
%\alias{salign,AASequence,AASequenceList-method}
\alias{show,AASequenceList-method}
\title{Class "AASequenceList" - Container}
\description{Container for objects of calss AASequence. The AASequences
  in the list can be accessed and manipulated like normal list elements.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AASequenceList", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} : list of
      AASequence objects.}
    \item{\code{info}:}{Object of class \code{"character"} : info character}
  }
}
\section{Extends}{
Class \code{"list"}, directly.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "AASequenceList")}: create subset AASequenceList }
    \item{[[<-}{\code{signature(x = "AASequenceList")}: assign AASequence}
    \item{[<-}{\code{signature(x = "AASequenceList")}: replace subset of
    AASequenceList}
    \item{readFasta}{\code{signature(object = "AASequenceList")}:
      read file in FASTA format. }
    \item{salign}{\code{signature(obj1 = "AASequenceList", obj2 =
	"NULL")}: see: \code{salign}}
    \item{salign}{\code{signature(obj1 = "AASequenceList", obj2 = "AASequence")}:  see: \code{salign} }
    \item{salign}{\code{signature(obj1 = "AASequence", obj2 = "AASequenceList")}:  see: \code{salign} }
  }
}
\references{\url{http://users.sourceforge.net}}
\author{ Wolski E. Witold \email{witek96@users.sourceforge.net}}
%\note{ ~~further notes~~ }

\seealso{
  \code{\link{salign-methods}}, or \code{\link{readFasta-methods}}
}
\examples{
}
\keyword{classes}
