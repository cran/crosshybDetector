\name{AAAlignment-class}
\docType{class}
\alias{AAAlignment-class}
\alias{show,AAAlignment-method}
\alias{summary,AAAlignment-method}

\title{"AAAlignment" result of pairwise sequence alignment.}
\description{Result of a pairwise sequence alignment.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AAAlignment", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{info1}:}{Object of class \code{"character"} :  name of
      first sequence. }
    \item{\code{info2}:}{Object of class \code{"character"} : name of
      second sequence }
    \item{\code{selfs1}:}{Object of class \code{"numeric"} : selfscore
      of sequence 1 }
    \item{\code{selfs2}:}{Object of class \code{"numeric"} : selfscore
      of sequence 2}
    \item{\code{score}:}{Object of class \code{"numeric"} : Smith
      Watermann score of alignment. }
    \item{\code{identity}:}{Object of class \code{"numeric"} : nr of
      identitacal AA in the alignment.}
    \item{\code{alignsimilarity}:}{Object of class \code{"numeric"} : nr
    of similarities.}
    \item{\code{lch1}:}{Object of class \code{"numeric"} : length of
      sequence 1.}
    \item{\code{lch2}:}{Object of class \code{"numeric"} : length of
      sequence 2 }
    \item{\code{alig1}:}{Object of class \code{"character"} : alignment string }
    \item{\code{alig2}:}{Object of class \code{"character"} : alignment string }
    \item{\code{beautify}:}{Object of class \code{"character"} :
      internal }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "AAAlignment")}: show }
    \item{summary}{\code{signature(object = "AAAlignment")}: summary }
  }
}
\references{\url{http://r4proteomics.sourceforge.net}}
\author{Witold E. Wolski \email{witek96@users.sourceforge.net}}
%\note{ ~~further notes~~ }

\seealso{
  \code{\link{salign-methods}} 
}
\examples{
data(EBLOSUM62)
seq1<-new("AASequence","MEDQVGFGFRPNDEEL",info="seq1")
seq2<-AASequence("seq1","VAISEVNICSYDPWNL")
res<-salign(seq1,seq2,EBLOSUM62,delta=-10,gapext=-1,alignment="global")
class(res)
summary(res)
}
\keyword{classes}
