\name{AAAlphabet-class}
\docType{class}
\alias{AAAlphabet-class}
\alias{levels,AAAlphabet-method}
\alias{show,AAAlphabet-method}

\title{Class "AAAlphabet" Amino Acid Alphabet }
\description{Chars alowed in a protein sequence. It is used by the
  method \code{initialize} of class AASequence to check if the chars in
  the sequence are known (recognized by the package).}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AAAlphabet", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"character"} Chars
      allowed in protein sequences}
    \item{\code{info}:}{Object of class \code{"character"} "AAAlphabet"}
  }
}
\section{Extends}{
Class \code{"character"}, directly.
}
\section{Methods}{
  \describe{
    \item{levels}{\code{signature(x = "AAAlphabet")}:  internal}
    \item{show}{\code{signature(object = "AAAlphabet")}:  show}
  }
}
\references{\url{http://r4proteomics.sourcefroge.net}}
\author{ Witold E. Wolski \email{witek96@users.sourceforge.net}}
%\note{ ~~further notes~~ }

\seealso{
}
\examples{
new("AAAlphabet")
}
\keyword{classes}
