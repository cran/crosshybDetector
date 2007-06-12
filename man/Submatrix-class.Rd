\name{Submatrix-class}
\docType{class}
\alias{Submatrix}
\alias{Submatrix-class}
%\alias{selfalign,AASequence,Submatrix-method}
\alias{show,Submatrix-method}
\alias{subFromEmboss,Submatrix,character-method}
\alias{subFromEmboss}

\title{Class "Submatrix" Substitution Matrix}
\description{Amino Acid Substitution Matrices}
\section{Objects from the Class}{
\code{new("Submatrix", ...)} or \code{data(EPAM100)}
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} substitution matrix }
    \item{\code{copyright}:}{Object of class \code{"character"} GNU GPL2 }
    \item{\code{info}:}{Object of class \code{"character"} name}
    \item{\code{head}:}{Object of class \code{"character"} description}
    \item{\code{alphabet}:}{Object of class \code{"character"} alphabet}
  }
}
\section{Extends}{
Class \code{"matrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{selfalign}{\code{signature(object = "AASequence", sub =
	"Submatrix")}: compute teh score for aligning the sequence with itself }
    \item{show}{\code{signature(object = "Submatrix")}: show}
    \item{subFromEmboss}{\code{signature(object = "Submatrix", path =
	"character")}: read Substitution matrices from EMBOSS package }
  }
}
\references{http://www.hgmp.mrc.ac.uk/Software/EMBOSS/overview.html}
\author{ Witold Eryk Wolski \url{http://r4proteomics.sourceforge.net} \email{witek96@users.sourceforge.net}}
\keyword{classes}
