\name{testalign-methods}
\docType{methods}
\alias{testalign}
\alias{testalign-methods}
\alias{testalign,AASequence,AASequence-method}
\title{Scoring of optimal alingments.}
\section{Methods}{\describe{
\item{obj1 = "AASequence", obj2 = "AASequence"}{Computes obtimal
  alignment of two sequences and returns score.}
}}
\description{
  Computes optimal pairwise alignment of 2 protein sequences. For
  internal use.
}
\arguments{
  \item{obj1}{see above}
  \item{obj2}{see above}
  \item{sub}{Substitution matrix}
  \item{delta}{gap opening penalty}
  \item{gapext}{gap extension penalty}
  \item{alignment}{type of alignment: c("global","local","overlap")}
  \item{scoring}{type of scoring, identity, similarity, score, scoreN, pozitive}
}
\value{
  For two objects of class \code{AASequence} returns an object of class
  \code{AAAlignment}.\\
  Otherwise it returns a score. If used with NULL as second parameter it
  returns a distance matrix. The scores are transformed into
  distances. For details see the vignette.
}
\references{\url{http://r4proteomics.sourceforge.net}}
\author{Witold E. Wolski}
\keyword{methods}
