\name{salign-methods}
\docType{methods}
\alias{salign}
\alias{salign-methods}
\alias{salign,AASequenceList,NULL-method}
\alias{salign,AASequenceList,AASequence-method}
\alias{salign,AASequence,AASequenceList-method}
\alias{salign,AASequence,AASequence-method}
\title{Optimal pairwise alignment of two protein sequences.}
\section{Methods}{\describe{
\item{obj1 = "AASequenceList", obj2 = "NULL"}{Compute distance matrix.}
\item{obj1 = "AASequenceList", obj2 = "AASequence"}{Compute scores.}
\item{obj1 = "AASequence", obj2 = "AASequenceList"}{Compute scores.}
\item{obj1 = "AASequence", obj2 = "AASequence"}{Compute alignment.}
}}
\description{
  Computes optimal pairwise alignment of 2 protein sequences.
}
\arguments{
  \item{obj1}{see above}
  \item{obj2}{see above}
  \item{sub}{Substitution matrix}
  \item{delta}{gap opening penalty}
  \item{gapext}{gap extension penalty}
  \item{alignment}{type of alignment: c("global","local","overlap")}
  \item{scoring}{type of scoring, identity, similarity, score, scoreN,
    AAAlignemnt, pozitive}
}
\value{
  For two objects of class \code{AASequence} returns an object of class
  \code{AAAlignment}.  Otherwise it returns a score. If used with \code{NULL} as second parameter it
  returns a distance matrix. The scores than are transformed into distances. For details see the vignette.
}
\references{Hilary S Booth, John H Maindonald, Ole M Nielsen, Susan R
  Wilson "Normalizing sequence alignment scores for composition bias"
  Recomb 2003 - Berlin}
\author{Witold E. Wolski \email{witek96@users.sourceforge.net} \url{http://r4proteomics.sourceforge.net}}
\examples{
}
\keyword{methods}

