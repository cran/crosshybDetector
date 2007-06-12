\name{listdist-methods}
\docType{methods}
\alias{listdist}
\alias{listdist-methods}
\alias{listdist,list-method}
\title{Compute distance matrix from a list}
\section{Methods}{
\describe{
\item{object = "list"}{Compute distance matrix from a list.}
}}
\description{
  Compute distance matrix from list. Internal use.
}
\arguments{
  \item{object}{object of class \code{list}}
  \item{FUN}{Function to apply to list entries}
  \item{diag}{default \code{FALSE}}
}
\value{
  Returns an object of class \code{dist}.
}
\author{Witold E. Wolski \email{witek96@users.sourceforge.net}}
\keyword{methods}
