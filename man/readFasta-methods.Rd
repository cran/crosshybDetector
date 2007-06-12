\name{readFasta-methods}
\docType{methods}
\alias{seqgrep}
\alias{infogrep}
\alias{readFasta}
\alias{readFasta-methods}
\alias{readFasta,AASequence-method}
\alias{readFasta,AASequenceList-method}
\title{reader for databases in Fasta file format.}
\section{Methods}{\describe{
\item{object = "AASequenceList"}{Read sequence database in FASTA format.}
}}
\description{
  Reads file in FASTA format. An sequence entry in the fasta dabase
  consists of a header line starting with the \">\" sign which contains
  the id of the sequence. It is followed by several lines of the protein
  sequence. Because the header line may vary the user can specify a
  function which will extract (grep) the id.
\code{
  infogrep <- function(x)
  {
    return(sub("^>([a-zA-Z0-9]+) .+","\\1",x,perl=TRUE))
  }
  }
}
\arguments{
  \item{object}{either \code{AASequence} or \code{AASequenceList}}
  \item{file}{Path to the file}
  \item{grepinfo}{Function which extracts the id from the fasta
    description line.}
  \item{grepseq}{Function to pre-process the sequence.}
}
\value{
  Returns an object of class \code{AASequenceList}.
}
\seealso{\code{\link{AASequenceList-class}}}
\references{\url{http://r4proteomics.sourceforge.net}}
\author{Witold E. Wolski \email{witek96@users.sourceforge.net}}
\examples{
mySequlist <- new("AASequenceList",info="my sequence list")
#mySequlist<-readFasta(mySequlist,"ex.fasta",grepinfo=infogrep,grepseq=seqgrep)
}
\keyword{methods}

