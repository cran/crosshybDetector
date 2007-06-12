\name{getSeqs}
\alias{getSeqs}
\title{ Read sequences }
\description{
  Reads sequences from a two columns file.
}
\usage{
getSeqs(file, start = NULL, stop = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ the text file containing the sequences. The files MUST have
  no header row and the following columns: i)probe names and ii) probe sequences }
  \item{start}{ integer. The first element of the sequence to keep.
  If is NULL (default), all the sequence will be kept}
  \item{stop}{ integer. The last element of the sequence to keep.
  If is NULL (default), all the sequence will be kept }
}
\value{
  A vector with the probe sequences
}
\author{ Paolo Uva }
\examples{
\dontrun{
# Keep sequences from position 1 to 40
seqs <- getSeqs(file="myfile.txt",start=1, stop=40)
}
}
\keyword{ misc }
