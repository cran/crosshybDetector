\name{list2ascii}
\alias{list2ascii}
\title{ Write R list }
\description{
  This function writes a R list to an ASCII file.
}
\usage{
list2ascii(x, file)
}
\arguments{
  \item{x}{ a list }
  \item{file}{ the output file }
}
\author{ Paolo Uva }
\examples{
\dontrun{
# Create list
mylist <- list(a=1:10, b=11:20)

# Write list to file
list2ascii(mylist, file="myfile.txt")
}
}
\keyword{ misc }
