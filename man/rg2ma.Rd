\name{rg2ma}
\alias{rg2ma}
\title{ Convert intensities to ratios }
\description{
  This function converts dual channel intensities to log2-intensity ratios
  (M-values) and log2-intensity averages (A-values).
}
\usage{
rg2ma(r, g)
}
\arguments{
  \item{r}{ a vector of red channel intensity }
  \item{g}{ a vector of green channel intensity }
}
\value{
  A two columns matrix containing log2-intensity ratios (M-values) and
  log2-intensity averages (A-values)
}
\author{ Paolo Uva }
\seealso{ \code{\link{crosshybMAplot}} }
\examples{
# Creates random data
r <- sample(1:65535,20000)
g <- sample(1:65535,20000)
ma <- rg2ma(r,g)
}
\keyword{ misc }
