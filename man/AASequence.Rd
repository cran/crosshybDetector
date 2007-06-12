\name{AASequence}
\alias{AASequence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{AASequence constructor}
\description{
  AASequence constructor
}
\usage{
AASequence(info, sequence)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{info}{info character}
  \item{sequence}{Protein sequence character}
}
\details{
}
\value{
  Returns object of class \code{AAsequence}
}
\references{ \url{http://r4proteomics.sourceforge.net}}
\author{Witold E. Wolski \email{witek96@users.sourceforge.net}}
%\note{ ~~further notes~~ }
\seealso{\code{\link{AASequence-class}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
AASequence("my seq","MEDQVGFGFR")
## The function is currently defined as
function(info,sequence)
  {
    if(!missing(info) && !missing(sequence))
      new("AASequence",sequence,info)
    else if(!missing(sequence))
      new("AASequence",sequence)
    else
      new("AASequnece")
  }
}
\keyword{misc}% at least one, from doc/KEYWORDS
