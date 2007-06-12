\name{computeTm}
\alias{computeTm}
\title{Computing of melting temperature}
\description{
  Computes the melting temperature (Tm) using the basic formula:

      Tm= (wA+xT) * 2 + (yG+zC) * 4

      where w,x,y,z are the number of the bases A,T,G,C in the sequence, respectively (from Marmur,J., and Doty,P. (1962) J Mol Biol 5:109-118).
}
\usage{
computeTm(seqs, plot = FALSE)
}
\arguments{
  \item{seqs}{a vector of nucleotide sequences }
  \item{plot}{if TRUE, plot the histogram of Tm }
}
\value{
  A vector with the same lenght of 'seqs' with the computed Tm
}
\references{ Marmur,J., and Doty,P. (1962) J Mol Biol 5:109-118 }
\author{ Paolo Uva }
\examples{
# Create a vector of 1000 random sequences of 30bp
code <- c("a", "c", "g", "t")
sequences <- vector("character", 1000)
for (i in 1:1000){
  sequences[i] <- paste(sample(code,30,replace=TRUE), collapse="")
}

# Compute the Tm and plot the corresponding histogram
tm <- computeTm(sequences, plot=TRUE)
table(tm)
}
\keyword{ misc }
