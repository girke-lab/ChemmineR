\name{conMA}
\alias{conMA}
\title{
Connection Table
}
\description{
Create connection matrix from \code{SDF} and \code{SDFset} objects.
}
\usage{
conMA(x, exclude = "none")
}
\arguments{
  \item{x}{
\code{SDF} or \code{SDFset} containers
}
  \item{exclude}{
if \code{exclude="none"}, then all atoms will be considered in the resulting connection table; if \code{exclude=c("H")}; then H atoms will be excluded. Argument supports any number of atoms.
}
}
\details{
...
}
\value{
If \code{x} is of class \code{SDF}, then a single connection \code{matrix} is returned. If \code{x} is of class \code{SDFset}, then a \code{list} of matrices is returned that has the same length as \code{x}.
}
\references{
...
}
\author{
Thomas Girke
}
\seealso{
Class: \code{SDF} and \code{SDFset}
}
\examples{
## Instances of SDFset class
data(sdfsample)
sdfset <- sdfsample

## Create connection matrix for first two molecules in sdfset
conMA(sdfset[1:2], exclude=c("H"))

}
\keyword{ utilities }
