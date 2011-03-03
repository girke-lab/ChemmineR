\name{header}
\alias{header}
\alias{header<-}
\alias{header-methods}
\title{
Return header block
}
\description{
Returns header block(s) from an object of class SDF or SDFset.
}
\usage{
header(x)
}
\arguments{
  \item{x}{
object of class \code{SDF} or \code{SDFset}
}
}
\details{
...
}
\value{
named \code{character} vector if \code{SDF} is provided or \code{list} of named \code{character} vectors if \code{SDFset} is provided
}
\references{
...
}
\author{
Thomas Girke
}
\seealso{
\code{atomblock}, \code{atomcount}, \code{bondblock}, \code{datablock}, \code{cid}, \code{sdfid} 
}
\examples{
## SDF/SDFset instances
data(sdfsample)
sdfset <- sdfsample
sdf <- sdfset[[1]]

## Extract header block
header(sdf)
header(sdfset[1:4])

## Replacement methods
sdfset[[1]][[1]][1] <- "test"
sdfset[[1]]
header(sdfset)[1] <- header(sdfset[2])  
view(sdfset[1:2])
}

\keyword{ utilities }
