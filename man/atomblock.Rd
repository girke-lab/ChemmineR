\name{atomblock}
\alias{atomblock}
\alias{atomblock<-}
\alias{atomblock-methods}
\title{
Return atom block
}
\description{
Returns atom block(s) from an object of class SDF or SDFset.
}
\usage{
atomblock(x)
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
\code{matrix} if \code{SDF} is provided or \code{list} of matrices if \code{SDFset} is provided
}
\references{
...
}
\author{
Thomas Girke
}

\seealso{
\code{header}, \code{atomcount}, \code{bondblock}, \code{datablock}, \code{cid}, \code{sdfid} 
}
\examples{
## SDF/SDFset instances
data(sdfsample)
sdfset <- sdfsample
sdf <- sdfset[[1]]

## Extract atome block
atomblock(sdf)
atomblock(sdfset[1:4])

## Replacement methods
sdfset[[1]][[2]][1,1] <- 999
sdfset[[1]]
atomblock(sdfset)[1:2] <- atomblock(sdfset)[3:4]
atomblock(sdfset[[1]]) == atomblock(sdfset[[3]]) 
view(sdfset[1:2])
}

\keyword{ utilities }
