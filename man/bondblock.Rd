\name{bondblock}
\alias{bondblock}
\alias{bondblock<-}
\alias{bondblock-methods}
\title{
Return bond block
}
\description{
Returns bond block(s) from an object of class SDF or SDFset.
}
\usage{
bondblock(x)
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
\code{header}, \code{atomcount}, \code{atomblock}, \code{datablock}, \code{cid}, \code{sdfid} 
}
\examples{
## SDF/SDFset instances
data(sdfsample)
sdfset <- sdfsample
sdf <- sdfset[[1]]

## Extract bond block
bondblock(sdf)
bondblock(sdfset[1:4])

## Replacement methods
sdfset[[1]][[3]][1,1] <- 999
sdfset[[1]]
bondblock(sdfset)[1:2] <- bondblock(sdfset)[3:4]
bondblock(sdfset[[1]]) == bondblock(sdfset[[3]]) 
view(sdfset[1:2])
}

\keyword{ utilities }

