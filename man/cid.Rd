\name{cid}
\alias{cid}
\alias{cid<-}
\title{
Return compound IDs
}
\description{
Returns the compound identifiers from the ID slot of an \code{SDFset} object. 
}
\usage{
cid(x)
}
\arguments{
  \item{x}{
object of class \code{SDFset} or \code{APset}
}
}
\details{
...
}
\value{
\code{character} vector}

\references{
...
}
\author{
Thomas Girke
}
\seealso{
\code{atomblock}, \code{atomcount}, \code{bondblock}, \code{datablock}, \code{header}, \code{sdfid} 
}
\examples{
## SDFset/APset instances
data(sdfsample)
sdfset <- sdfsample
apset <- sdf2ap(sdfset[1:4])

## Extract compound IDs from SDFset/APset
cid(sdfset[1:4])
cid(apset[1:4])

## Extract IDs defined in SD file
sdfid(sdfset[1:4])

## Assigning compound IDs and keeping them unique
unique_ids <- makeUnique(sdfid(sdfset))
cid(sdfset) <- unique_ids 
cid(sdfset[1:4])

## Replacement Method
cid(sdfset) <- as.character(1:100)

}

\keyword{ utilities }



