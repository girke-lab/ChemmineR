\name{sdfid}
\alias{sdfid}
\title{
Return SDF compound IDs
}
\description{
Returns the compound identifiers from the header block of \code{SDF} or \code{SDFset} objects. 
}
\usage{
sdfid(x, tag = 1)
}
\arguments{
\item{x}{
object of class \code{SDFset} or \code{SDF}
}
\item{tag}{
values from 1-4 to extract different header block fields; SDF ID is in first one (default)
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
\code{atomblock}, \code{atomcount}, \code{bondblock}, \code{datablock}, \code{header}, \code{cid} 
}
\examples{
## SDF/SDFset instances
data(sdfsample)
sdfset <- sdfsample
sdf <- sdfset[[1]]

## Extract IDs from header block
sdfid(sdf, tag=1)
sdfid(sdfset[1:4])

## Extract compound IDs from ID slot in SDFset container
cid(sdfset[1:4])

## Assigning compound IDs and keeping them unique
unique_ids <- makeUnique(sdfid(sdfset))
cid(sdfset) <- unique_ids 
cid(sdfset[1:4])

}

\keyword{ utilities }






