\name{datablock}
\alias{datablock}
\alias{datablock<-}
\alias{datablocktag}
\alias{datablock-methods}
\title{
Return data block
}
\description{
Returns data block(s) from an object of class SDF or SDFset.
}
\usage{
datablock(x)

datablocktag(x, tag)
}
\arguments{
\item{x}{
object of class \code{SDF} or \code{SDFset}
}
\item{tag}{
\code{numeric} position (index) or \code{character} name of entry in data block vector
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
\code{atomblock}, \code{atomcount}, \code{bondblock}, \code{header}, \code{cid}, \code{sdfid} 
}
\examples{
## SDF/SDFset instances
data(sdfsample)
sdfset <- sdfsample
sdf <- sdfset[[1]]

## Extract data block
datablock(sdf)
datablock(sdfset[1:4])
datablocktag(sdfset, tag="PUBCHEM_OPENEYE_CAN_SMILES")

## Replacement methods
sdfset[[1]][[1]][1] <- "test"
sdfset[[1]]
datablock(sdfset)[1] <- datablock(sdfset[2])  
view(sdfset[1:2])

## Example for injecting a custom matrix/data frame into the data block of an
## SDFset and then writing it to an SD file
props <- data.frame(MF=MF(sdfset), MW=MW(sdfset), atomcountMA(sdfset))
datablock(sdfset) <- props
view(sdfset[1:4])
# write.SDF(sdfset[1:4], file="sub.sdf", sig=TRUE, cid=TRUE)

}

\keyword{ utilities }




