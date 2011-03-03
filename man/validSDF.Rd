\name{validSDF}
\alias{validSDF}
\title{
Validity check of SDFset
}
\description{
Performs a validity check of SDFs stored in \code{SDFset} objects. Currently, the function tests whether the atom block and the bond block in each \code{SDF} component of an \code{SDFset} have at least three columns. The function returns a logical vector with \code{TRUE} values for valid compounds and \code{FALSE} values for invalid ones.
}
\usage{
validSDF(x)
}
\arguments{
  \item{x}{
\code{x} object of class \code{SDFset}
}
}
\details{
The function is important to remove invalid compounds from \code{SDFset} containers.
}
\value{
\code{logical} vector of length \code{x} with \code{TRUE} for valid compounds and \code{FALSE} for invalid compounds.
}
\references{
...
}
\author{
Thomas Girke
}
\seealso{
Functions: \code{read.SDFset}
}
\examples{
## SDFset instance
data(sdfsample)
sdfset <- sdfsample

## Detect and remove invalid SDFs in SDFset. 
valid <- validSDF(sdfset)
which(!valid) # Returns index for invalid SDFs
sdfset <- sdfset[valid] # Returns only valid SDFs.
}
\keyword{ utilities }
