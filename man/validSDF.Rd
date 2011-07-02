\name{validSDF}
\alias{validSDF}
\title{
Validity check of SDFset
}
\description{
Performs validity check of SDFs stored in \code{SDFset} objects. Currently, the function tests whether the atom block and the bond block in each \code{SDF} component of an \code{SDFset} have at least \code{Nabcol} and \code{Nbbcol} columns (default is 3 for both). The function returns a logical vector with \code{TRUE} values for valid compounds and \code{FALSE} values for invalid ones.
}
\usage{
validSDF(x, Nabcol = 3, Nbbcol = 3, logic = "&")
}
\arguments{
  \item{x}{
\code{x} object of class \code{SDFset}
}
  \item{Nabcol}{
minimum number of columns in atom block
}
  \item{Nbbcol}{
minimum number of columns in bond block
}
  \item{logic}{
logical connection (& or |) among Nabcol and Nbbcol cutoffs 
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
