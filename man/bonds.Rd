\name{bonds}
\alias{bonds}
\title{
Bonds, charges and missing hydrogens 
}
\description{
Returns information about bonds, charges and missing hydrogens in \code{SDF} and \code{SDFset} objects. 
}
\usage{
bonds(x, type = "bonds")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{SDF} or \code{SDFset} containers
}
  \item{type}{
If \code{type="bonds"} (default), a \code{data.frame} is returned with columns: \code{atom} (atom labels), \code{Nbondcount} (observed bond count), \code{Nbondrule} (bond count according to position in periodic table) and \code{charge} (charge of each atom). 

If \code{type="charge"}, all charged atoms are returned and if \code{type="addNH"}, the number of missing hydrogens are returned for each molecule. 
}
}
\details{
...
}
\value{
If \code{x} is of class \code{SDF}, then a single \code{data.frame} or \code{vector} is returned. If \code{x} is of class \code{SDFset}, then a \code{list} of \code{data.frames} or \code{vecotors} is returned that has the same length and order as \code{x}. 
}
\references{
...
}
\author{
Thomas Girke
}
\seealso{
Functions: \code{conMA} 

Class: \code{SDF} and \code{SDFset}
}
\examples{
## Instances of SDFset class
data(sdfsample)
sdfset <- sdfsample

## Returns data frames with bonds and charges 
bonds(sdfset[1:2], type="bonds")

## Returns charged atoms in each molecule
bonds(sdfset[1:2], type="charge")

## Returns the number of missing hydrogens in each molecule
bonds(sdfset[1:2], type="addNH")

}
\keyword{ utilities }
