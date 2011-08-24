\name{groups}
\alias{groups}
\title{
Enumeration of Functional Groups and Atom Neighbors
}
\description{
Returns frequency information of functional groups in molecules provided as \code{SDF} or \code{SDFset} objects. Alternatively, the function can return for each atom its atom/bond neighbor information.
}
\usage{
groups(x, groups = "fctgroup", type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{SDF} or \code{SDFset} containers
}
  \item{groups}{
if \code{groups="fctgroup"}, frequencies of functional groups are returned; if \code{groups="neighbors"}, atom/bond neighbor information is returned. 
}
  \item{type}{
if \code{type="all"}, then the complete neighbor information is generated for each atom in a molecule; if \code{type="count"}, the neighbors are enumerated in a list and if \code{type="countMA"}, then the counts of atom neighbors or functional groups are returned in a frequency matrix. 
}
}
\details{
At this point this function is in an experimental stage.
}
\value{
...
}
\references{
...
}
\author{
Thomas Girke
}
\seealso{
...
}
\examples{
## Instances of SDFset class
data(sdfsample)
sdfset <- sdfsample

## Enumerate functional groups
groups(sdfset[1:20], groups="fctgroup", type="countMA") 

## Report atom/bond neighbors
groups(sdfset[1:4], groups="neighbors", type="countMA")
groups(sdfset[1:4], groups="neighbors", type="count")
groups(sdfset[1:4], groups="neighbors", type="all")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
