\name{rings}
\alias{rings}
\title{
Ring and Aromaticity Perception
}
\description{
Identifies all possible rings in molecules using the exhaustive ring perception
algorithm from Hanser et al (1996). In addition, the function can return all
smallest possible rings as well as aromaticity information for each ring.
}
\usage{
rings(x, upper = Inf, type = "all", arom = FALSE, inner = FALSE)
}
\arguments{
  \item{x}{
\code{SDF} or \code{SDFset} containers
}
  \item{upper}{
allows to specify an upper length limit for ring predictions. The default setting \code{upper=Inf} will return all possible rings. Smaller length limits will reduce the search space resulting in shortened compute times. 
}
  \item{type}{
if \code{type="all"}, the function returns each ring of a compound as character vector of atom symbols that are numbered by their position in the atom block of an \code{SDF/SDFset} object. Note: the example below shows how to plot structures with the same numbering information for visual inspection. If \code{type="arom"}, only aromatic rings are returned, while \code{type="count"} returns the ring and/or aromaticity counts for each compound in a matrix. 
}
  \item{arom}{
if \code{arom="TRUE"}, ring aromaticity information will be computed. If \code{type="all"}, the output is a logical vector where 'TRUE' values indicate aromatic rings in the associated ring list. If \code{type="arom"}, then the function returns only aromatic rings. A ring is considered aromatic if it meets the following requirements: (i) all atoms in the ring need to be sp2 hybridized. This means each atom has to have a double bond or at least one lone electron pair and it needs to be attached to an sp2 hybridized atom. (ii) In addition, Hueckel's rule '4n + 2' needs to be true, where 'n' is either zero or any positive integer.   
}
  \item{inner}{
if \code{inner="TRUE"}, only inner (smallest possible) rings will be returned. They are identified by first computing all possible rings and then selecting only the inner rings. Note: this requires the setting \code{upper=Inf}. If only rings below a certain size limit (e.g. 6) are of interest, then it will be more time efficient to set this limmit under the \code{upper} argument than identifying all smallest rings.  
}
}
\details{
...
}
\value{
The settings \code{type="all"} and \code{type="arom"} return \code{lists}, and \code{type="count"} returns a \code{matrix}. 
}
\references{
Hanser, Jauffret and Kaufmann (1996) A New Algorithm for Exhaustive Ring Perception in a Molecular Graph.
Journal of Chemical Information and Computer Sciences, 36: 1146-1152.
URL: http://pubs.acs.org/doi/abs/10.1021/ci960322f
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

## Return all possible rings for a single compound 
rings(sdfset[1], upper=Inf, type="all", arom=FALSE, inner=FALSE)
plot(sdfset[1], print=FALSE, atomnum=TRUE, no_print_atoms="H") 

## Return all possible rings for several compounds plus their 
## aromaticity information
rings(sdfset[1:4], upper=Inf, type="all", arom=TRUE, inner=FALSE)

## Return rings with no more than 6 atoms
rings(sdfset[1:4], upper=6, type="all", arom=TRUE, inner=FALSE)

## Return rings with no more than 6 atoms that are also armomatic
rings(sdfset[1:4], upper=6, type="arom", arom=TRUE, inner=FALSE)

## Return shortest possible rings (no complex rings)
rings(sdfset[1:4], upper=Inf, type="all", arom=TRUE, inner=TRUE)

## Count shortest possible rings
rings(sdfset[1:4], upper=Inf, type="count", arom=TRUE, inner=TRUE)
}
\keyword{ utilities }
