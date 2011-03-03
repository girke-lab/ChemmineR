\name{SDFset2list}
\alias{SDFset2list}
\alias{SDFset2list-methods}
\title{
\code{SDFset} to \code{list}
}
\description{
Returns object of class \code{SDFset} as \code{list} where each component conists of a \code{list} of the four SDF sub-components: header block, atom block, bond block and data block.
}
\usage{
SDFset2list(x)
}
\arguments{
  \item{x}{
object of class \code{SDFset}
}
}
\details{
...
}
\value{
\item{list}{containing one or many lists each with following components:}
\item{character}{SDF header block}
\item{matrix}{SDF bond block}
\item{matrix}{SDF atom block}
\item{character}{SDF data block}
}
\references{
SDF format definition: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp 
}
\author{
Thomas Girke
}

\seealso{
Functions: \code{sdfstr2list}, \code{sdf2str}, \code{sdf2list}, \code{SDFset2SDF}
}
\examples{
## Instance of SDFset class
data(sdfsample); sdfset <- sdfsample
sdfset 

## Returns sdfset as list
SDFset2list(sdfset[1:4])
as(sdfset, "list")[1:4] # similar result
}

\keyword{ utilities }
