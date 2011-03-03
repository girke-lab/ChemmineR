\name{SDFset2SDF}
\alias{SDFset2SDF}
\alias{SDFset2SDF<-}
\alias{SDFset2SDF-methods}
\title{
\code{SDFset} to list with many \code{SDF}
}
\description{
Returns object of class \code{SDFset} as \code{list} were each component consists of an \code{SDF} object.
}
\usage{
SDFset2SDF(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object of class \code{SDFset}
}
}
\details{
...
}
\value{
\item{list}{containing one or many \code{SDF} objects}
}
\references{
SDF format definition: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp 
}
\author{
Thomas Girke
}

\seealso{
Functions: \code{sdfstr2list}, \code{sdf2str}, \code{sdf2list}, \code{SDFset2list}
}
\examples{
## Instance of SDFset class
data(sdfsample); sdfset <- sdfsample
sdfset 

## Returns sdfset as list
SDFset2SDF(sdfset[1:4])
as(sdfset, "SDF")[1:4] # similar result
view(sdfset[1:4]) # same result

}

\keyword{ utilities }
