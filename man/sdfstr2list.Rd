\name{sdfstr2list}
\alias{sdfstr2list}
\alias{sdfstr2list<-}
\alias{sdfstr2list-methods}
\title{
\code{SDFstr} to \code{list}
}
\description{
Returns objects of class \code{SDFstr} as \code{list}.
}
\usage{
sdfstr2list(x)
}
\arguments{
  \item{x}{
object of class \code{SDFstr}
}
}
\details{
...
}
\value{
\item{list}{with many of the following components:}
\item{character}{SDF content of one molecule vectorized line by line}
}
\references{
SDF format definition: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
}
\author{
Thomas Girke
}

\seealso{
Functions: \code{sdf2list}, \code{sdf2str}, \code{SDFset2list}, \code{SDFset2SDF}
}
\examples{
## Instance of SDFstr class
data(sdfsample); sdfset <- sdfsample
sdfstr <- as(sdfset, "SDFstr")

## Return as list
sdfstr2list(sdfstr)
as(sdfstr, "list") # similar result
}
\keyword{ utilities }
