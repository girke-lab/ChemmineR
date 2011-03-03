\name{sdf2list}
\alias{sdf2list}
\title{
\code{SDF} to \code{list}
}
\description{
Returns objects of class \code{SDF} as \code{list}.
}
\usage{
sdf2list(x)
}
\arguments{
  \item{x}{
object of class \code{SDF}
}
}
\details{
...
}
\value{
\item{list}{with following components:}
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
Functions: \code{sdfstr2list}, \code{sdf2str}, \code{SDFset2list}, \code{SDFset2SDF}
}
\examples{
## Instance of SDF class
data(sdfsample); sdfset <- sdfsample
sdf <- sdfset[[1]]

## Return as list
sdf2list(sdf)
as(sdf, "list") # similar result
}
\keyword{ utilities }
