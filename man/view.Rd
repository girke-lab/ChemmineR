\name{view}
\alias{view}
\alias{view-methods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Viewing of complex objects
}
\description{
Convenience function for viewing the content of complex objects like SDFset and APset containers. The function
is a shorthand wrapper for \code{as(sdfset, "SDF")} and \code{as(apset, "AP")}.
}
\usage{
view(x)
}
\arguments{
  \item{x}{ object of class \code{SDFset} or \code{APset}
}
}
\details{
...
}
\value{
List populated with \code{SDF} and \code{AP} components.
}
\references{
...
}
\author{
Thomas Girke
}

\seealso{
Classes: SDF, SDFset, AP, APset
}
\examples{
## Viewing content of SDFset 
data(sdfsample); sdfset <- sdfsample
view(sdfset[1:4])

## Viewing content of APset 
apset <- sdf2ap(sdfset[1:10])
view(apset)
}
\keyword{ utilities }
