\name{ap}
\alias{ap}
\alias{ap-methods}
\title{
Return atom pair component of \code{AP/APset}
}
\description{
Returns atom pair component of objects of class \code{AP} or \code{APset} as list of vectors.
}
\usage{
ap(x)
}
\arguments{
  \item{x}{
Object of class \code{AP} and \code{APset}
}
}
\details{
...
}
\value{
\item{List}{ with one to many of following components:}
\item{numeric}{ atom pairs }
}
\references{
Chen X and Reynolds CH (2002). "Performance of similarity measures
in 2D fragment-based similarity searching: comparison of structural descriptors
and similarity coefficients", J Chem Inf Comput Sci.
}
\author{
Thomas Girke
}
\seealso{
Functions: \code{SDF2apcmp}, \code{apset2descdb}, \code{cmp.search}, \code{cmp.similarity}
}
\examples{

## Instance of SDFset class
data(sdfsample)
sdfset <- sdfsample[1:50]
sdf <- sdfset[[1]]

## Compute atom pair library
ap <- sdf2ap(sdf)
(apset <- sdf2ap(sdfset))
view(apset[1:4])

## Return main components of APset object
cid(apset[1:4]) # compound IDs
ap(apset[1:4]) # atom pair descriptors

## Return atom pairs in human readable format
db.explain(apset[1]) 

}

\keyword{ utilities }
