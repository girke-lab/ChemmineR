\name{SDF2apcmp}
\alias{SDF2apcmp}
\title{
\code{SDF} to \code{list} for AP generation
}
\description{
Returns \code{SDF} class as \code{list} containing the components for generating atom pair descriptors.
}

\usage{
SDF2apcmp(SDF)
}

\arguments{
  \item{SDF}{\code{SDF}}
}

\details{
...
}

\value{
\item{list}{ with atom and bond components}
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
Functions: \code{sdf2ap}, \code{apset2descdb}, \code{cmp.search}, \code{cmp.similarity}
}

\examples{

## Instances of SDFset class
data(sdfsample)
sdf <- sdfsample[[1]]

## Return list 
SDF2apcmp(sdf)

}

\keyword{ utilities }
