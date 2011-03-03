\name{apset}
\alias{apset}
\docType{data}
\title{
Atom pairs stored in \code{APset} object
}
\description{
Atom pairs for 100 molecules stored in \code{sdfsample}.
}
\usage{data(apset)}
\format{
  Object of class \code{apset}
}
\details{
Object stores atom pairs of 100 molecules.
}
\source{
\code{apset <- sdf2ap(sdfsample)}
}
\references{
Chen X and Reynolds CH (2002). "Performance of similarity measures
in 2D fragment-based similarity searching: comparison of structural descriptors
and similarity coefficients", J Chem Inf Comput Sci.
}
\examples{
data(apset)
apset[1:4]
view(apset[1:4])
}
\keyword{datasets}
