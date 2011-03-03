\name{read.SDFset}
\alias{read.SDFset}
\title{
SD file to \code{SDFset} 
}
\description{
Imports one or many molecules from an SD/MOL file and stores it in an \code{SDFset} container. 
}
\usage{
read.SDFset(sdfstr = sdfstr, ...)
}
\arguments{
  \item{sdfstr}{
path/name to an SD file; alternatively an \code{SDFstr} object can be provided
}
  \item{\dots}{
option to pass on additional arguments
}
}
\details{
...
}
\value{
\item{SDFset}{ for details see ?"SDFset-class"}
}
\references{
SDF format defintion: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
}
\author{
Thomas Girke
}

\seealso{
Functions: \code{read.SDFstr}
}
\examples{
## Write instance of SDFset class to SD file
data(sdfsample); sdfset <- sdfsample
# write.SDF(sdfset[1:4], file="sub.sdf")

## Import SD file 
# read.SDFset("sub.sdf")

## Pass on SDFstr object
sdfstr <- as(sdfset, "SDFstr")
read.SDFset(sdfstr) 
}
\keyword{ utilities }
