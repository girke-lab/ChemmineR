\name{read.SDFstr}
\alias{read.SDFstr}
\title{
SD file to \code{SDFstr} 
}
\description{
Imports one or many molecules from an SD/MOL file and stores it in an \code{SDFstr} container. 
}
\usage{
read.SDFstr(sdfstr)
}
\arguments{
  \item{sdfstr}{
path/name to an SD file
}
}
\details{
...
}
\value{
\item{SDFstr}{ for details see ?"SDFstr-class"}
}
\references{
SDF format defintion: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
}
\author{
Thomas Girke
}

\seealso{
Functions: \code{read.SDFset}
}
\examples{
## Write instance of SDFstr class to SD file
data(sdfsample); sdfset <- sdfsample
sdfstr <- as(sdfset, "SDFstr")
# write.SDF(sdfset[1:4], file="sub.sdf")

## Import SD file 
# read.SDFstr("sub.sdf")

## Pass on SDFstr object
sdfstr <- as(sdfset, "SDFstr")
read.SDFset(sdfstr) 
}

\keyword{ utilities }
