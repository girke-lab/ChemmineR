\name{sdfsample}
\alias{sdfsample}
\docType{data}
\title{
SD file in \code{SDFset} object
}
\description{
First 100 compounds from PubChem SD file: Compound_00650001_00675000.sdf.gz
}
\usage{data(sdfsample)}
\format{
  Object of class \code{sdfset}
}
\details{
Object stores 100 molecules from a sample SD file.
}
\source{
ftp://ftp.ncbi.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_00650001_00675000.sdf.gz
}
\references{
SDF format definition: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
}
\examples{
data(sdfsample)
sdfset <- sdfsample
view(sdfset[1:4])
}
\keyword{datasets}
