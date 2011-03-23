\name{conMA}
\alias{conMA}
\title{
Bond Matrices
}
\description{
Creates a bond matrix from \code{SDF} and \code{SDFset} objects. The matrix contains the atom labels in the row and column titles and the bond types are given in the data part as follows: 0 is no connection, 1 is a single bond, 2 is a double bond and 3 is a triple bond.
}
\usage{
conMA(x, exclude = "none")
}
\arguments{
  \item{x}{
\code{SDF} or \code{SDFset} containers
}
  \item{exclude}{
if \code{exclude="none"}, then all atoms will be considered in the resulting connection table; if \code{exclude=c("H")}, then the H atoms will be excluded. Any number of atom labels to be excluded can be passed on to this argument in form of a \code{character} vector.
}
}
\details{
...
}
\value{
If \code{x} is of class \code{SDF}, then a single bond \code{matrix} is returned. If \code{x} is of class \code{SDFset}, then a \code{list} of matrices is returned that has the same length as \code{x}.
}
\references{
...
}
\author{
Thomas Girke
}
\seealso{
Class: \code{SDF} and \code{SDFset}
}
\examples{
## Instances of SDFset class
data(sdfsample)
sdfset <- sdfsample

## Create bond matrix for first two molecules in sdfset
conMA(sdfset[1:2], exclude=c("H"))

## Return bond matrix for first molecule and plot its structure with atom numbering
conMA(sdfset[[1]], exclude=c("H"))
plot(sdfset[1], atomnum = TRUE, noHbonds=FALSE , no_print_atoms = "", atomcex=0.8)

## Return number of non-H bonds for each atom 
rowSums(conMA(sdfset[[1]], exclude=c("H")))
}
\keyword{ utilities }
