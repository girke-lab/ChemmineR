\name{atomsubset}
\alias{atomsubset}
\title{
Subset SDF/SDFset Objects by Atom Index to Obtain Substructure
}
\description{
Function to obtain a substructure from SDF/SDFset objects by providing a row index for the atom block in an SDF referencing the atoms of interest. The function subsets both the atom and bond block(s) accordingly.
}
\usage{
atomsubset(x, atomrows, type="new", datablock = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object of class \code{SDFset} or \code{SDF}
}
  \item{atomrows}{
The argument \code{atomrows} can be assigned a numeric index referencing the atoms in the atom block of \code{x}. If \code{x} is of class \code{SDF}, the index needs to be provided as \code{vector}. If \code{x} is of class \code{SDFset}, the same number of index vectors as molecules stored in \code{x} need to be passed on in a list with component names identical to the component (molecule) names stored in \code{x}.  
}
  \item{type}{
The argument \code{type="new"} assigns new atom numbers to a subsetted SDF, while \code{type="old"} maintains the numbering of the source SDF. 
}
  \item{datablock}{
By default the data block(s) in \code{SDF/SDFset} objects are removed after atom subsetting. The setting \code{datablock=TRUE} will maintain the data block information in the subsetted result.
}
}
\details{
...
}
\value{
object of class \code{SDF} or \code{SDFset}
}
\references{
...
}
\author{
Thomas Girke
}
\seealso{
...
}
\examples{
## Instance of SDFset class
data(sdfsample)
sdfset <- sdfsample

## Subset one or more molecules with atom index(es) to obtain substructure(s)
atomsubset(sdfset[[1]], atomrows=1:18) 
indexlist <- list(1:18, 1:12)
names(indexlist) <- cid(sdfset[1:2])
atomsubset(sdfset[1:2], atomrows=indexlist)
}
\keyword{ utilities }
