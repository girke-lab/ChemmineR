\name{datablock2ma}
\alias{datablock2ma}
\alias{splitNumChar}
\title{
SDF data blocks to matrix
}
\description{
Convert data blocks in \code{SDFset} to character matrix with
\code{datablock2ma}, then store its numeric columns as numeric matrix and its
character columns as character matrix.  }
\usage{
datablock2ma(datablocklist = datablock(sdfset), cleanup = " \\\(.*", ...)   

splitNumChar(blockmatrix = blockmatrix)

}
\arguments{
  \item{datablocklist}{
\code{list} of data block vectors; can be created with \code{datablock(sdfset)}
}
  \item{blockmatrix}{
\code{matrix} returned by \code{datablock2ma} 
}
  \item{cleanup}{
\code{character} pattern to be used to clean up the name fields of the data block vectors; the exact pattern matches are replaced by nothing (deleted).
}
  \item{\dots}{
option to pass on additional arguments
}
}
\details{
...
}
\value{
\item{datablock2ma}{\code{character matrix}}
\item{splitNumChar}{\code{list} with two components, a numeric matrix and a character matrix}
}
\references{
...
}
\author{
Thomas Girke
}

\seealso{
Classes: \code{SDFset}
}
\examples{

## SDFset instance
data(sdfsample)
sdfset <- sdfsample

# Convert data block to matrix  
blockmatrix <- datablock2ma(datablocklist=datablock(sdfset)) 
blockmatrix[1:4, 1:4]

# Split matrix to numeric matrix and character matrix
numchar <- splitNumChar(blockmatrix=blockmatrix)
names(numchar)
numchar[[1]][1:4,] 
numchar[[2]][1:4,] 
}

\keyword{ utilities }
