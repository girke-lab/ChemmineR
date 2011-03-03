\name{makeUnique}
\alias{makeUnique}
\title{
Uniquify CMP names 
}
\description{
Creates unique CMP names by appending a counter to each duplicatation set. The function can be used for any character vector.
}
\usage{
makeUnique(x, silent = FALSE)
}
\arguments{
  \item{x}{
\code{character} vector
}
  \item{silent}{
\code{silent = TRUE} suppresses message about duplicate count
}
}
\details{
The function is important to maintain unique compound names in the ID slot of \code{SDFset} containers.
}
\value{
\item{character}{ of same length as x but without duplications}
}
\references{
...
}
\author{
Thomas Girke
}

\seealso{
Functions: \code{cid}, \code{sdfid}
}
\examples{
## SDFset instance
data(sdfsample)
sdfset <- sdfsample

## Create unique compound IDs 
unique_ids <- makeUnique(sdfid(sdfset))
cid(sdfset) <- unique_ids 
cid(sdfset[1:4])

}

\keyword{ utilities }
