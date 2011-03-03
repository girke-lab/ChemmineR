\name{grepSDFset}
\alias{grepSDFset}
\title{
String search in \code{SDFset}
}
\description{
Convenience grep function for string searching in \code{SDFset} containers.
}
\usage{
grepSDFset(pattern, x, field = "datablock", mode = "subset", ignore.case = TRUE, ...)
}
\arguments{
  \item{pattern}{
search pattern
}
  \item{x}{
\code{SDFset}
}
  \item{field}{ delimits search to specific section in SDF; can be
   \code{header}, \code{atomblock}, \code{bondblock} or \code{datablock} 
}
  \item{mode}{
if \code{mode = "index"}, then the match positions are returned as vector; if \code{mode = "subset"}, a 
\code{list} with \code{SDF} components is returned where every entry has at least one query match
}
  \item{ignore.case}{
\code{TRUE} turns off case sensitivity 
}
  \item{\dots}{
option to pass on additional arguments
}
}
\details{
...
}
\value{
\item{numeric}{ index vector where the name field contains the component positions in the \code{SDFset} and the values the row positions in each sub-component. }
\item{list}{ if \code{mode = "subset"}}
}
\references{
...
}
\author{
Thomas Girke
}

\seealso{
Class: \code{SDFset}
}
\examples{

## Instances of SDFset class
data(sdfsample)
sdfset <- sdfsample

## String Searching in SDFset
q <- grepSDFset("65000", sdfset, field="datablock", mode="subset") 
as(q, "SDFset")
grepSDFset("65000", sdfset, field="datablock", mode="index") 

}

\keyword{ utilities }
