\name{atomcount}
\alias{atomcount}
\alias{atomcountMA}
\alias{MW}
\alias{MF}
\title{
Molecular property functions
}
\description{
Functions to compute molecular properties: weight, formula, atom frequencies, etc. 
}
\usage{
atomcount(x, addH = FALSE, ...)

atomcountMA(x, ...)

MW(x, mw=atomprop, ...)

MF(x, ...)
}
\arguments{
\item{x}{
object of class \code{SDFset} or \code{SDF}
}
\item{mw}{
\code{data.frame} with atomic weights; imported by default with data(atomprop); supports custom data sets
}
\item{addH}{
'addH = TRUE' should be passed on to any of these function to add hydrogens that are often not specified in SD files
}
\item{\dots}{
Arguments to be passed to/from other methods.
}
}
\details{
...
}
\value{
\item{named vector}{\code{MW} and \code{MF}}
\item{list}{\code{atomcount}}
\item{matrix}{\code{atomcountMA}}
}
\references{
Standard atomic weights (2005) from: http://iupac.org/publications/pac/78/11/2051/
}
\author{
Thomas Girke
}

\seealso{
Functions: \code{datablock}, \code{datablocktag}
}
\examples{

## Instance of SDFset class
data(sdfsample)
sdfset <- sdfsample

## Compute properties; to consider missing hydrogens, set 'addH = TRUE'
MW(sdfset[1:4], addH = FALSE)
MF(sdfset[1:4], addH = FALSE)
atomcount(sdfset[1:4], addH = FALSE)
propma <- atomcountMA(sdfset[1:4], addH = FALSE)
boxplot(propma, main="Atom Frequency")

## Example for injecting a custom matrix/data frame into the data block of an
## SDFset and then writing it to an SD file
props <- data.frame(MF=MF(sdfset), MW=MW(sdfset), atomcountMA(sdfset))
datablock(sdfset) <- props
view(sdfset[1:4])
# write.SDF(sdfset[1:4], file="sub.sdf", sig=TRUE, cid=TRUE)

}

\keyword{ utilities }
