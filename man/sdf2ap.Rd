\name{sdf2ap}
\alias{sdf2ap}
\title{
Atom pair library
}
\description{
Creates from a \code{SDFset} a searchable atom pair library that is stored in a container of class \code{APset}.
}
\usage{
sdf2ap(sdfset)
}
\arguments{
  \item{sdfset}{
Objects of classes \code{SDFset} or \code{SDF}
}
}
\details{
...
}
\value{
\item{APset}{ if input is \code{SDFset}}
\item{AP}{ if input is \code{SDF}}
}
\references{
Chen X and Reynolds CH (2002). "Performance of similarity measures
in 2D fragment-based similarity searching: comparison of structural descriptors
and similarity coefficients", J Chem Inf Comput Sci.
}
\author{
Thomas Girke
}

\seealso{
Functions: \code{SDF2apcmp}, \code{apset2descdb}, \code{cmp.search}, \code{cmp.similarity}

Related classes: SDF, SDFset, SDFstr, APset. 
}
\examples{

## Instance of SDFset class
data(sdfsample)
sdfset <- sdfsample[1:50]
sdf <- sdfsample[[1]]

## Compute atom pair library
ap <- sdf2ap(sdf)
(apset <- sdf2ap(sdfset))
view(apset[1:4])

## Return main components of APset object
cid(apset[1:4]) # compound IDs
ap(apset[1:4]) # atom pair descriptors

## Return atom pairs in human readable format
db.explain(apset[1]) 

## Coerce APset to other objects 
apset2descdb(apset) # returns old list-style AP database
tmp <- as(apset, "list") # returns list
as(tmp, "APset") # converst list back to APset

## Compound similarity searching with APset
cmp.search(apset, apset[1], type=3, cutoff=0.2) 
plot(sdfset[names(cmp.search(apset, apset[6], type=2, cutoff=0.4))])

## Identify compounds with identical AP sets 
cmp.duplicated(apset, type=2)

## Structure similarity clustering 
cmp.cluster(db=apset, cutoff = c(0.65, 0.5))[1:20,]

}

\keyword{ utilities }
