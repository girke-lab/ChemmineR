\name{AP-class}
\Rdversion{1.1}
\docType{class}
\alias{AP-class}
\alias{ap,AP-method}
\alias{show,AP-method}

\title{Class "AP"}
\description{
Container for storing the atom pair descriptors of a single compound as numeric vector. The atom pairs are used as structural similarity measures and for compound similarity searching.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AP", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{AP}:}{Object of class \code{"numeric"} }
  }
}
\section{Methods}{
  \describe{
    \item{ap}{\code{signature(x = "AP")}: returns atom pairs as numeric vector }
    \item{coerce}{\code{signature(from = "APset", to = "AP")}: \code{as(apset, "AP")} }
    \item{show}{\code{signature(object = "AP")}: prints summary of \code{AP} }
	 }
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
Related classes: SDF, SDFset, SDFstr, APset. 

Functions: \code{SDF2apcmp}, \code{apset2descdb}, \code{cmp.search}, \code{cmp.similarity}
}
\examples{
showClass("AP")

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

\keyword{classes}
