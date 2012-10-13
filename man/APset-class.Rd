\name{APset-class}
\Rdversion{1.1}
\docType{class}
\alias{APset-class}
\alias{[,APset-method}
\alias{[[,APset-method}
\alias{[[<-,APset-method}
\alias{[<-,APset-method}
\alias{ap,APset-method}
\alias{c,APset-method}
\alias{cid,APset-method}
\alias{cid<-,APset-method}
\alias{coerce,APset,AP-method}
\alias{coerce,APset,list-method}
\alias{coerce,list,APset-method}
\alias{length,APset-method}
\alias{show,APset-method}
\alias{view,APset-method}

\title{Class "APset"}
\description{
List-like container for storing the atom pair descriptors of a many compounds as objects of class \code{AP}. This container is used for structure similarity searching of compounds.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("APset", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{AP}:}{Object of class \code{"list"} }
    \item{\code{ID}:}{Object of class \code{"character"} }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "APset")}: subsetting of class with bracket operator }
    \item{[[}{\code{signature(x = "APset")}: returns single component as \code{AP} object }
    \item{[[<-}{\code{signature(x = "APset")}: replacement method for single \code{AP} component }
    \item{[<-}{\code{signature(x = "APset")}: replacement method for several \code{AP} components }
    \item{ap}{\code{signature(x = "APset")}: returns atom pair list from AP slot }
    \item{c}{\code{signature(x = "APset")}: concatenates two \code{APset} containers }
    \item{cid}{\code{signature(x = "APset")}: returns all compound identifiers from ID slot }
    \item{cid<-}{\code{signature(x = "APset")}: replacement method for compound identifiers in ID slot }
    \item{coerce}{\code{signature(from = "APset", to = "AP")}: \code{as(apset, "AP")} }
    \item{coerce}{\code{signature(from = "APset", to = "list")}: \code{as(apset, "list")}  }
    \item{coerce}{\code{signature(from = "list", to = "APset")}: \code{as(list, "APset")}  }
    \item{length}{\code{signature(x = "APset")}: returns number of entries stored in object }
    \item{show}{\code{signature(object = "APset")}: prints summary of \code{APset}  }
    \item{view}{\code{signature(x = "APset")}: prints extended summary of \code{APset} }
	 }
}
\references{
Chen X and Reynolds CH (2002). "Performance of similarity measures
in 2D fragment-based similarity searching: comparison of structural descriptors
and similarity coefficients", in J Chem Inf Comput Sci.
}
\author{
Thomas Girke
}

\seealso{
Related classes: SDF, SDFset, SDFstr, AP, FPset, FP. 

Functions: \code{SDF2apcmp}, \code{apset2descdb}, \code{cmp.search}, \code{cmp.similarity}
}
\examples{
showClass("APset")

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
