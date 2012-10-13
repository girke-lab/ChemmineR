\name{FPset-class}
\Rdversion{1.1}
\docType{class}
\alias{FPset-class}
\alias{[,FPset-method}
\alias{[[,FPset-method}
\alias{[<-,FPset-method}
\alias{as.character,FPset-method}
\alias{as.matrix,FPset-method}
\alias{c,FPset-method}
\alias{cid,FPset-method}
\alias{cid<-,FPset-method}
\alias{coerce,FPset,FP-method}
\alias{coerce,matrix,FPset-method}
\alias{coerce,character,FPset-method}
\alias{length,FPset-method}
\alias{show,FPset-method}
\alias{view,FPset-method}

\title{Class \code{"FPset"}}
\description{
Container for storing fingerprints of a many compounds. This container is used for structure similarity searching of compounds.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("FPset", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{fpma}:}{Object of class \code{"matrix"} with compound identifiers stored in row names}
    
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "FPset")}: subsetting of class with bracket operator }
    \item{[[}{\code{signature(x = "FPset")}: returns single component as \code{FP} object }
    \item{[<-}{\code{signature(x = "FPset")}: replacement method for several components }
    \item{as.character}{\code{signature(x = "FPset")}: returns content as named character vector }
    \item{as.matrix}{\code{signature(x = "FPset")}: returns content as numeric matrix }
    \item{c}{\code{signature(x = "FPset")}: concatenates two \code{FPset} containers }
    \item{cid}{\code{signature(x = "FPset")}: returns all compound identifiers from row names }
    \item{cid<-}{\code{signature(x = "FPset")}: replacement method for compound identifiers  }
    \item{coerce}{\code{signature(from = "FPset", to = "FP")}: \code{as(fpset, "FP")} }
    \item{coerce}{\code{signature(from = "matrix", to = "FPset")}: \code{as(fpma, "FPset")} }
    \item{coerce}{\code{signature(from = "character", to = "FPset")}: \code{as(fpchar, "FPset")} }
    \item{length}{\code{signature(x = "FPset")}: returns number of entries stored in object }
    \item{show}{\code{signature(object = "FPset")}: prints summary of \code{FPset} }
    \item{view}{\code{signature(x = "FPset")}: prints extended summary of \code{FPset} }
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
Related classes: SDF, SDFset, SDFstr, AP, APset, FP. 
}
\examples{
showClass("FPset")

## Instance of FPset class
data(apset)
(fpset <- desc2fp(apset))

## Class usage 
fpset[1:4] # behaves like a list
fpset[[1]] # returns FP object
length(fpset) # number of compounds
cid(fpset) # returns compound ids
fpset[1] <- 0 # replacement
cid(fpset) <- 1:length(fpset) # replaces compound ids
c(fpset[1:4], fpset[11:14]) # concatenation

## Coerce FPset from/to other objects
fpma <- as.matrix(fpset) # coerces to matrix
fpchar <- as.character(fpset) # coerces to character strings
as(fpma, "FPset")
as(fpchar, "FPset")

## Compound similarity searching with FPset 
fpSim(x=fpset[1], y=fpset, method="Tanimoto", cutoff=0.4, top=4) 
}
\keyword{classes}
