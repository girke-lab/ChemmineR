\name{SDFstr-class}
\Rdversion{1.1}
\docType{class}
\alias{SDFstr-class}
\alias{[,SDFstr-method}
\alias{[[,SDFstr-method}
\alias{[[<-,SDFstr-method}
\alias{[<-,SDFstr-method}
\alias{coerce,character,SDFstr-method}
\alias{coerce,list,SDFstr-method}
\alias{coerce,SDFstr,list-method}
\alias{coerce,SDFstr,SDFset-method}
\alias{length,SDFstr-method}
\alias{sdfstr2list,SDFstr-method}
\alias{sdfstr2list<-,SDFstr-method}
\alias{show,SDFstr-method}

\title{Class "SDFstr"}
\description{
List-like container for storing one or many molecules from an SD (or MOL) file.
Each component of an \code{SDFstr} object stores the SD data line by line from
a single molecule in a character vector. The \code{SDFstr} class is an
intermediate container to import SD files into the more important \code{SDFset}
object or to export the data back from an \code{SDFset} container to a valid SD
file.  }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SDFstr", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{a}:}{Object of class \code{"list"} with \code{character} components}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "SDFstr")}: subsetting of class with bracket operator }
    \item{[[}{\code{signature(x = "SDFstr")}: returns single component as character vector}
    \item{[[<-}{\code{signature(x = "SDFstr")}: replacement method for single \code{SDFstr} component}
    \item{[<-}{\code{signature(x = "SDFstr")}: replacement method for several \code{SDFstr} components}
    \item{coerce}{\code{signature(from = "character", to = "SDFstr")}: \code{as(character, "SDFstr")} }
    \item{coerce}{\code{signature(from = "list", to = "SDFstr")}: \code{as(list, "SDFstr")}}
    \item{coerce}{\code{signature(from = "SDF", to = "SDFstr")}: \code{as(sdf, "SDFstr")}}
    \item{coerce}{\code{signature(from = "SDFset", to = "SDFstr")}: \code{as(sdfset, "SDFstr")}}
    \item{coerce}{\code{signature(from = "SDFstr", to = "list")}: \code{as(sdfstr, "list")}}
    \item{coerce}{\code{signature(from = "SDFstr", to = "SDFset")}: \code{as(sdfstr, "SDFset")}}
    \item{length}{\code{signature(x = "SDFstr")}: returns length of \code{SDFstr}}
    \item{sdfstr2list}{\code{signature(x = "SDFstr")}: accessor method to return \code{SDFstr} as \code{list}}
    \item{sdfstr2list<-}{\code{signature(x = "SDFstr")}: replacement method for several \code{SDFstr} components }
    \item{show}{\code{signature(object = "SDFstr")}: prints summary of \code{SDFstr}}
	 }
}
\references{
SDF format definition: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
}
\author{
Thomas Girke
}
\seealso{
Related classes: SDFset, AP, APset

Import function: read.SDFstr("some_SDF_file")
}
\examples{
showClass("SDFstr")

## Instances of SDFstr class
data(sdfsample); sdfset <- sdfsample
sdfstr <- as(sdfset, "SDFstr")
sdfstr[1:4] # print summary of container content 
sdfstr[[1]] # returns character vector

## Import: sdfstr <- read.SDFstr("some_SDF_file")
## Export: write.SDF(sdfstr, "some_file.sdf")

}
\keyword{classes}
