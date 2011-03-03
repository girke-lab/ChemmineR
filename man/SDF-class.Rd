\name{SDF-class}
\Rdversion{1.1}
\docType{class}
\alias{SDF-class}
\alias{[,SDF-method}
\alias{[[,SDF-method}
\alias{[[<-,SDF-method}
\alias{[<-,SDF-method}
\alias{atomblock,SDF-method}
\alias{atomcount,SDF-method}
\alias{bondblock,SDF-method}
\alias{coerce,character,SDF-method}
\alias{coerce,list,SDF-method}
\alias{coerce,SDF,character-method}
\alias{coerce,SDF,list-method}
\alias{coerce,SDF,SDFset-method}
\alias{coerce,SDF,SDFstr-method}
\alias{datablock,SDF-method}
\alias{datablocktag,SDF-method}
\alias{header,SDF-method}
\alias{plot,SDF-method}
\alias{sdf2list,SDF-method}
\alias{sdf2str,SDF-method}
\alias{sdfid,SDF-method}
\alias{show,SDF-method}

\title{Class "SDF"}
\description{
Container for storing every element of a single molecule defined in an SD/MOL
file without information loss in a list-like container. The import occurs via
the \code{SDFstr} container class. The header block is stored as named
character vector, the atom/bond blocks as matrices and the data block as named
character vector.  }

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SDF", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{header}:}{Object of class \code{"character"} }
    \item{\code{atomblock}:}{Object of class \code{"matrix"} }
    \item{\code{bondblock}:}{Object of class \code{"matrix"} }
    \item{\code{datablock}:}{Object of class \code{"character"} }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "SDF")}:  subsetting of class with bracket operator }
    \item{[[}{\code{signature(x = "SDF")}: returns one of the four object components }
    \item{[[<-}{\code{signature(x = "SDF")}: replacement method for the four sub-components }
    \item{[<-}{\code{signature(x = "SDF")}: replacement method for the four sub-components }
    \item{atomblock}{\code{signature(x = "SDF")}: returns atom block as matrix}
    \item{atomcount}{\code{signature(x = "SDF")}: returns atom frequency }
    \item{bondblock}{\code{signature(x = "SDF")}: returns bond block as matrix }
    \item{coerce}{\code{signature(from = "character", to = "SDF")}: \code{as(character, "SDF")} }
    \item{coerce}{\code{signature(from = "list", to = "SDF")}: \code{as(list, "SDF")} }
    \item{coerce}{\code{signature(from = "SDF", to = "character")}: \code{as(sdf, "character")} }
    \item{coerce}{\code{signature(from = "SDF", to = "list")}: \code{as(sdf, "list")} }
    \item{coerce}{\code{signature(from = "SDF", to = "SDFset")}: \code{as(sdf, "SDFset")} }
    \item{coerce}{\code{signature(from = "SDF", to = "SDFstr")}: \code{as(SDF, "SDFstr")} }
    \item{coerce}{\code{signature(from = "SDFset", to = "SDF")}: \code{as(sdfset, "SDF")} }
    \item{datablock}{\code{signature(x = "SDF")}: returns data block as named character vector }
    \item{datablocktag}{\code{signature(x = "SDF")}: returns data block as named character vector with subsetting support }
    \item{header}{\code{signature(x = "SDF")}: returns header block as named character vector }
    \item{plot}{\code{signature(x = "SDF")}: plots molecule structure for \code{SDF} object}
    \item{sdf2list}{\code{signature(x = "SDF")}: returns \code{SDF} object as \code{list}}
    \item{sdf2str}{\code{signature(sdf = "SDF")}: returns \code{SDF} object as \code{character} vector}
    \item{sdfid}{\code{signature(x = "SDF")}: returns molecule ID field from header block }
    \item{show}{\code{signature(object = "SDF")}: prints summary of \code{SDF} }
	 }
}
\references{
SDF format definition: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
}
\author{
Thomas Girke
}

\seealso{
Related classes: SDFset, SDFstr, AP, APset
}
\examples{
showClass("SDF")

## Instances of SDF class
data(sdfsample); sdfset <- sdfsample
(sdf <- sdfset[[1]]) # returns first molecule in sdfset as SDF object

## Accessing SDF components
header(sdf); atomblock(sdf); bondblock(sdf); datablock(sdf)
sdfid(sdf)

## Plot molecule structure of SDF 
plot(sdf) # plots to R graphics device
# sdf.visualize(sdf) # viewing in browser
}
\keyword{classes}
