\name{SDFset-class}
\Rdversion{1.1}
\docType{class}
\alias{SDFset-class}
\alias{[,SDFset-method}
\alias{[[,SDFset-method}
\alias{[[<-,SDFset-method}
\alias{[<-,SDFset-method}
\alias{atomblock,SDFset-method}
\alias{atomcount,SDFset-method}
\alias{bondblock,SDFset-method}
\alias{c,SDFset-method}
\alias{cid,SDFset-method}
\alias{cid<-,SDFset-method}
\alias{header<-,SDFset-method}
\alias{atomblock<-,SDFset-method}
\alias{bondblock<-,SDFset-method}
\alias{datablock<-,SDFset-method}
\alias{coerce,list,SDFset-method}
\alias{coerce,SDFset,list-method}
\alias{coerce,SDFset,SDF-method}
\alias{coerce,SDFset,SDFstr-method}
\alias{datablock,SDFset-method}
\alias{datablocktag,SDFset-method}
\alias{header,SDFset-method}
\alias{length,SDFset-method}
\alias{plot,SDFset-method}
\alias{sdfid,SDFset-method}
\alias{SDFset2list,SDFset-method}
\alias{SDFset2SDF,SDFset-method}
\alias{SDFset2SDF<-,SDFset-method}
\alias{show,SDFset-method}
\alias{view,SDFset-method}

\title{Class "SDFset"}
\description{
List-like container for storing one or many objects of class \code{SDF} each containing the structure definition information of molecules provided by an SD/MOL file. The \code{SDFset} is the most important class in the ChemmmineR package for accessing and manipulating information stored in SD files.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SDFset", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{SDF}:}{Object of class \code{"list"} storing \code{SDF} components }
    \item{\code{ID}:}{Object of class \code{"character"} storing compound identifiers }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "SDFset")}: subsetting of class with bracket operator }
    \item{[[}{\code{signature(x = "SDFset")}: returns single component as \code{SDF} object }
    \item{[[<-}{\code{signature(x = "SDFset")}: replacement method for single \code{SDF} component }
    \item{[<-}{\code{signature(x = "SDFset")}: replacement method for several \code{SDF} components }
    \item{atomblock}{\code{signature(x = "SDFset")}: returns all atom blocks as list }
    \item{atomcount}{\code{signature(x = "SDFset")}: returns all atom frequencies as list }
    \item{bondblock}{\code{signature(x = "SDFset")}: returns all bond blocks as list }
    \item{c}{\code{signature(x = "SDFset")}: concatenates two \code{SDFset} containers }
    \item{cid}{\code{signature(x = "SDFset")}: returns all compound identifiers from ID slot }
    \item{header<-}{\code{signature(x = "SDFset")}: replacement method for header block }
    \item{atomblock<-}{\code{signature(x = "SDFset")}: replacement method for atom block }
    \item{bondblock<-}{\code{signature(x = "SDFset")}: replacement method for bond block }
    \item{datablock<-}{\code{signature(x = "SDFset")}: replacement method for data block }
    \item{coerce}{\code{signature(from = "list", to = "SDFset")}: \code{as(list, "SDFset")} }
    \item{coerce}{\code{signature(from = "SDF", to = "SDFset")}: \code{as(sdf, "SDFset")} }
    \item{coerce}{\code{signature(from = "SDFset", to = "list")}: \code{as(sdfset, "list")} }
    \item{coerce}{\code{signature(from = "SDFset", to = "SDF")}: \code{as(sdfset, "SDF")} }
    \item{coerce}{\code{signature(from = "SDFset", to = "SDFstr")}: \code{as(sdfset, "SDFstr")} }
    \item{coerce}{\code{signature(from = "SDFstr", to = "SDFset")}: \code{as(sdfstr, "SDFset")} }
    \item{datablock}{\code{signature(x = "SDFset")}: returns all data blocks as list }
    \item{datablocktag}{\code{signature(x = "SDFset")}: returns all data blocks as named as list with subsetting support }
    \item{header}{\code{signature(x = "SDFset")}: returns all header blocks as list }
    \item{length}{\code{signature(x = "SDFset")}: returns number of entries stored in object }
    \item{plot}{\code{signature(x = "SDFset")}: plots one or many molecule structures from \code{SDFset} object }
    \item{sdfid}{\code{signature(x = "SDFset")}: returns molecule ID field from header block }
    \item{SDFset2list}{\code{signature(x = "SDFset")}: returns \code{SDFset} object as \code{list} }
    \item{SDFset2SDF}{\code{signature(x = "SDFset")}: returns \code{SDFset} object as \code{list} with \code{SDF} components}
    \item{SDFset2SDF<-}{\code{signature(x = "SDFset")}: replacement method for \code{SDFset} component in \code{SDFset} using accessor method }
    \item{show}{\code{signature(object = "SDFset")}: prints summary of \code{SDFset} }
    \item{view}{\code{signature(x = "SDFset")}: prints extended summary of \code{SDFset} }
	 }
}
\references{
SDF format definition: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
}
\author{
Thomas Girke
}

\seealso{
Related classes: SDF, SDFstr, AP, APset

Import function: read.SDFset("some_SDF_file") 

Export function: write.SDF(sdfset, "some_file.sdf")
}
\examples{
showClass("SDFset")

## Instances of SDFset class
data(sdfsample); sdfset <- sdfsample
sdfset; view(sdfset[1:4])
sdfset[[1]]

## Import and store SD File in SDFset container
# sdfset <- read.SDFset("some_SDF_file")

## Miscellaneous accessor methods
header(sdfset[1:4])
atomblock(sdfset[1:4])
atomcount(sdfset[1:4])
bondblock(sdfset[1:4])
datablock(sdfset[1:4])

## Assigning compound IDs and keeping them unique
cid(sdfset); sdfid(sdfset)
unique_ids <- makeUnique(sdfid(sdfset))
cid(sdfset) <- unique_ids

## Convert data block to matrix
blockmatrix <- datablock2ma(datablocklist=datablock(sdfset)) # Converts data block to matrix  
numchar <- splitNumChar(blockmatrix=blockmatrix) # Splits to numeric and character matrix
numchar[[1]][1:4,]; numchar[[2]][1:4,]

## Compute atom frequency matrix, molecular weight and formula
propma <- data.frame(MF=MF(sdfset), MW=MW(sdfset), atomcountMA(sdfset))
propma[1:4, ]

## Assign matrix data to data block
datablock(sdfset) <- propma 
view(sdfset[1:4])

## String Searching in SDFset
grepSDFset("650001", sdfset, field="datablock", mode="subset") # To return index, set mode="index")

## Export SDFset to SD file
# write.SDF(sdfset[1:4], file="sub.sdf", sig=TRUE)

## Plot molecule structure of SDF 
plot(sdfset[1:4]) # plots to R graphics device
# sdf.visualize(sdfset[1:4]) # viewing in browser

}
\keyword{classes}
