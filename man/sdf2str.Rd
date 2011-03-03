\name{sdf2str}
\alias{sdf2str}
\alias{sdf2str-methods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\code{SDF} to \code{SDFstr}
}
\description{
Converts \code{SDF} to \code{SDFstr}. Its main use is to facilitate the export to SD files. It contains optional arguments to generate custom SDF output.
}
\usage{
sdf2str(sdf, head, ab, bb, db, cid = NULL, sig = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdf}{
object of class \code{SDF}
}
  \item{head}{
optional \code{character} vector to supply custom header block
}
  \item{ab}{
optional \code{matrix} to supply custom atom block 
}
  \item{bb}{
optional \code{matrix} to supply custom bond block 
}
  \item{db}{
optional \code{character} vector to supply custom data block 
}
  \item{cid}{
\code{character} can be provided to inject custom compound ID into header block
}
  \item{sig}{
\code{if = TRUE} then the ChemmineR signature will be injected into the header block for tracking purposes
}
  \item{\dots}{
option to pass on additional arguments
}
}
\details{
If the export function \code{write.SDF} is supplied with an \code{SDFset} object, then \code{sdf2str} is used internally to customize the export of many molecules to a single SD file using the same optional arguments.
}
\value{
\item{sdfstr}{SDF data of one molecule collapsed to character vector}
}
\references{
SDF format definition: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp 
}
\author{
Thomas Girke
}

\seealso{
Coerce functions: \code{sdfstr2list}, \code{sdf2str}, \code{SDFset2list}, \code{SDFset2SDF}

Export function: \code{write.SDF}
}
\examples{
## Instance of SDF class
data(sdfsample); sdfset <- sdfsample
sdf <- sdfset[[1]]

## Customize SDF blocks for export to SD file
sdf2str(sdf=sdf, sig=TRUE, cid=TRUE) # uses default SDF components
sdf2str(sdf=sdf, head=letters[1:4], db=NULL) # uses custom components for header and datablock

## The same arguments can be supplied to the write.SDF function for
## batch export of custom SDFs
# write.SDF(sdfset[1:4], file="sub.sdf", sig=TRUE, cid=TRUE, db=NULL)
}

\keyword{ utilities }
