\name{write.SDF}
\alias{write.SDF}
\title{
SDF export function
}
\description{
Writes one or many molecules stored in a \code{SDFset}, \code{SDFstr} or \code{SDF} object to SD file.
}
\usage{
write.SDF(sdf, file, cid = FALSE, ...)
}
\arguments{
  \item{sdf}{
object of class \code{SDFset}, \code{SDFstr} or \code{SDF}
}
  \item{file}{
name of SD file to write to
}
  \item{cid}{
if \code{cid = TRUE} and an \code{SDFset} object is provide as input, then the compound IDs in the ID slot of the \code{SDFset} are used for compound naming
}
  \item{\dots}{
the optional arguments of the \code{sdf2str} function can be provided here, including \code{head, ab, bb, db}; details are provided in the help page for the \code{sdf2str} function
}
}
\details{
If the \code{write.SDF} function is supplied with an \code{SDFset} object, then it uses internally the \code{sdf2str} function to allow customizing the resulting SD file. For this all optional arguments of the \code{sdf2str} function can be passed on to \code{write.SDF}.
}

\references{
SDF format definition: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
}
\author{
Thomas Girke
}

\seealso{
Import function: \code{read.SDFset}, \code{read.SDFstr}
}
\examples{
## Instance of SDFset class
data(sdfsample); sdfset <- sdfsample

## Write objects of classes SDFset/SDFstr/SDF to file
# write.SDF(sdfset[1:4], file="sub.sdf")

## Example for writing customized SDFset to file containing
## ChemmineR signature, IDs from SDFset and no data block
# write.SDF(sdfset[1:4], file="sub.sdf", sig=TRUE, cid=TRUE, db=NULL)

## Example for injecting a custom matrix/data frame into the data block of an
## SDFset and then writing it to an SD file
props <- data.frame(MF=MF(sdfset), MW=MW(sdfset), atomcountMA(sdfset))
datablock(sdfset) <- props
view(sdfset[1:4])
# write.SDF(sdfset[1:4], file="sub.sdf", sig=TRUE, cid=TRUE)
}

\keyword{ utilities }
