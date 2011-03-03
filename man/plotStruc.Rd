\name{plotStruc}
\alias{plotStruc}
\alias{plot}
\title{
Plot compound structures
}
\description{
Plots compound structure(s) for molecules stored in SDF and SDFset containers.
}
\usage{
## Convenience plot method
# plot(x, griddim, print_cid=cid(x), print=TRUE, ...)

## Less important for user 
plotStruc(sdf, atomcex = 1.2, atomnum = FALSE, no_print_atoms = c("C"),
          noHbonds = TRUE, bondspacer = 0.08, ...)
}

\arguments{
  \item{sdf}{
Object of class \code{SDF} 
}
  \item{atomcex}{
Font size for atom labels
}
  \item{atomnum}{
If \code{TRUE}, then the atom numbers are included in the plot. They are the position numbers of each atom in the atom block of an SDF.
}
  \item{no_print_atoms}{
Excludes specified atoms from being plotted.
}
  \item{noHbonds}{
If \code{TRUE}, then the C-hydrogens and their bonds - explicitly defined in an SDF - are excluded from the plot.
}
  \item{bondspacer}{
Numeric value specifying the plotting distance for double/triple bonds.
}
\item{\dots}{
Arguments to be passed to/from other methods.
}
}
\details{
The function \code{plotStruc} depicts a single 2D compound structure based
on the XY-coordinates specified in the atom block of an SDF. The generic method
\code{plot} can be used as a convenient shorthand to plot one or many
structures at once. Both functions depend on the availability of the
XY-coordinates in the source SD file and only 2D (not 3D) representations are plotted
correctly. 

Additional arguments that can only be passed on to the \code{plot} function when supplied with
an SDFset object:

   \code{griddim}: numeric vector of length two to define the dimensions for arranging several structures in one plot. 
   
   \code{print_cid}: character vector for printing custom compound labels. Default is \code{print_cid=cid(sdfset)}.

   \code{print}: if \code{print=TRUE}, then a summary of the SDF content for each supplied compound is printed to the screen.
          This behavior is turned off with \code{print=TRUE}.
} 
\value{
Prints summary of SDF/SDFset to screen and plots their structures to graphics device.
}
\references{
...
}
\author{
Thomas Girke
}
\note{
The compound depictions created by this function are not as pretty as the structure representations generated with the \code{sdf.visualize} function. This will be improved in the future. 
}

\seealso{
\code{sdf.visualize}
}
\examples{
## Import SDFset sample set
data(sdfsample)
(sdfset <- sdfsample)

## Plot single compound structure
plotStruc(sdfset[[1]])

## Plot several compounds structures
plot(sdfset[1:4])

## Customize plot 
plot(sdfset[1:4], griddim=c(2,2), print_cid=letters[1:4], print=FALSE, noHbonds=FALSE)
}

\keyword{ aplot  }
\keyword{ utilities }
