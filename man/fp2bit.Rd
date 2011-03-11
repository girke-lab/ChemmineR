\name{fp2bit}
\alias{fp2bit}
\title{
Convert base 64 fingerprints to binary
}
\description{
The function converts the base 64 encoded PubChem fingerprints to a binary \code{matrix} or a \code{character} vector. If applied to a \code{SDFset} object, then its data block needs to contain the PubChem fingerprint information.
}
\usage{
fp2bit(x, type = 2, fptag = "PUBCHEM_CACTVS_SUBSKEYS")
}
\arguments{
  \item{x}{
Object of class \code{SDFset} or \code{matrix}
}
  \item{type}{
If set to \code{1}, the results are returned as binary \emph{matrix}. If set to \code{2} (default), the results are returned as \code{character} strings in a named \emph{vector}. 
}
  \item{fptag}{
Name tag in SDF data block where the PubChem fingerprints are stored. Default is set to "PUBCHEM_CACTVS_SUBSKEYS".
}
}
\details{
...
}
\value{
Returns results as binary \code{matrix} when \code{type=2} or as \code{character} strings stored in a named \code{vector} when \code{type=1}.
}
\references{
See PubChem fingerprint specification at: ftp://ftp.ncbi.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
}
\author{
Thomas Girke
}
\seealso{
Functions: \code{fpSim}
}
\examples{
## Load PubChem SDFset sample
data(sdfsample); sdfset <- sdfsample
cid(sdfset) <- sdfid(sdfset)

## Convert base 64 encoded fingerprints to character vector or binary matrix
fpset <- fp2bit(x=sdfset, type=1)
fpset <- fp2bit(x=sdfset, type=2)

## Pairwise compound structure comparisons
fpSim(x=fpset[1,], y=fpset[2,]) 

## Structure similarity searching: x is query and y is fingerprint database
fpSim(x=fpset[1,], y=fpset) 

## Compute fingerprint based Tanimoto similarity matrix 
simMA <- sapply(rownames(fpset), function(x) fpSim(x=fpset[x,], fpset)) 

## Hierarchical clustering with simMA as input
hc <- hclust(as.dist(simMA), method="single")

## Plot hierarchical clustering tree
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE)
}
\keyword{utilities}
