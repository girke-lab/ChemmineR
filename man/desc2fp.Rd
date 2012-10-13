\name{desc2fp}
\alias{desc2fp}
\title{
Fingerprints from descriptor vectors
}
\description{
Generates fingerprints from descriptor vectors such as atom pairs stored in \code{APset} or \code{list} containers. The obtained fingerprints can be used for structure similarity comparisons, searching and clustering. Due to their compact size, computations on fingerprints are often more time and memory efficient than on their much more complex atom pair counterparts.
}
\usage{
desc2fp(x, descnames=1024, type = "FPset")
}
\arguments{
  \item{x}{
Object of classe \code{APset} or \code{list} of vectors
}
  \item{descnames}{
Descriptor set to consider for fingerprint encoding. If a single value from 1-4096 is provided then the function uses the corresponding number of the most frequent atom pairs stored in the \code{apfp} data set provided by the package. Alternatively, one can provide here any custom atom pair selection in form of a \code{character} vector.
}
  \item{type}{
return fingerprint set as \code{FPset}, \code{matrix} or \code{character} vector 
}
}
\details{
...
}
\value{
\code{matrix} or \code{character} vectors
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
Functions: \code{sdf2ap}, \code{SDF2apcmp}, \code{apset2descdb}, \code{cmp.search}, \code{cmp.similarity}

Related classes: SDF, SDFset, SDFstr, APset. 
}
\examples{
## Instance of SDFset class
data(sdfsample)
sdfset <- sdfsample[1:10]

## Compute atom pair library
apset <- sdf2ap(sdfset)

## Compute atom pair fingerprint matrix using internal atom pair
## selection containing 4096 most common atom pairs in DrugBank.
## For details see ?apfp. The following example uses from this 
## set the 1024 most frequent atom pairs: 
fpset <- desc2fp(x=apset, descnames=1024, type="FPset")

## Alternatively, one can provide any custom atom pair selection. Here
## 1024 most common ones in apset object.
fpset1024 <- names(rev(sort(table(unlist(as(apset, "list")))))[1:1024])
fpset2 <- desc2fp(x=apset, descnames=fpset1024, type="FPset")

## A more compact way of storing fingerprints is as character values
fpchar <- desc2fp(x=apset, descnames=1024, type="character")

## Convert character fingerprints back to FPset or matrix
fpset <- as(fpchar, "FPset")
fpma <- as.matrix(fpset)

## Similarity searching returning Tanimoto similarity coefficients
fpSim(x=fpset[1], y=fpset)

## Clustering example
simMAap <- sapply(cid(fpset), function(x) fpSim(x=fpset[x], fpset, sorted=FALSE)) 
hc <- hclust(as.dist(1-simMAap), method="single")
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE)
}
\keyword{ utilities }
