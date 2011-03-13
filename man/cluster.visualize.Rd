\name{cluster.visualize}
\alias{cluster.visualize}
\title{visualize clustering result using multi-dimensional scaling}
\description{
	'cluster.visualize' takes clustering result returned by 'cmp.cluster' and generate multi-dimensional scaling plot for visualization purpose.
}
\usage{
cluster.visualize(db, cls, size.cutoff, distmat=NULL, color.vector=NULL, non.interactive="", cluster.result=1, dimensions=2, quiet=FALSE, highlight.compounds=NULL, highlight.color=NULL, ...)
}
\arguments{
  \item{db}{The desciptor database, in the format returned by 'cmp.parse'.}
  \item{cls}{The clustering result returned by 'cmp.cluster'.}
  \item{size.cutoff}{The cutoff size for clusters considered in this visualization. Clusters of size smaller than the cutoff will not be considered.}
  \item{distmat}{A distance matrix that corresponds to the 'db'. If not provided, it will be computed on-the-fly in an efficient manner.}
  \item{color.vector}{Colors to be used in the plot. If the number of colors in the vector is not enough for the plot, colors will be reused. If not provided, color will be generated and randomly sampled from 'rainbow'.}
  \item{non.interactive}{If provided, will enable the non-interactive mode, and the plot will be in an eps file named after this value.}
  \item{cluster.result}{Used to select the clustering result if multiple clustering results are present in 'cls'.}
  \item{dimensions}{Dimensionality to be used in visualization. See details.}
  \item{quiet}{Whether to supress the progress bar.}
  \item{highlight.compounds}{A vector of compound IDs, corresponding to compounds to be highlighted in the plot. A highlighted compound is represented as a filled circle.}
  \item{highlight.color}{Color used for highlighted compounds. If not set, a highlighted compounds will have the same color as that used for other compounds in the same cluster.}
  \item{...}{Further arguments will be passed to 'cmp.similarity' to calculate similarity matrix.}
}
\details{
	'cluster.visualize' internally calls the 'cmdscale' function to
	generate a set of points in 2-D for the compounds in selected clusters.
	Note that for compounds in clusters smaller than the cutoff size, they
	will not be considered in this calculation - their entries in 'distmat'
	will be discarded if 'distmat' is provided, and distances involving
	them will not be computed if 'distmat' is not provided.

	To determine the value for 'size.cutoff', you can use
	'cluster.sizestat' to see the size distribution of clusters.

	Because 'cmp.cluster' function allows you to perform multiple clustering
	processes simultaneously with different cutoff values, the 'cls' parameter
	may point to a data frame containing multiple clustering results. The user 
	can use 'cluster.result' to specify which result to use. By default, this
	is set to 1, and the first clustering result will be used in visualization.
	Whatever the value is, in interactive mode (described below), all clustering
	result will be displayed when a compound is selected in the interactive plot.

	If the colors provided in 'color.vector' are not enough to distinguish
	clusters by colors, the function will silently reuse the colors,
	resulting multiple clusters colored in the same color. We suggest you
	use 'cluster.sizestat' to see how many clusters will be selected using
	your 'size.cutoff', or simply provide no 'color.vector'.

	If 'non.interative' is not set, the final plot is interactive. You will
	be able to select points by clicking them. When you click on any point,
	information about the compound represented by that point will be displayed.
	This includes the cluster ID, cluster size, compound index in the SDF
	and compound name if any. You can then perform another
	selection. To exit this process, right click on X11 device or press ESC
	in non-X11 device (Quartz and Windows).

	By default, 'dimensions' is set to 2, and the built-in 'plot' function will
	be used for plotting. If you need to do 3-Dimensional plotting, set 'dimensions'
	to 3, and pass the returned value to 3D plot utilities, such as 'scatterplot3d'
	or 'rggobi'. This package does not perform 3D plot on its own.
}
\value{
	This function returns a data frame of MDS coordinates and clustering
	result. This value can be passed to 3D plot utilities such as
	'scatterplot3d' and 'rggobi'.

	The last column of the output gives whether the compounds have been
clicked in the interactive mode.
}
\author{Y. Eddie Cao}
\seealso{\code{\link{cmp.parse}}, \code{\link{cmp.cluster}}, \code{\link{cluster.sizestat}}}
\examples{
## Load sample SD file
# data(sdfsample); sdfset <- sdfsample

## Generate atom pair descriptor database for searching
# apset <- sdf2ap(sdfset) 

## Loads same atom pair sample data set provided by library
data(apset) 
db <- apset

## cluster db with 2 cutoffs
clusters <- cmp.cluster(db, cutoff=c(0.5, 0.4))

## Return size stats
sizestat <- cluster.sizestat(clusters)

## Visualize results, using a cutoff of 3, write to file 'test.eps'
coord <- cluster.visualize(db, clusters, 2, non.interactive="test.eps")

\dontrun{
## visualize it in interactive mode, using a cutoff of 3 and the 2nd clustering result
coord <- cluster.visualize(db, clusters, cluster.result=2, 3)

## 3D visualization with scatterplot3d
coord <- cluster.visualize(db, clusters, 3, dimensions=3)
library(scatterplot3d)
scatterplot3d(coord)

}
}
\keyword{utilities}
