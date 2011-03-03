\name{cmp.parse1}
\alias{cmp.parse1}
\title{Parsing an SDF file and calculate the descriptor for one compound}
\description{
	Read SDF information from an SDF file or connection, parse the first
	compound, and calculate the descriptor for that compound. The returned
	descriptor can be added to database returned by 'cmp.parse' or be used
	as the query structure when calling 'search'. This function will only
	parse one compound and return only the descriptor. To parse all
	compounds in an SDF file, use 'cmp.parse'.
}
\usage{
cmp.parse1(filename)
}
\arguments{
  \item{filename}{The file name of the SDF file or a URL or a connection.}
}
\details{
	'cmp.parse1' can take a file name or a URL or a
	connection. When a connection is used, the current line must be the
	first line of SDF of the compound to be parsed. 'cmp.parse1'  will skip
	the header and parse from the 4th line.
	Therefore, the compound ID information will be skipped. After the
	parsing is done, if 'filename' is a connection, it will then point to
	the line after the connection table of SDF. You can use some other
	procedure to parse the annotation block.
}
\value{
	Return the descriptor, which is encoded as a vector.
}
\references{Chen X and Reynolds CH (2002). "Performance of similarity measures
in 2D fragment-based similarity searching: comparison of structural descriptors
and similarity coefficients", in \emph{J Chem Inf Comput Sci}.}
\author{Y. Eddie Cao, Li-Chang Cheng}
\seealso{\code{\link{cmp.parse}}, \code{\link{cmp.search}}, \code{\link{cmp.cluster}},
    \code{\link{cmp.similarity}}}
\examples{
# load an SDF file from web and parse it
\dontrun{structure <- cmp.parse1("http://bioweb.ucr.edu/ChemMineV2/compound/Aurora/b32:NNQS2MBRHAZTI===/sdf")}
}
\keyword{utilities}
