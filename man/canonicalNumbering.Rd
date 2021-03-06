\name{canonicalNumbering}
\alias{canonicalNumbering}
\title{
	Canonical Numbering
}
\description{
 Computes a re-arrangement required to transform the atom numbering
      of the given compound into the canonical atom numbering. This
      function uses the OBGraphSym and CanonicalLabels classes of Open
      Babel to compute the re-arrangement.

}
\usage{
canonicalNumbering(sdf)
}
\arguments{
  \item{sdf}{
	  Any sdfset object.
	}
}
\value{
   A list of vectors of index values. Each item in the list corresponds to one of the given
    compounds. The values of a list item are the re-arrangement of the atoms.
    For example, if the value in item 1, column 1 is  25, that means that
    atom number 1 in the original compound should become atom number 25
    in the canonical version of that compound.
}
\references{
	http://openbabel.org/api/2.3/canonical_code_algorithm.shtml
}
\author{
	Kevin Horan
}

\seealso{
\code{\link{canonicalize}}
}
\examples{
	\dontrun{
	data(sdfsample)
	labels = canonicalNumbering(sdfsample[1])
	}
}
\keyword{canonical}
\keyword{morgan}
