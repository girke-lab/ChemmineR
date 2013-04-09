\name{getCompounds}
\alias{getCompounds}
\title{
	Get Compounds From Database
}
\description{
	Create SDF objects from the given set of compound ids. Id numbers can be found using the findCompounds function.
}
\usage{
   getCompounds(conn,compoundIds,filename=NA)
}
\arguments{
  \item{conn}{
	  A database connection object, such as is returned by \code{\link{initDb}}.
   }
  \item{compoundIds}{
     A vector of compound ids, as returned by \code{\link{loadSdf}} or \code{\link{findCompounds}}.
   }
   \item{filename}{
      If given, writes the compounds directly to the file named.
   }
}
\value{
   An SDFset with the requested compounds or nothing if \code{filename} was specified. 
   A warning will be raised if not all compounds could be found.
}
\author{
   Kevin Horan
}
\seealso{
   \code{\link{loadSdf}} 
   \code{\link{findCompounds}}.
}
\examples{
   #create and initialize a new SQLite database
   conn = initDb("test.db")

	data(sdfsample)

	#just load the data with no features or descriptors
	ids=loadSdf(conn,sdfsample)

   #returns a SDFset with 3 compounds
   getCompounds(conn, ids[1:3])

	unlink("test.db")

}
