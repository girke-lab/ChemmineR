\name{findCompoundsByName}
\alias{findCompoundsByName}
\title{
	Find compound by name
}
\description{
	Find the ids of compounds given the names.
}
\usage{
findCompoundsByName(conn, names, keepOrder = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{conn}{
	  A database connection object, such as is returned by \code{\link{initDb}}.
   }
  \item{names}{
     A list of names of compounds to search for. The names are those that would be
     returned by \code{\link{sdfid}}. An error will be raised if any names
     are not found.

   }
  \item{keepOrder}{
     If true, the order of the output compound ids will be the same as the 
     input names. This imposes a performance hit that can be significant for
     large datasets, thus it should be left FALSE unless needed.
   }
}
\value{
   Returns the compound ids for compounds with the given name. The output
   order is not guaranteed unless keepOrder is set to TRUE.
   An error will be raised if any name cannot be found.
}
\author{
   Kevin Horan
}

\examples{
   #create and initialize a new SQLite database
   conn = initDb("test4.db")

	data(sdfsample)

	#just load the data with no features or descriptors
	ids=loadSdf(conn,sdfsample)

   # find id of compound 650003
   findCompoundsByName(conn,c("650003"))
	unlink("test4.db")
}

