\name{findCompounds}
\alias{findCompounds}
\title{
	Find Compounds in Database
}
\description{
   Searches the SQL database using features computed at load time. Each feature
   used should be specified in the \code{featureNames} parameter. Then a set of filters
   can be given to search for specific compounds.
}
\usage{
findCompounds(conn, featureNames, tests)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{conn}{
	  A database connection object, such as is returned by \code{\link{initDb}}.
   }
  \item{featureNames}{
     A list of all feature names used in any test.
   }
  \item{tests}{
     A vector of filters that must all be true for a compound to be returned. For example:
     c("MW <= 400","RINGS > 3") would return all compounds with a molecular weight of 400 or
     less and a more than 3 rings, assuming these features exist in the database. 
     The syntax for each test is "<feature name> <SQL operator> <value>".
     These tests will simply be concatenated together with " AND " in-between them and tacked on the
     end of a WHERE clause of an SQL statement. So any SQL that will work in that context is fine.
   }
}
\value{
   Returns a list of compound ids. The actual compounds can be fetched with \code{\link{getCompounds}}.
}
\author{
   Kevin Horan
}
\seealso{
\code{\link{getCompounds}}
}
\examples{
   #create and initialize a new SQLite database
   conn = initDb("test1.db")

	data(sdfsample)

	#load data and compute 3 features: molecular weight, with the MW function, 
	# and counts for RINGS and AROMATIC, as computed by rings, which returns a data frame itself.
	ids=loadSdf(conn,sdfsample,
			  function(sdfset) 
					data.frame(MW = MW(sdfset),  rings(sdfset,type="count",upper=6, arom=TRUE))
			 )
   #search for compounds with molecular weight less than 200
   lightIds = findCompounds(conn,"MW",c("MW < 200"))
   MW(getCompounds(conn,lightIds)) # should find one compound with weight 140
	unlink("test1.db")

}
