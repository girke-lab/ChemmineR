\name{loadSdf}
\alias{loadSdf}
\title{
	Load SDF Data
}
\description{
	Load an SDF formatted file or SDFSet objects into the database. This will also load arbitrary features 
	from the data as well as descriptor data. The \code{fct} parameter can be used to specify a function
	which will compute features which will then be indexed and stored in the database. These features can later
	be used to quickly search for compounds. Descriptors can also be computed and stored in another table. 
}
\usage{
loadSdf(conn, sdfFile, fct = function(x) data.frame(), descriptors=function(x) data.frame(descriptor=c(),descriptor_type=c()),
				Nlines = 10000, startline = 1, restartNlines = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{conn}{
	  A database connection object, such as is returned by \code{\link{initDb}}.
	}
  \item{sdfFile}{
	  Either the filename of an SDF formated file, or and SDFSet object.
	}
  \item{fct}{
	  A function to extract features from the data. It will be handed an SDFSet generated from the data being loaded.
	  This may be done in batches, so there is no guarantee that the given SDFSset will contain the whole dataset. This
	  function should return a data frame with a column for each feature and a row for each compound given, and in the
	  same order. Each of these features will become a new, indexed, table in the database which can be used later to
	  search for compounds. The column name will become the feature name.
		If not given, no features are computed.
	}
	\item{descriptors}{
		This function will also be given an SDFSet object, which may be done in batches. It should return a data frame
		with the following two columns: "descriptor" and "descriptor_type".  The "descriptor"  column should contain
		a string representation of the descriptor, and "descriptor_type" is the type of the descriptor. Our convention
		for atom pair is "ap" and "fp" for finger print. The order should be maintained.
			If not given no descriptors are computed.
	}
  \item{Nlines}{
	  When reading data from a file, the number of lines to read at a time. See also \code{\link{sdfStream}}.
	}
  \item{startline}{
	  When reading data from a file, the line number to start reading from.See also \code{\link{sdfStream}}
  }
  \item{restartNlines}{
	  When reading data from a file and startline > 1, the number of lines to look forward to find the start of the
	  next compound. See also \code{\link{sdfStream}}
  }
}
\details{
	New features can also be added using this function. However, all compounds must have all features so if new
	features are added to a new set of compounds, all existing features must be computable by the \code{fct} function 
	given. If new features are detected, all existing compounds will be run through \code{fct} in order to compute
	the new features for them as well.

	For example, if dataset X is loaded with features F1 and F2, and then at a later time we load dataset Y with
	new feature F3, the \code{fct} function used to load dataset Y must compute and return features F1, F2, and F3. 
	\code{loadSdf} will call \code{fct} with both datasets X and Y so that all features are available for all
	compounds. If any features are missing an error will be raised.

	If just new features are being added, but no new compounds, use the \code{\link{addNewFeatures}} function.
}
\value{
	Returns the compound id numbers of each compound loaded. These can be used to retrieve compounds later.
	These are id numbers computed by the database and are not extracted from the compound data itself.
}
\author{
	Kevin Horan
}


\seealso{
	\code{\link{sdfStream}}
}
\examples{
	  
   #create and initialize a new SQLite database
   conn = initDb("test.db")

	data(sdfsample)

	#just load the data with no features or descriptors
	ids=loadSdf(conn,sdfsample)
	unlink("test.db")

   conn = initDb("test.db")
	#load data and compute 3 features: molecular weight, with the MW function, 
	# and counts for RINGS and AROMATIC, as computed by rings, which returns a data frame itself.
	ids=loadSdf(conn,sdfsample,
			  function(sdfset) 
					data.frame(MW = MW(sdfset),  rings(sdfset,type="count",upper=6, arom=TRUE))
			 )
	unlink("test.db")
}
