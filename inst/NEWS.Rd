\name{NEWS}
\title{ChemmineR News}
\section{CHANGES IN VERSION 3.40.0}{
   \subsection{NEW FEATURES}{
		\itemize{
			\item Prior to this version, V3000 sdf files could be read, but would be written back
			out in V2000 format. Now these files will be written back out as V3000 files.
		}
	}
}

\section{CHANGES IN VERSION 2.32.0}{
   \subsection{NEW FEATURES}{
		\itemize{
			\item SDFDataTable function, allows viewing compound image and
			data in a web browser.
			\item read.SDFset can now skip over compounds with syntax
			errors in sdf files.
			\item The functions 'getIds', 'searchString', and 'searchSim'
			have been modified to query pubchem directly rather than
			going through the intermediate web service 'ChemmineTools'.
		}
	}
}
\section{CHANGES IN VERSION 2.30.0}{

   \subsection{NEW FEATURES}{
		\itemize{
			\item Plots can be generated using OpenBabel library (requires ChemmineOB package and OpenBabel)
			\item Functions to query data from PubChem directly. 
		}
	}
}

\section{CHANGES IN VERSION 2.20.0}{

   \subsection{NEW FEATURES}{
		\itemize{
			\item Fingerprint search now provides E-values
			\item New SDF plotting function, draw_sdf
		}
	}
}
\section{CHANGES IN VERSION 2.18.0}{

   \subsection{NEW FEATURES}{
		\itemize{
			\item Faster atompair generation
			\item Support for viewing features from SQL database
			\item Added function to generate 3D coordinates
		}
	}
}

\section{CHANGES IN VERSION 2.16.0}{

   \subsection{NEW FEATURES}{
		\itemize{
			\item SMARTS Search availible through ChemmineOB
			\item Folding of FPset objects
			\item Support for SQL database updates by compound name
			\item Coordinate re-generation via ChemmineOB improves structure rendering
		}
	}
}
\section{CHANGES IN VERSION 2.14.0}{

   \subsection{NEW FEATURES}{
      \itemize{
		  \item Integration of OpenBabel functionalities via new ChemmineOB add-on package
		  \item Improved SMILES support via new SMIset object class and SMILES import/export functions
		  \item Many other compound formats are now supported via ChemmineOB
      }
   }
}

\section{CHANGES IN VERSION 2.12.0}{

   \subsection{NEW FEATURES}{
      \itemize{
         \item Accelerated similarity searching of large small molecule data sets via new eiR add-on package
         \item Jarvis-Patrick clustering of large small molecule data sets 
         \item SQLite support for small molecule management 
      }
   }
}

\section{CHANGES IN VERSION 2.10.0}{

   \subsection{NEW FEATURES}{
      \itemize{
          \item  Streaming functionality for SDFs enables processing of millions 
            of molecules on a laptop
          \item  Fast and memory efficient fingerprint searches with atom pair 
            fingerprints or PubChem fingerprints
          \item  Flexible maximum common substructure (MCS) search support provided by 
            new fmcs.R add-on package
      }
   }
}

