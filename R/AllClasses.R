################################################
## Class and Method Definitions for ChemmineR ##
################################################
## SDF format definition: 
        # http://www.epa.gov/ncct/dsstox/MoreonSDF.html
	# http://www.symyx.com/downloads/public/ctfile/ctfile.jsp

#################################################
## (1) Class and Method Definitions for SDFstr ##
#################################################
## Download PubChem SDF samples for code testing 
## This function is not intended to be used by users.
.sdfDownload <- function(mypath="ftp://ftp.ncbi.nih.gov/pubchem/Compound/CURRENT-Full/SDF/", myfile="Compound_00650001_00675000.sdf.gz") {
	system(paste("wget ", mypath, myfile, sep=""))
	system(paste("gunzip ", myfile))
}
# .sdfDownload(mypath="ftp://ftp.ncbi.nih.gov/pubchem/Compound/CURRENT-Full/SDF/", myfile="Compound_00650001_00675000.sdf.gz") 

## Import SD File and Return SDFstr Class (list-like) 
read.SDFstr <- function(sdfstr) { 
	mysdf <- readLines(sdfstr) # reads file line-wi se into vector
        y <- regexpr("^\\${4,4}", mysdf, perl=TRUE) # identifies all fields that start with a '$$$$' sign
        index <- which(y!=-1)
        indexDF <- data.frame(start=c(1, index[-length(index)]+1), end=index)
        mysdf_list <- apply(indexDF, 1, function(x) mysdf[seq(x[1], x[2])])
        if(class(mysdf_list) != "list") { mysdf_list <- list(as.vector(mysdf_list))}
	names(mysdf_list) <- 1:length(mysdf_list)
        mysdf_list <- new("SDFstr", a=mysdf_list)
	return(mysdf_list)
}

## Define SDFstr class
setClass("SDFstr", representation(a = "list"))

## Accessor method for SDFstr
setGeneric(name="sdfstr2list", def=function(x) standardGeneric("sdfstr2list"))
setMethod(f="sdfstr2list", signature="SDFstr", definition=function(x) {return(x@a)}) 

## Replacement method for SDFstr using accessor method
setGeneric(name="sdfstr2list<-", def=function(x, value) standardGeneric("sdfstr2list<-"))
setReplaceMethod(f="sdfstr2list", signature="SDFstr", definition=function(x, value) {
	x@a <- value 
	return(x)
})

## Replacement method for SDFstr using "[" operator 
## It doesn't provide here full set of expected functionalities.
setReplaceMethod(f="[", signature="SDFstr", definition=function(x, i, j, value) {
	x@a[i] <- sdfstr2list(value)
	return(x)
})

## Define behavior of "[" operator for SDFstr 
setMethod(f="[", signature="SDFstr", definition=function(x, i, ..., drop) {
		x@a <- x@a[i]                   
		return(x)
})

## Behavior of "[[" operator to convert single SDFstr component to character vector
setMethod(f="[[", signature="SDFstr",
	definition=function(x, i, ..., drop) {
		return(x@a[[i]])                 
})

## Replacement method for SDFstr using "[" operator 
setReplaceMethod(f="[", signature="SDFstr", definition=function(x, i, value) {
	x@a[i] <- value
	return(x)
})

## Replacement method for SDFstr using "[[" operator 
setReplaceMethod(f="[[", signature="SDFstr", definition=function(x, i, value) {
	x@a[[i]] <- value
	return(x)
})

## Define print behavior for SDFstr
setMethod(f="show", signature="SDFstr",
	definition=function(object) { 
		cat("An instance of ", "\"", class(object), "\"", " with ", length(sdfstr2list(object)), " molecules", "\n", sep="")
})

## Length function
setMethod(f="length", signature="SDFstr",
    definition=function(x) {
    	return(length(sdfstr2list(x)))
})

## Write SDFstr Class to SD File 
write.SDF <- function(sdf, file, cid=FALSE, ...) {
	if(class(sdf)=="SDF") sdfstr <- as(sdf, "SDFstr")
	if(class(sdf)=="SDFstr") sdfstr <- sdf
	if(class(sdf)=="SDFset") {
		if(cid==TRUE) {	sdflist <- lapply(cid(sdf), function(x) sdf2str(sdf=sdf[[x]], cid=x, ...)) }	
		if(cid==FALSE) { sdflist <- lapply(cid(sdf), function(x) sdf2str(sdf=sdf[[x]], ...)) } 
		sdfstr <- as(sdflist, "SDFstr") 
	}
	cat(unlist(sdfstr2list(sdfstr)), sep="\n", file=file)
} 

## Coerce Methods for SDFstr Class
## SDFstr to list
setAs(from="SDFstr", to="list", 
	def=function(from) { 
		sdfstr2list(from)
})

## List to SDFstr
setAs(from="list", to="SDFstr", 
	def=function(from) {
		new("SDFstr", a=from)
})

## List to SDFstr
setAs(from="character", to="SDFstr", 
	def=function(from) {
		new("SDFstr", a=list(from))
})

################################################
## (2) Class and Method Definitions for SDF ##
################################################
setClass("SDF", representation(header="character", atomblock="matrix", bondblock="matrix", datablock="character"))

## Convert SDFstr to SDF Class
## SDFstr Parser Function
.sdfParse <- function(sdf, datablock=TRUE, tail2vec=TRUE, ...) {
	countpos <- grep("V\\d\\d\\d\\d$", sdf, perl=TRUE)
	if(length(countpos)==0) { countpos <- grep("V {0,}\\d\\d\\d\\d$", sdf, perl=TRUE) }
	if(length(countpos)==0) { countpos <- 4 }
	countline <- sdf[countpos]
        if(nchar(gsub("\\d| ", "", substring(countline, 1, 6))) != 0) { countline <- "  0  0" } # Create dummy countline if it contains non-numeric values 
	Natom <- as.numeric(substring(countline, 1, 3))
	Nbond <- as.numeric(substring(countline, 4, 6))
	start <- c(header=1, atom=countpos+1, bond=countpos+Natom+1, extradata=countpos+Natom+Nbond+1)
	index <- cbind(start=start, end=c(start[2:3], start[4], length(sdf)+1)-1)
	
	## Header block
	header <- sdf[index["header",1]:index["header",2]]
	if(length(header)==4) names(header) <- c("Molecule_Name", "Source", "Comment", "Counts_Line")	
	
	## Atom block
	ab2matrix <- function(ct=sdf[index["atom",1]:index["atom",2]]) {
		if((index["atom","end"] - index["atom","start"]) < 1) {
                        ctma <- matrix(rep(0,2), 1, 2, dimnames=list("0", c("C1", "C2"))) # Creates dummy matrix in case there is none.
                } else {
                        ct <- gsub("^ {1,}", "", ct)
		        ctlist <- strsplit(ct, " {1,}")
		        ctma <- matrix(unlist(ctlist), ncol=length(ctlist[[1]]), nrow=length(ct), byrow=TRUE)
		        myrownames <- paste(ctma[,4], 1:length(ctma[,4]), sep="_")
		        Ncol <- length(ctlist[[1]])
		        ctma <- matrix(as.numeric(ctma[,-4]), nrow=length(ct), ncol=Ncol-1, dimnames=list(myrownames, paste("C", c(1:3, 5:Ncol), sep="")))	
		}
                return(ctma)
	}
	atomblock <- ab2matrix(ct=sdf[index["atom",1]:index["atom",2]])
	
	## Bond block
	bb2matrix <- function(ct=sdf[index["bond",1]:index["bond",2]]) {
		#if((index["bond","end"] - index["bond","start"]) < 1) {
		if(((index["bond","end"] - index["bond","start"])+1) < 1) {
                        ctma <- matrix(rep(0,2), 1, 2, dimnames=list("0", c("C1", "C2"))) # Creates dummy matrix in case there is none.
                } else {
                    ct <- gsub("^(...)(...)(...)(...)(...)(...)(...)", "\\1 \\2 \\3 \\4 \\5 \\6 \\7", ct)
                    ct <- gsub("^ {1,}", "", ct)
                    ctlist <- strsplit(ct, " {1,}")
                    ctma <- matrix(unlist(ctlist), ncol=length(ctlist[[1]]), nrow=length(ct), byrow=TRUE)
                    Ncol <- length(ctlist[[1]])
                    ctma <- matrix(as.numeric(ctma), nrow=length(ct), ncol=Ncol, dimnames=list(1:length(ct), paste("C", 1:Ncol, sep="")))	
		}
                return(ctma)
	}
	bondblock <- bb2matrix(ct=sdf[index["bond",1]:index["bond",2]])
	
	## SDF name/value block
	ex2vec <- function(extradata=sdf[index["extradata",1]:index["extradata",2]]) {
                exstart <- grep("^>", extradata)
		if(length(exstart)==0) { 
                        exvec <- vector("character", length=0) 
                } else {
                    names(exstart) <- gsub("^>.*<|>", "", extradata[exstart])
                    exindex <- cbind(start=exstart, end=c(exstart[-1], length(extradata)+1)-1)
                    exvec <- sapply(rownames(exindex), function(x) paste(extradata[(exindex[x,1]+1):(exindex[x,2]-1)], collapse=" __ "))
		}
                return(exvec)
	}
	
	if(tail2vec==TRUE) {
		extradata <- ex2vec(extradata=sdf[index["extradata",1]:index["extradata",2]])
	} else {	
		extradata <- sdf[index["extradata",1]:index["extradata",2]]
	}

	## Assemble components in object of class SDF
	if(datablock==TRUE) {
		sdf <- new("SDF", header=header, atomblock=atomblock, bondblock=bondblock, datablock=extradata)
	} else {
		sdf <- new("SDF", header=header, atomblock=atomblock, bondblock=bondblock)
	}
	return(sdf)
}

## Accessor methods for SDF class
setGeneric(name="sdf2list", def=function(x) standardGeneric("sdf2list"))
setMethod(f="sdf2list", signature="SDF", definition=function(x) {return(list(header=header(x), atomblock=atomblock(x), bondblock=bondblock(x), datablock=datablock(x)))}) 
setGeneric(name="header", def=function(x) standardGeneric("header"))
setMethod(f="header", signature="SDF", definition=function(x) {return(x@header)}) 
setGeneric(name="sdfid", def=function(x, tag=1) standardGeneric("sdfid"))
setMethod(f="sdfid", signature="SDF", definition=function(x, tag=1) {return(x@header[tag])}) 
setGeneric(name="atomblock", def=function(x) standardGeneric("atomblock"))
setMethod(f="atomblock", signature="SDF", definition=function(x) {return(x@atomblock)}) 
setGeneric(name="atomcount", def=function(x, addH=FALSE, ...) standardGeneric("atomcount"))
setMethod(f="atomcount", signature="SDF", definition=function(x, addH=FALSE, ...) {
	if(addH==TRUE) { 
		return(table(c(gsub("_.*", "", rownames(x@atomblock)), rep("H", bonds(x, type="addNH")))))
	} else {
		return(table(gsub("_.*", "", rownames(x@atomblock))))
	}
})
setGeneric(name="bondblock", def=function(x) standardGeneric("bondblock"))
setMethod(f="bondblock", signature="SDF", definition=function(x) {return(x@bondblock)}) 
setGeneric(name="datablock", def=function(x) standardGeneric("datablock"))
setMethod(f="datablock", signature="SDF", definition=function(x) {return(x@datablock)}) 
setGeneric(name="datablocktag", def=function(x, tag) standardGeneric("datablocktag"))
setMethod(f="datablocktag", signature="SDF", definition=function(x, tag) {return(x@datablock[tag])}) 

## Define print behavior for SDF
setMethod(f="show", signature="SDF",                
   definition=function(object) {
         cat("An instance of ", "\"", class(object), "\"", "\n", sep="")
         cat("\n<<header>>", "\n", sep="")
         print(header(object))
         cat("\n<<atomblock>>", "\n", sep="")
         if(length(atomblock(object)[,1])>=5) {
             print(as.data.frame(rbind(atomblock(object)[1:2,], ...=rep("...", length(atomblock(object)[1,])), 
             atomblock(object)[(length(atomblock(object)[,1])-1):length(atomblock(object)[,1]),])))
         } else {
             print(atomblock(object))}
         cat("\n<<bondblock>>", "\n", sep="")
         if(length(bondblock(object)[,1])>=5) {
             print(as.data.frame(rbind(bondblock(object)[1:2,], ...=rep("...", length(bondblock(object)[1,])), 
             bondblock(object)[(length(bondblock(object)[,1])-1):length(bondblock(object)[,1]),])))
         } else {
             print(bondblock(object))}
         cat("\n<<datablock>> (", length(datablock(object)), " data items)", "\n", sep="")
         if(length(datablock(object))>=5) {
         	print(c(datablock(object)[1:4], "..."))
         } else {
         	print(datablock(object))}
})

## Behavior of "[" operator for SDF
setMethod(f="[", signature="SDF",
	definition=function(x, i, ..., drop) {
		return(sdf2list(x)[i]) 
})

## Replacement method for SDF using "[" operator 
setReplaceMethod(f="[", signature="SDF", definition=function(x, i, j, value) {
	if(i==1) x@header <- value
	if(i==2) x@atomblock <- value 
	if(i==3) x@bondblock <- value
	if(i==4) x@datablock <- value
	if(i=="header") x@header <- value
	if(i=="atomblock") x@atomblock <- value 
	if(i=="bondblock") x@bondblock <- value
	if(i=="datablock") x@datablock <- value
	return(x)
})

## Behavior of "[[" operator for SDF
setMethod(f="[[", signature="SDF",
	definition=function(x, i, ..., drop) {
		return(sdf2list(x)[[i]]) 
})

## Replacement method for SDF using "[[" operator 
setReplaceMethod(f="[[", signature="SDF", definition=function(x, i, j, value) {
	if(i==1) x@header <- value
	if(i==2) x@atomblock <- value 
	if(i==3) x@bondblock <- value
	if(i==4) x@datablock <- value
	if(i=="header") x@header <- value
	if(i=="atomblock") x@atomblock <- value 
	if(i=="bondblock") x@bondblock <- value
	if(i=="datablock") x@datablock <- value
	return(x)
})

## Coerce Methods for SDF Class 
## Character vector to SDF
setAs(from="character", to="SDF", 
	def=function(from) {
		.sdfParse(sdf=from) 
})

## SDF to list
setAs(from="SDF", to="list", 
	def=function(from) {
		list(header=header(from), atomblock=atomblock(from), bondblock=bondblock(from), datablock=datablock(from))
})

## SDF to SDFstr
setAs(from="SDF", to="SDFstr", 
	def=function(from) {
		new("SDFstr", a=list(as(from, "character")))		
})

## list (w. SDF components) to SDF  
setAs(from="list", to="SDF", 
	def=function(from) {
		new("SDF", bondblock=from$bondblock, header=from$header, atomblock=from$atomblock, datablock=from$datablock)
})

################################################# 
## (3) Class and Method Definitions for SDFset ##
#################################################
setClass("SDFset", representation(SDF="list", ID="character"))

## Store many SDFs in SDFset Object
read.SDFset <- function(sdfstr=sdfstr, ...) {
	## If a file name is provided run read.SDFstr function first
	if(is.character(sdfstr)) { 
		sdfstr <- read.SDFstr(sdfstr=sdfstr) 
	}
	## Iterate over SDFstr components	
	sdfset <- lapply(seq(along=sdfstr@a), function(x) .sdfParse(sdfstr2list(sdfstr)[[x]], ...))
	sdfset <- new("SDFset", SDF=sdfset, ID=paste("CMP", seq(along=sdfset), sep=""))
        ## Validity check of SDFs based on atom/bond block column numbers
        badsdf <- sum(!validSDF(sdfset))
        if(sum(badsdf)!=0) warning(paste(c(sum(badsdf), " invalid SDFs detected. To fix, run: valid <- validSDF(sdfset); sdfset <- sdfset[valid]")))
        return(sdfset)
}

## Accessor methods for SDFset class
setGeneric(name="SDFset2list", def=function(x) standardGeneric("SDFset2list"))
setMethod(f="SDFset2list", signature="SDFset", definition=function(x) {
	SDFlist <- x@SDF
	charlist <- lapply(seq(along=SDFlist), function(x) as(SDFlist[[x]], "list"))
	names(charlist) <- x@ID
	return(charlist)
}) 
setGeneric(name="SDFset2SDF", def=function(x) standardGeneric("SDFset2SDF"))
setMethod(f="SDFset2SDF", signature="SDFset", definition=function(x) { tmp <- x@SDF; names(tmp) <- x@ID; return(tmp)}) 
setMethod(f="header", signature="SDFset", definition=function(x) { return(lapply(SDFset2SDF(x), header))}) 
setMethod(f="sdfid", signature="SDFset", definition=function(x, tag=1) {return(as.vector(sapply(SDFset2SDF(x), sdfid, tag)))}) 
setGeneric(name="cid", def=function(x) standardGeneric("cid"))
setMethod(f="cid", signature="SDFset", definition=function(x) {return(x@ID)}) 
setMethod(f="atomblock", signature="SDFset", definition=function(x) {return(lapply(SDFset2SDF(x), atomblock))}) 
setMethod(f="atomcount", signature="SDFset", definition=function(x, addH, ...) {
	atomcounts <- lapply(SDFset2SDF(x), function(y) atomcount(y, addH, ...))
	return(atomcounts)
}) 
setMethod(f="bondblock", signature="SDFset", definition=function(x) {return(lapply(SDFset2SDF(x), bondblock))}) 
setMethod(f="datablock", signature="SDFset", definition=function(x) {return(lapply(SDFset2SDF(x), datablock))}) 
setMethod(f="datablocktag", signature="SDFset", definition=function(x, tag) {return(as.vector(sapply(SDFset2SDF(x), datablocktag, tag)))}) 

## Replacement method for SDF component of SDFset using accessor methods
setGeneric(name="SDFset2SDF<-", def=function(x, value) standardGeneric("SDFset2SDF<-"))
setReplaceMethod(f="SDFset2SDF", signature="SDFset", definition=function(x, value) {
	x@SDF <- value 
	return(x)
})

## Replacement method for ID component of SDFset using accessor methods
setGeneric(name="cid<-", def=function(x, value) standardGeneric("cid<-"))
setReplaceMethod(f="cid", signature="SDFset", definition=function(x, value) {
	x@ID <- value 
	if(any(duplicated(x@ID))) { 
		warning("The values in the CMP ID slot are not unique anymore. To fix this, run: cid(sdfset) <- makeUnique(cid(sdfset))")
	}
	return(x)
})

## Replacement method for SDFset using "[" operator 
## It doesn't provide here full set of expected functionalities.
setReplaceMethod(f="[", signature="SDFset", definition=function(x, i, j, value) {
	x@SDF[i] <- SDFset2SDF(value)
	return(x)
})

## Behavior of "[" operator for SDFset 
setMethod(f="[", signature="SDFset", definition=function(x, i, ..., drop) {
	if(is.logical(i)) {
                i <- which(i)
        }
        if(is.character(i)) { 
		ids <-seq(along=x@ID); names(ids) <- x@ID
		i <- ids[i]
	}
	x@SDF <- x@SDF[i]                   
	x@ID <- x@ID[i]                   
	if(any(duplicated(i))) {
		warning("The values in the CMP ID slot are not unique anymore. To fix this, run: cid(sdfset) <- makeUnique(cid(sdfset))")
	}
	return(x)
})

## Replacement method for SDFset using "[" operator 
setReplaceMethod(f="[", signature="SDFset", definition=function(x, i, value) {
	x@SDF[i] <- value
	return(x)
})

## Replacement method for SDFset using "[[" operator 
setReplaceMethod(f="[[", signature="SDFset", definition=function(x, i, value) {
	x@SDF[[i]] <- value
	return(x)
})

## Behavior of "[[" operator for SDFset to convert single SDFset component to SDF 
setMethod(f="[[", signature="SDFset", definition=function(x, i, ..., drop) {
	if(is.character(i)) { i <- which(x@ID %in% i) }
	return(x@SDF[[i]])                 
})

## Batch replacement of header, atomblock, bondblock and datablock sections for
## one to all SDF objects in an SDFset.
.gsubsdfsec <- function(sdfset, what, secdata) {
	if(class(secdata) %in% c("data.frame", "matrix")) { 
		secdata <- as.matrix(secdata)
		tmp <- as.list(seq(along=secdata[,1]))
		for(i in seq(along=tmp)) {
			vec <- as.character(secdata[i,])
			names(vec) <- colnames(secdata)
			tmp[[i]] <- vec
		}
		secdata <- tmp
	}
	if(length(sdfset) != length(secdata)) { stop("The length of the two data sets needs to match.") }
	sdfsec <- list(header=header(sdfset), atomblock=atomblock(sdfset), bondblock=bondblock(sdfset), datablock=datablock(sdfset))
	sdfsec[[what]] <- secdata
	sdflist <- lapply(seq(along=sdfsec[[1]]), function(x) as(list(header=sdfsec$header[[x]], atomblock=sdfsec$atomblock[[x]], bondblock=sdfsec$bondblock[[x]], datablock=sdfsec$datablock[[x]]), "SDF")) 
	names(sdflist) <- cid(sdfset)
	return(as(sdflist, "SDFset"))
}

setGeneric(name="header<-", def=function(x, value) standardGeneric("header<-"))
setReplaceMethod(f="header", signature="SDFset", definition=function(x, value) {
	sdfset <- .gsubsdfsec(sdfset=x, what=1, secdata=value) 
	return(sdfset)
})

setGeneric(name="atomblock<-", def=function(x, value) standardGeneric("atomblock<-"))
setReplaceMethod(f="atomblock", signature="SDFset", definition=function(x, value) {
	sdfset <- .gsubsdfsec(sdfset=x, what=2, secdata=value) 
	return(sdfset)
})

setGeneric(name="bondblock<-", def=function(x, value) standardGeneric("bondblock<-"))
setReplaceMethod(f="bondblock", signature="SDFset", definition=function(x, value) {
	sdfset <- .gsubsdfsec(sdfset=x, what=3, secdata=value) 
	return(sdfset)
})

setGeneric(name="datablock<-", def=function(x, value) standardGeneric("datablock<-"))
setReplaceMethod(f="datablock", signature="SDFset", definition=function(x, value) {
	sdfset <- .gsubsdfsec(sdfset=x, what=4, secdata=value) 
	return(sdfset)
})


## Length function
setMethod(f="length", signature="SDFset",
    definition=function(x) {
    	return(length(SDFset2SDF(x)))
})

## Print behavior for SDFset
setMethod(f="show", signature="SDFset",
	definition=function(object) { 
		cat("An instance of ", "\"", class(object), "\"", " with ", length(object), " molecules", "\n", sep="")
		if(any(duplicated(object@ID))) {
			warning("The values in the CMP ID slot are not unique. To fix this, run: cid(sdfset) <- makeUnique(cid(sdfset))")
		}
})

## Concatenate function for SDFset
## Note: is currently limited to 2 arguments!
setMethod(f="c", signature="SDFset", definition=function(x, y) {
	sdflist1 <- as(x, "SDF")
	sdflist2 <- as(y, "SDF")
	sdflist <- c(sdflist1, sdflist2)
	sdfset <- as(sdflist, "SDFset")
	if(any(duplicated(cid(sdfset)))) {
		warning("The values in the CMP ID slot are not unique anymore, makeUnique() can fix this!")
	}
	return(sdfset)
})

## Convert SDF to SDFstr Class for Export to File
## Function allows to customize output via optional arguments 
setGeneric(name="sdf2str", def=function(sdf, head, ab, bb, db, cid=NULL, sig=FALSE, ...) standardGeneric("sdf2str"))
setMethod(f="sdf2str", signature="SDF", definition = function(sdf, head, ab, bb, db, cid=NULL, sig=FALSE, ...) {
	## Checks
	if(class(sdf)!="SDF") stop("Function expects molecule object of class SDF as input!")	
	
	## Header
	if(missing(head)) {
		head <- as.character(sdf[[1]])
		if(sig==TRUE) head[2] <- paste("ChemmineR-", format(Sys.time(), "%m%d%y%H%M"), "XD", sep="")	
		if(length(cid)==1) head[1] <- cid
	}
	
	## Atom block
	if(missing(ab)) {
		ab <- sdf[[2]]
		ab <- cbind(Indent="", format(ab[,1:3], width=9, justify="right"), A=format(gsub("_.*", "", rownames(ab)), width=1, justify="right"), Space="", format(ab[,-c(1:3)], width=2, justify="right"))
		ab <- sapply(seq(along=ab[,1]), function(x) paste(ab[x, ], collapse=" "))
	}

	## Bond block
	if(missing(bb)) {
		bb <- sdf[[3]]
		bb <- cbind(Indent="", format(bb, width=3, justify="right"))
		bb <- sapply(seq(along=bb[,1]), function(x) paste(bb[x, ], collapse=""))
	}

	## Data block
	if(missing(db)) {
		db <- sdf[[4]]
		if(length(db)>0) {
			dbnames <- paste("> <", names(db), ">", sep="")
			dbvalues <- as.character(db)
			db <- as.vector(rbind(dbnames, dbvalues, ""))
		} else {
			db <- NULL
		}
	}
	
	## Assemble in character vector
	sdfstrvec <- c(head, ab, bb, "M  END", db, "$$$$")
        return(sdfstrvec)
})

## Coerce Methods for SDFset Class 
## SDFset to list with lists of SDF sub-components
setAs(from="SDFset", to="list", 
	def=function(from) {
		SDFset2list(from)
})

## SDFset to list with many SDF objects (for summary view)
setAs(from="SDFset", to="SDF", 
	def=function(from) {
		SDFset2SDF(from)
})
setGeneric(name="view", def=function(x) standardGeneric("view"))
setMethod(f="view", signature="SDFset", definition=function(x) { as(x, "SDF") })

## SDFstr to SDFset
setAs(from="SDFstr", to="SDFset", 
	def=function(from) {
		read.SDFset(sdfstr=from) 
})

## SDF to SDFset of length one
setAs(from="SDF", to="SDFset", 
	def=function(from) {
		new("SDFset", SDF=list(from), ID="CMP") 
})

## SDF to SDFstr
setAs(from="SDF", to="SDFstr", 
	def=function(from) {
		as(as(from, "SDFset"), "SDFstr") 
})

## SDFset to character for SDFstr coercion
setAs(from="SDF", to="character", 
	def=function(from) {
       		sdf2str(sdf=from)
})

## SDFset to SDFstr
setAs(from="SDFset", to="SDFstr", 
	def=function(from) {
		from <- lapply(seq(along=from), function(x) as(from[[x]], "character"))
		new("SDFstr", a=from)		
})

## List of SDFs to SDFset
setAs(from="list", to="SDFset", 
	def=function(from) {
		new("SDFset", SDF=from, ID=names(from))
})


#######################################################
## (4) Class and Method Definitions for AP and APset ##
#######################################################
## Function to coerce SDF class to old non-S4 AP CMP object
SDF2apcmp <- function(SDF) { 
	atoms <- gsub("_.*", "", rownames(atomblock(SDF)))
	u <- as.numeric(bondblock(SDF)[,1])
	v <- as.numeric(bondblock(SDF)[,2])
	t <- as.numeric(bondblock(SDF)[,3])
	n_atoms <- length(atoms)
	n_bonds <- length(u)
	return(list(atoms=atoms, bonds=list(u=u, v=v, t=t), n_atoms=n_atoms, n_bonds=n_bonds))
}

## Define AP/APset S4 classes for single AP vector and AP list
setClass("AP", representation(AP="numeric"))
setClass("APset", representation(AP="list", ID="character"))

## Create instance of APset form SDFset 
sdf2ap <- function(sdfset) {
	if(!class(sdfset) %in% c("SDF", "SDFset")) stop("Functions expects input of classes SDF or SDFset.")
	if(class(sdfset)=="SDF") {
		return(new("AP", AP=.gen_atom_pair(SDF2apcmp(sdfset))))
	}
	if(class(sdfset)=="SDFset") {
		aplist <- as.list(seq(along=sdfset))
		for(i in seq(along=aplist)) {
			aplist[[i]] <- .gen_atom_pair(SDF2apcmp(sdfset[[i]]))
		}
		return(new("APset", AP=aplist, ID=cid(sdfset)))
	}
}

## Accessor methods for APset class
setGeneric(name="ap", def=function(x) standardGeneric("ap"))
setMethod(f="ap", signature="AP", definition=function(x) { return(x@AP) })
setMethod(f="ap", signature="APset", definition=function(x) { tmp <- x@AP; names(tmp) <- x@ID; return(tmp) })
setMethod(f="cid", signature="APset", definition=function(x) { return(x@ID) })

## Replacement method for ID component of APset using accessor methods.
## Note: generic for cid() is defined under SDFset class section.
setReplaceMethod(f="cid", signature="APset", definition=function(x, value) {
	x@ID <- value 
	if(any(duplicated(x@ID))) { 
		warning("The values in the CMP ID slot are not unique anymore. To fix this, run: cid(apset) <- makeUnique(cid(apset))")
	}
	return(x)
})

## Replacement method for APset using "[" operator 
## It doesn't provide here full set of expected functionalities.
setReplaceMethod(f="[", signature="APset", definition=function(x, i, j, value) {
	x@AP[i] <- ap(value)
	x@ID[i] <- cid(value)
	return(x)
})

## Behavior of "[" operator for APset 
setMethod(f="[", signature="APset", definition=function(x, i, ..., drop) {
	if(is.logical(i)) {
                i <- which(i)
        }
	if(is.character(i)) { 
		ids <- seq(along=x@ID)
		names(ids) <- x@ID
		i <- ids[i]
	}
	x@AP <- x@AP[i]                   
	x@ID <- x@ID[i]                   
	if(any(duplicated(i))) {
		warning("The values in the CMP ID slot are not unique anymore. To fix this, run: cid(apset) <- makeUnique(cid(apset))")
	}
	return(x)
})

## Replacement method for APset using "[[" operator 
setReplaceMethod(f="[[", signature="APset", definition=function(x, i, value) {
	x@AP[[i]] <- value
	return(x)
})

## Behavior of "[[" operator for SDFset to convert single SDFset component to SDF 
setMethod(f="[[", signature="APset", definition=function(x, i, ..., drop) {
	if(is.character(i)) { i <- which(x@ID %in% i) }
	return(new("AP", AP=x@AP[[i]]))                 
})

## Length function
setMethod(f="length", signature="APset",
    definition=function(x) {
    	return(length(ap(x)))
})

## Print behavior for APset
setMethod(f="show", signature="APset",
	definition=function(object) { 
		cat("An instance of ", "\"", class(object), "\"", " with ", length(object), " molecules", "\n", sep="")
		if(any(duplicated(object@ID))) {
			warning("The values in the CMP ID slot are not unique. To fix this, run: cid(apset) <- makeUnique(cid(apset))")
		}
})

## Concatenate function for APset
## Note: is currently limited to 2 arguments!
setMethod(f="c", signature="APset", definition=function(x, y) {
	aplist1 <- as(x, "list")
	aplist2 <- as(y, "list")
	aplist <- c(aplist1, aplist2)
	apset <- as(aplist, "APset")
	if(any(duplicated(cid(apset)))) {
		warning("The values in the CMP ID slot are not unique anymore, makeUnique() can fix this!")
	}
	return(apset)
})

## Define print behavior for AP
setMethod(f="show", signature="AP",                
   definition=function(object) {
         cat("An instance of ", "\"", class(object), "\"", "\n", sep="")
         cat("<<atom pairs>>", "\n", sep="")
         if(length(object@AP)>=5) {
             cat(c(object@AP[1:5], "... length:", length(object@AP), "\n"))
         } else {
             print(object@AP)}
})

## Coerce Methods for APset Class
## APset to list with many SDF objects
setAs(from="APset", to="list",
        def=function(from) {
                ap(from)
})

## List of APs to APset
setAs(from="list", to="APset", 
        def=function(from) {
                new("APset", AP=from, ID=names(from))
})

## APset to list with many AP objects (for summary view)
setAs(from="APset", to="AP", 
        def=function(from) {
                tmp <- lapply(seq(along=from), function(x) from[[x]])
		names(tmp) <- cid(from)
		return(tmp)
})
setMethod(f="view", signature="APset", definition=function(x) { as(x, "AP") })

## Coerce APset to old list-style descriptor database used by search/cluster functions
apset2descdb <- function(apset) {
                list(descdb=ap(apset), cids=cid(apset), sdfsegs=NULL, source="SDFset", type="SDFset")
}

###################
## (5) Utilities ##
###################

#################################################
## (5.1) Detect Invalid SDFs in SDFset Objects ##
#################################################
validSDF <- function(x, Nabcol = 3, Nbbcol = 3, logic="&") {
        if(class(x)!="SDFset") warning("x needs to be of class SDFset")
	ab <- atomblock(x); abcol <- sapply(names(ab), function(x) length(ab[[x]][1,]))
        bb <- bondblock(x); bbcol <- sapply(names(bb), function(x) length(bb[[x]][1,]))
        if(logic=="|") { validsdf <- abcol >= Nabcol | bbcol >= Nbbcol }
	if(logic=="&") { validsdf <- abcol >= Nabcol & bbcol >= Nbbcol }
        return(validsdf)
}

######################################################################
## (5.2) Create Unique CMP Names by Appending a Counter to Duplates ##
######################################################################
makeUnique <- function(x, silent=FALSE) {
	if(all(!duplicated(x))) {
		if(silent!=TRUE) print("No duplicates detected!")
		return(x)
	} else {
		count <- table(x); count <- count[x]
		dupids <- count[count>1]; dupids <- dupids[!duplicated(names(dupids))]
		for(i in seq(along=dupids)) {
			names(count)[names(count) %in% names(dupids[i])] <- paste(names(dupids)[i], seq(1, dupids[i]), sep="_")
		}
		if(silent!=TRUE) print(paste("Counter appended to", length(dupids), "duplicates!"))
		return(names(count))
	}
}

###############################
## (5.3) Molecule Properties ##
###############################
## (5.3.1) Atom count matrix
atomcountMA <- function(x, ...) {
        if(class(x)=="SDF") x <- as(x, "SDFset")
	atomcountlist <- atomcount(x, ...) 	
	columns <- unique(unlist(lapply(seq(along=atomcountlist), function(x) names(atomcountlist[[x]]))))
        myMA <- matrix(NA, length(atomcountlist), length(columns), dimnames=list(NULL, columns))
        for(i in seq(along=atomcountlist)) myMA[i, names(atomcountlist[[i]])] <- atomcountlist[[i]]
	myMA[is.na(myMA)] <- 0
	rownames(myMA) <- names(atomcountlist)
	return(myMA)
}

## (5.3.2) Molecular weight (MW data from http://iupac.org/publications/pac/78/11/2051/)
data(atomprop); atomprop <- atomprop # Import MW data frame from /data into workspace.
MW <- function(x, mw=atomprop, ...) {
        if(class(x)=="SDF") x <- as(x, "SDFset")
	## Create MW vector with atom symbols in name slot
	AW <- mw$Atomic_weight; names(AW) <- mw$Symbol
	
	## Calculate MW
	propma <- atomcountMA(x, ...)
	MW <- rowSums(t(t(propma) * AW[colnames(propma)]), na.rm = TRUE) 
	return(MW)
}

## (5.3.3) Molecular formula
MF <- function(x, ...) {
        if(class(x)=="SDF") x <- as(x, "SDFset")
	propma <- atomcountMA(x, ...)
	propma <- propma[c(1, seq(along=propma[,1])),] # Duplicates first row to support processing of single molecule with same code
	hillorder <- colnames(propma); names(hillorder) <- hillorder
	hillorder <- na.omit(unique(hillorder[c("C", "H", sort(hillorder))]))
	propma <- propma[, hillorder]
	propma[propma==1] <- ""
	MF <- paste(colnames(propma), t(propma), sep="")
	propma <- matrix(MF, nrow=length(propma[,1]), ncol=length(propma[1,]), dimnames=list(rownames(propma), colnames(propma)), byrow=TRUE)
	MF <- seq(along=propma[,1]); names(MF) <- rownames(propma)
	zeroma <- matrix(grepl("[\\*A-Za-z]0$", propma), nrow=length(propma[,1]), ncol=length(propma[1,]), dimnames=list(rownames(propma), colnames(propma)))
	propma[zeroma] <- ""
	for(i in seq(along=MF)) { MF[i] <-  paste(propma[i,], collapse="") }
	return(MF[-1]) # Minus one to remove duplicated entry in first row of propma
}

##############################################################
## (5.4) Convert SDF Tail to Numeric and Character Matrices ##
##############################################################
## (5.4.1) Store everything in one character matrix
datablock2ma <- function(datablocklist=datablock(sdfset), cleanup=" \\(.*", ...) {
	if(exists("cleanup")) for(i in seq(along=datablocklist)) names(datablocklist[[i]]) <- gsub(cleanup, "", names(datablocklist[[i]])) # Required if name tags contain compound ids
        columns <- unique(unlist(lapply(seq(along=datablocklist), function(x) names(datablocklist[[x]]))))
        myMA <- matrix(NA, length(datablocklist), length(columns), dimnames=list(NULL, columns))
        for(i in seq(along=datablocklist)) myMA[i, names(datablocklist[[i]])] <- datablocklist[[i]]
	rownames(myMA) <- names(datablocklist)
	return(myMA)
}
# Usage:
# blockmatrix <- datablock2ma(datablocklist=datablock(sdfset))

## (5.4.2) Split SDF tail matrix into character and numeric matrices
splitNumChar <- function(blockmatrix=blockmatrix) {
	# Define function to check for valid numeric values in a character vector
	numberAble <- function(myvec, type=c("single", "vector"), extras = c(".", "NA")) {
	    type <- match.arg(type)
	    old <- options(warn = -1)
	    on.exit(options(old))
	    myvec <- gsub(" ", "", myvec)
	    myvecsub <- myvec[!myvec %in% c("", extras)]
	    isnum <- !any(is.na(as.numeric(myvecsub)))
	    if (type == "single") 
		isnum
	    else if (isnum) 
		as.numeric(myvec)
	    else myvec
	}
        colindex <- sapply(seq(along=blockmatrix[1,]), function(x) numberAble(blockmatrix[,x]))
        numMA <- blockmatrix[ , colindex]
        storage.mode(numMA) <- "numeric"
        charMA <- blockmatrix[ , !colindex]
        return(list(numMA=numMA, charMA=charMA))
}
# Usage:
# numchar <- splitNumChar(blockmatrix=blockmatrix)

#########################
## (5.5) Bond Matrices ##
#########################

## (5.5.1) Generate bond matrix from SDFset or SDF objects
conMA <- function(x, exclude="none") {
        ## Function for SDF object 
        .conMA <- function(x, exclude=exclude) {
            atoms <- rownames(atomblock(x))
            conma <- matrix(0, length(atoms), length(atoms), dimnames=list(atoms, atoms))
            bondblock <- bondblock(x)
            for(i in seq(along=bondblock[,1])) conma[bondblock[i,1], bondblock[i,2]] <- bondblock[i,3]
            for(i in seq(along=bondblock[,1])) conma[bondblock[i,2], bondblock[i,1]] <- bondblock[i,3]
            index <- !gsub("_.*", "", rownames(conma)) %in% exclude
            conma <- conma[index,index]
            return(conma)
        }   
        ## Run on SDF objects
        if(class(x)=="SDF") {
                conma <- .conMA(x, exclude)
                return(conma)
        }
        ## Run on SDFset objects containing one or many molecules
        if(class(x)=="SDFset") {
                conma_set <- lapply(seq(along=x), function(y) .conMA(x[[y]], exclude))
                names(conma_set) <- cid(x)
                return(conma_set)
        }
}
# Usage:
# conma <- conMA(sdfset[1:2], exclude=c("H"))

## (5.5.2) Compute bond/charge count for each atom in SDFset or SDF objects
## This is used to add hydrogens with methods/functions atomcount, atomcountMA, MW and MF
bonds <- function(x, type="bonds") {
	## Input checks
        if(!any(c("SDF", "SDFset") %in% class(x))) stop("x needs to be of class SDF or SDFset")
	if(!any(c("bonds", "charge", "addNH") %in% type)) stop("type can only be assigned: bonds, charge or addNH") 
	
	## Compute bonds, charges and missing hydrogens
	.bonds <- function(x, type=type) {
		atomMA <- atomblock(x)
		atoms <- gsub("_.*", "", rownames(atomMA))
		bondMA <- bondblock(x)
		Nbonds1 <- cbind(atoms=c(bondMA[,1], bondMA[,2]), bonds=c(bondMA[,3], bondMA[,"C3"]))
		Nbonds1 <- tapply(Nbonds1[, "bonds"], Nbonds1[, "atoms"], sum)
		Nbonds <- rep(0, length(atomMA[,1])); names(Nbonds) <- seq(along=atomMA[,1]); Nbonds[names(Nbonds1)] <- Nbonds1 	
		
		## Valence related to position in periodic table (following octet rule) 
		val <- c("1"=1, "17"=1, "2"=2, "16"=2, "13"=3, "15"=3, "14"=4)
		group <- as.numeric(atomprop$Group); names(group) <- as.character(atomprop$Symbol)
		Nbondrule <- val[as.character(group[atoms])]
		Nbondrule[is.na(Nbondrule)] <- 0 # Atoms with undefined Nbondrule (NAs) are assigned zero
		Nbondrule[Nbondrule < Nbonds] <- Nbonds[Nbondrule < Nbonds] # Set Nbondrule to Nbonds values, where latter is larger 

		## Charges
		charge <- c("0"=0, "1"=3, "2"=2, "3"=1, "4"=0, "5"=-1, "6"=-2, "7"=-3) # 4 is "doublet_radical"
		charge <- charge[as.character(atomMA[,5])]
		Nbonds <- data.frame(atom=atoms, Nbondcount=Nbonds, Nbondrule=Nbondrule, charge=charge) 
		
		## Data type to return
		if(type=="bonds") { return(Nbonds) }
		if(type=="charge") {
			chargeindex <- Nbonds[, "charge"] != 0 
			if(sum(chargeindex) == 0) {
				return(NULL)
			} else {
				chargeDF <- Nbonds[chargeindex, ]
				charge <- chargeDF[, "charge"]
				names(charge) <- chargeDF[, "atom"] 
				return(charge)
			}
		}
		if(type=="addNH") { 
			Nbonds[Nbonds[,"Nbondcount"] >= Nbonds[,"Nbondrule"], "charge"] <- 0 # Ignore charge where Nbondcount greater or equal than Nbondrule
			Nbonds[Nbonds[,"Nbondcount"] == 0, c("Nbondrule", "charge")] <- 0 # Ignore atoms with zero bonds
			NH <- sum((Nbonds[, "Nbondrule"] + Nbonds[, "charge"]) - Nbonds[, "Nbondcount"])
			if(NH < 0) NH <- 0 # Count should not be negative
			return(NH) 
		}
	}
        ## Run on SDF objects
        if(class(x)=="SDF") {
                bonds <- .bonds(x, type)
                return(bonds)
        }
        ## Run on SDFset objects containing one or many molecules
        if(class(x)=="SDFset") {
                bonds_set <- lapply(seq(along=x), function(y) .bonds(x[[y]], type))
                names(bonds_set) <- cid(x)
		if(type=="bonds") {
                	return(bonds_set)
		}
		if(type=="charge") {
                	return(bonds_set)
		}
		if(type=="addNH") {
                	return(unlist(bonds_set))
		}
        }
}
# Usage:
# bondDF <- bonds(sdfset[1], type="df")) 
# bondcount <- bonds(sdfset[1], type="addNH")) 

#################################
## (5.6) String Search Method ##
#################################

## (5.6.1) String search function for SDFset
grepSDFset <- function(pattern, x, field="datablock", mode="subset", ignore.case=TRUE, ...) {
	## Generate search vector and index for desired field in SDFset
	if(field=="header" | field==1) {
		searchfield <- header(x)
		searchstr <- as.character(unlist(searchfield))
		index <- sapply(searchfield, length)
		index <- unlist(sapply(index, function(x) seq(1, x)))
		indexnames <- rep(seq(along=searchfield), sapply(searchfield, length))
		names(index) <- indexnames
	}
	if(field=="atomblock" | field==2) { 
		searchfield <- atomblock(x)
		searchstr <- as.character(unlist(sapply(searchfield, rownames)))
		searchstr <- gsub("_.*", "", searchstr)
		index <- sapply(searchfield, function(x) length(x[,1]))
		index <- unlist(sapply(index, function(x) seq(1, x)))
		indexnames <- rep(seq(along=searchfield), sapply(searchfield, function(x) length(x[,1])))
		names(index) <- indexnames
	}
	if(field=="bondblock" | field==3) { # Not very useful for string searching.
		searchfield <- bondblock(x)
		searchstr <- as.character(unlist(sapply(searchfield, rownames)))
		searchstr <- gsub("_.*", "", searchstr)
		index <- sapply(searchfield, function(x) length(x[,1]))
		index <- unlist(sapply(index, function(x) seq(1, x)))
		indexnames <- rep(seq(along=searchfield), sapply(searchfield, function(x) length(x[,1])))
		names(index) <- indexnames
	}
	if(field=="datablock" | field==4) { 
		searchfield <- datablock(x)
		searchstr <- paste(names(unlist(searchfield)), as.character(unlist(searchfield)), sep="___")
		index <- sapply(searchfield, length)
		index <- unlist(sapply(index, function(x) seq(1, x)))
		indexnames <- rep(seq(along=searchfield), sapply(searchfield, length))
		names(index) <- indexnames
	}
	## Search with grep
	matchpos <- grep(pattern, searchstr, ignore.case=ignore.case, ...)
	if(mode=="index") {
		return(index[matchpos])

	}
	if(mode=="subset") {
		xpos <- as.numeric(unique(names(index[matchpos])))
		return(as(x[xpos], "SDF"))
	}
}

##########################
## (6) Plotting Methods ##
##########################
## Plot single CMP Structure
plotStruc <- function(sdf, atomcex=1.2, atomnum=FALSE, no_print_atoms=c("C"), noHbonds=TRUE, bondspacer=0.12, ...) {
	toplot <- list(atomblock=atomblock(sdf)[,1:2], bondblock=bondblock(sdf)[,1:3])
	## Create empty plot with proper dimensions
	plot(toplot[[1]], type="n", axes=F, xlab="", ylab="", ...)
	## Remove C-hydrogens including their bonds 
	if(noHbonds==TRUE) {
		CHbondindex <- sapply(seq(toplot[[2]][,1]), function(x) paste(sort(gsub("_.*", "", rownames(toplot[[1]]))[toplot[[2]][x,1:2]]), collapse="") == "CH")
		toplot[[1]] <- toplot[[1]][sort(unique(as.numeric(toplot[[2]][!CHbondindex,1:2]))), ]
		toplot[[2]] <- toplot[[2]][!CHbondindex,] 
	}	
	## Plot bonds
	for(i in seq(along=toplot[[2]][,1])) {
		x <- toplot[[1]][gsub("*.*_", "", rownames(toplot[[1]])) %in% toplot[[2]][i,1:2],1]
		y <- toplot[[1]][gsub("*.*_", "", rownames(toplot[[1]])) %in% toplot[[2]][i,1:2],2]
		## Plot single bonds
		if(toplot[[2]][i,3]==1) {
			lines(x=x, y=y, lty=1, lwd=3) 
		} 
		## Plot double bonds
		if(toplot[[2]][i,3]==2) {
			rslope <- (atan(diff(y)/diff(x))*180/pi)/90
			lines(x=x-rslope*bondspacer, y=y+(1-abs(rslope))*bondspacer, lty=1, lwd=3)
			lines(x=x+rslope*bondspacer, y=y-(1-abs(rslope))*bondspacer, lty=1, lwd=3)
		}
		## Plot triple bonds
		if(toplot[[2]][i,3]==3) {
			rslope <- (atan(diff(y)/diff(x))*180/pi)/90
			bondspacer <- bondspacer * 2
			lines(x=x, y=y, lty=1, lwd=3) 
			lines(x=x-rslope*bondspacer, y=y+(1-abs(rslope))*bondspacer, lty=1, lwd=3)
			lines(x=x+rslope*bondspacer, y=y-(1-abs(rslope))*bondspacer, lty=1, lwd=3)
		}
	}
	## Exclude certain atoms from being printed
	exclude <- paste("(^", no_print_atoms, "_)", sep="", collapse="|")
	labelMA <- toplot[[1]][!grepl(exclude, rownames(toplot[[1]])), ] 
	if(is.vector(labelMA)) labelMA <- matrix(labelMA, 1, 2, byrow=TRUE, dimnames=list(rownames(toplot[[1]])[!grepl(exclude, rownames(toplot[[1]]))], c("C1", "C2")))
	if(is.matrix(labelMA) & length(labelMA[,1])>=1) {
		atomcol <- gsub("_.*", "", rownames(labelMA)); atomcol[!grepl("N|C|O|H", atomcol)] <- "any"; mycol <- c(C="black", H="black", N="blue", O="red", any="green"); atomcol <- mycol[atomcol]
		## Overplot nodes to display atom labels
		points(x=labelMA[,1], y=labelMA[,2], col="white", pch=16, cex=2.8)
		## Plot atom labels
		if(atomnum==TRUE) {
			text(x=labelMA[,1], y=labelMA[,2], gsub("_", "", rownames(labelMA)), cex=atomcex, col=atomcol) 
		} else {
			text(x=labelMA[,1], y=labelMA[,2], gsub("_.*", "", rownames(labelMA)), cex=atomcex, col=atomcol) 
		}
	}
}
# Usage:
# plotStruc(sdf=sdfset[[2]], atomcex=1.2, atomnum=F, no_print_atoms=c("C"), noHbonds=TRUE, bondspacer=0.08)
# par(mfrow=c(2,3)); for(i in 1:6) plotStruc(sdf=sdfset[[i]], atomcex=1.8, atomnum=F, no_print_atoms=c("C"), noHbonds=TRUE, bondspacer=0.08)

## Plot method for single SDF object
setMethod(f="plot", signature="SDF", definition=function(x, print=TRUE, ...) { 
		plotStruc(sdf=x, ...)
		if(print==TRUE) { return(x) } 
	}
)

## Plot method for multiple SDF objects in SDFset
setMethod(f="plot", signature="SDFset",
	definition=function(x, griddim, print_cid=cid(x), print=TRUE, ...) {
		if(missing(griddim)) {
			mydim <- ceiling(sqrt(length(x)))
			griddim <- c(mydim, mydim)
		}
		par(mfrow=griddim)
		for(i in 1:length(x)) { plotStruc(sdf=x[[i]], main=print_cid[i], ...) }
		if(print==TRUE) { return(SDFset2SDF(x)) }
})


############################################################################
## Overview of Objects, Methods and Functions for SDF and AP Classes in R ##
############################################################################
## Defined Classes
 	# SDFstr: list-like container for many SDFs where each one is stored line-by-line as character vector
	# SDF: single SDF stored in list with 4 components: header, atom block, bond block and sdf tail
	# SDFset: list-like container containing one to many SDF class objects and the same number of ids in ID slot 

## (A) Import SD File to SDFstr/SDFset Containers
	# library(ChemmineR)
	# sdfstr <- read.SDFstr("sdf_libs/MS_Spectrum.sdf")
	# sdfset <- read.SDFset(sdfstr)
	# sdfset <- read.SDFset("sdf_libs/MS_Spectrum.sdf")
	# data(sdfsample); sdfset <- sdfsample

## (B) Compound IDs
	# sdfid(sdfset[1:4]) # Retrieves CMP IDs from Molecule Name field in header block.
	# cid(sdfset[1:4]) # Retrieves CMP IDs from ID slot in SDFset.
	# unique_ids <- makeUnique(sdfid(sdfset)) # Creates unique IDs by appending a counter to duplicates.
	# cid(sdfset) <- unique_ids # Assigns uniquified IDs to ID slot.

## (C) Subsetting, Replacement and Related Methods
	# sdfstr[1:4]
	# sdfstr[[1]] # returns one sdf component as character vector
	# view(sdfset[1:4]) # Summary view of several SDFset components
	# sdf <- sdfset[[1]] # returns object of class SDF 
	# atomblock(sdf); sdf[[2]]; sdf[["atomblock"]] # all these methods return the same component
	# sdf[[2]][1] <- 999
	# sdfstr[1] <- sdfstr[1]
	# length(sdfset)
	# myids <- cid(sdfset)[sample(1:100, 10)] 
	# sdfset[myids]	

## (D) Miscellaneous Functions
        ## Batch retrieval of SDF components
	# header(sdfset[1:4]); atomblock(sdfset[1:4]); atomcount(sdfset[1:4]); bondblock(sdfset[1:4]); datablock(sdfset[1:4])
	
	## Create atom count matrix
	# propma <- atomcountMA(sdfset=sdfset)
	# boxplot(propma, main="Atom Frequency"); x11(); boxplot(rowSums(propma), main="All Atom Frequency")	

	## Compute MW and formula
	# MW(sdfset[1:4]); MF(sdfset[1:4])
	# propma <- data.frame(MF=MF(sdfset), MW=MW(sdfset), atomcountMA(sdfset)); propma[1:4,]
        # datablock(sdfset) <- propma # Works with all SDF components 
	# test <- apply(propma[1:4,], 1, function(x) data.frame(col=colnames(propma), value=x)); sdf.visualize(sdfset[1:4], extra = test)

	## Convert SDF tail to matrix
	# blockmatrix <- datablock2ma(datablocklist=datablock(sdfset)) # Converts data block to matrix	
	# numchar <- splitNumChar(blockmatrix=blockmatrix); numchar[[1]][1:4,] # Splits matrix to numeric matrix and character matrix
	
	## String Searching of SDFset
	# grepSDFset("Ampicillin", sdfset, field="datablock", mode="subset") # To return index, set mode="index")

## (E) Class Coercions
	## From SDFstr to list, SDF and SDFset
	# as(sdfstr[1:2], "list")  
	# as(sdfstr[[1]], "SDF") 
	# as(sdfstr[1:2], "SDFset") 
	## From SDF to SDFstr, SDFset, list with SDF sub-components
	# sdfcomplist <- as(sdf, "list")
	# sdfcomplist <- as(sdfset[1:4], "list"); as(sdfcomplist[[1]], "SDF")
	# sdflist <- as(sdfset[1:4], "SDF"); as(sdflist, "SDFset")
	# as(sdfset[[1]], "SDFstr")
	# as(sdfset[[1]], "SDFset")
	## From SDFset to lists with components consisting of SDF or sub-components
	# as(sdfset[1:4], "SDF")
	# as(sdfset[1:4], "list")
	# as(sdfset[1:4], "SDFstr")
	## Concatenation of several SDFsets
	# c(sdfset[1:4], sdfset[5:8]) # Note: c() is currently limited to two arguments!!

## (F) Export Custom SDFs
	## sdf2str(sdf=sdfset[[1]], sig=TRUE, cid=TRUE) # uses default components
	## sdf2str(sdf=sdfset[[1]], head=letters[1:4], db=NULL) # uses custom components for header and datablock

## (G) Write SDF, SDFset or SDFstr Classes to File
	# write.SDF(sdfset[1:4], file="sub.sdf", sig=TRUE, cid=TRUE, db=NULL)
	# write.SDF(sdfstr[1:4], file="sub.sdf")
	# cat(unlist(as(sdfstr[1:4], "list")), file="sub.sdf", sep="\n") 

## (H) Plot CMP Structures in R and on ChemMine Tools
	# plot(sdfset[[1]])
	# plot(sdfset[1:4])
	# sdf.visualize(sdfset[1:4]
	# extra <- apply(propma[1:4,], 1, function(x) data.frame(Property=colnames(propma), Value=x))
	# sdf.visualize(sdfset[1:4], extra=extra)

## (I) Usage of APset Class
	# source("sdfFct.R"); source("simFct.R"); source("clusterFct.R")
	# sdfset <- read.SDFset("sdf_libs/MS_Spectrum.sdf")
	# apset <- sdf2ap(sdfset[1:100])
	# cid(apset[1:4]); ap(apset[1:4])
	# tmp <- as(apset, "list"); as(tmp, "APset")
	# view(apset[1:4])
	# apset2descdb(apset)
	# cmp.search(apset, apset[1], type=3, cutoff=0.2) 
	# plot(sdfset[names(cmp.search(apset, apset[6], type=2, cutoff=0.4))])
	# db.explain(apset[2])
	# clusters <- cmp.cluster(db=apset[1:100], cutoff = c(0.65, 0.5), save.distances="distmat.rda") # Distance matrix can be loaded with: load("distmat.rda")) 
	# cmp.duplicated(apset, type=2)
