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

## User interface to SDFset() constructor
SDFset <- function(SDFlist=list(), ID=character()) {
	new("SDFset", SDF=SDFlist, ID=ID)
}

setClass("SDFset", representation(SDF="list", ID="character"))


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

## (5.3.4) Ring Perception and Aromaticity Assignment
## Implements with some modifications the exhaustive ring perception algorithm 
## from Hanser et al (1996). URL: http://pubs.acs.org/doi/abs/10.1021/ci960322f

## (a) Iterative removal of atoms with single non hydrogen bonds.
## Returns from molecule only its rings and their inter-connections.
.cyclicCore <- function(x) {
	if(length(x) != 1) stop("x needs to be a single molecule")
	if(class(x) == "SDFset") x <- x[[1]]
	path <- conMA(x, exclude="H")
	noconnect <- rowSums(path) != 0 # Removes atoms with no connections
	path <- path[noconnect, noconnect]
	if(all(dim(path) == 0)) { return(path) } 
	term <- which(rowSums(path > 0)==1)
	while(length(term) > 0) {
		path <- path[-term,-term]
		if(any(dim(path) == 0) | is.vector(path)) { break() }
		term <- which(rowSums(path > 0)==1)
	}
	return(path)
}

## (b) Function to return the longest possible linear bond paths where:
##     - internal atoms have only two heavy atom neighbors
##     - terminal atoms are atoms with more than two heavy atom neighbors or they are ring closures
.linearCon <- function(x) {
	secatoms <- rowSums(x > 0) == 2
	secatoms <- names(which(secatoms))
	.linearCon <- function(vertex, x=x) {
		con <- as.list(names(which(x[vertex, ] > 0)))
		for(i in seq(along=con)) {
			termatom <- con[[i]][length(con[[i]])]
			termcon <- names(which(x[termatom, ] > 0))
			termcon <- termcon[!termcon %in% vertex]
			while(length(termcon) == 1 & !any(duplicated(con[[i]]))) { 
				con[[i]] <- c(con[[i]], termcon)						
				termatom <- con[[i]][length(con[[i]])]
				termcon <- names(which(x[termatom, ] > 0))
				termcon <- termcon[!termcon %in% con[[i]][length(con[[i]])-1]]
			}
		}
		if(paste(sort(unique(con[[1]])), collapse="") == paste(sort(unique(con[[2]])), collapse="")) {
			return(con[[1]])
		} else {
			return(c(rev(con[[1]]), vertex, con[[2]]))
		}
	}
	linearconlist <- lapply(secatoms, function(y) .linearCon(vertex=y, x=x))
	nodups <- !duplicated(sapply(linearconlist, function(y) paste(sort(unique(y)), collapse=""))) 
	linearconlist <- linearconlist[nodups]
	return(linearconlist)
}

## (c) Assemble intermediate results in a list
.update <- function(con, path) {
	## Remove non-terminal atoms in each path
	center_atoms <- unique(unlist(lapply(path, function(x) x[-c(1, length(x))])))
	con1 <- con[!rownames(con) %in% center_atoms, !colnames(con) %in% center_atoms]
	## Add atom pairs with three neighbors to connection list
	if(is.matrix(con1)) {
		remainbonds <- con1
		remainbonds[lower.tri(remainbonds)] <- 0
		remainbonds <- lapply(rownames(remainbonds), function(x) names(which(remainbonds[x,] > 0)))
		names(remainbonds) <- rownames(con1)
		remainbonds <- cbind(rep(names(remainbonds), sapply(remainbonds, length)), unlist(remainbonds, use.names=F))
		remainbonds <- split(remainbonds, seq(along=remainbonds[,1]))
		path <- c(path, remainbonds)
		names(path) <- seq(along=path)
	}
	## Collect complete rings and remove them from path object
	index <- unlist(lapply(path, function(y) any(duplicated(y))))
	rings <- path[index]
	path <- path[!index]
	names(path) <- seq(along=path)
	## Connection list for path component
	conpath <- t(sapply(path, function(x) x[c(1, length(x))]))
	ends <- unique(as.vector(conpath))
	conpath <- lapply(ends, function(x) as.numeric(names(which(rowSums(conpath==x) > 0))))
	names(conpath) <- ends
	conpath <- conpath[sapply(conpath, length) > 1] # removes ends that occur only once 
	## Assemble results in list
	return(list(con=con, conpath=conpath, path=path, rings=rings))
}

## (d) Return rings from cyclist object
.rings <- function(cyclist, upper=Inf) {
	## Define data containers
	pathlist <- cyclist$path
	conpath <- cyclist$conpath
	pathlistnew <- list() 
	rings <- list()
	## Loop to join linear paths/fragments stored in pathlist
	for(i in names(conpath)) {
		if(length(conpath) == 0 | !any(names(conpath) == i)) { next() }
		pos <- t(combn(conpath[[i]], m=2))
		for(j in seq(along=pos[,1])) { 
			p1 <- pathlist[[pos[j,1]]]
			p2 <- pathlist[[pos[j,2]]]
			if(sum(p1[-c(1,length(p1))] %in% p2[-c(1,length(p2))]) > 0) {
				next()
			}
			if(p1[1] == i & p2[1] == i) { # matching s1:s2 
				pathlistnew[[length(pathlistnew)+1]] <- c(rev(p2[-1]), p1)
			}
			if(p1[length(p1)] == i & p2[length(p2)] == i) { # matching e1:e2 
				pathlistnew[[length(pathlistnew)+1]] <- c(p1, rev(p2[-length(p2)]))
			}
			if(p1[1] == i & p2[length(p2)] == i) { # matching s1:e2
				pathlistnew[[length(pathlistnew)+1]] <- c(p2, p1[-1])
			}
			if(p1[length(p1)] == i & p2[1] == i) { # matching e1:s2
				pathlistnew[[length(pathlistnew)+1]] <- c(p1, p2[-1])
			}
		}
		## Various postprocessing routines for joined fragments follow
		if(length(pathlistnew) == 0) { next() }
		## Remove duplicates
		dups <- duplicated(sapply(pathlistnew, function(x) paste(sort(unique(x)), collapse="_")))
		pathlistnew <- pathlistnew[!dups]
		## Set maximum ring size; improves time performance if outer rings are not needed
		if(upper != Inf) {   
			l <- sapply(pathlistnew, length)
			pathlistnew <- pathlistnew[l <= upper]
			if(length(pathlistnew) == 0) { next() }
		} 
		## Collect complete rings and remove them from path object
		index <- unlist(lapply(pathlistnew, function(y) any(duplicated(y[c(1, length(y))]))))
		rings[[length(rings)+1]] <- pathlistnew[index]
		pathlistnew <- pathlistnew[!index]
		## Remove paths with internal duplicates 
		if(length(pathlistnew) > 0) {
			index <- unlist(lapply(pathlistnew, function(y) any(duplicated(y))))
			pathlistnew <- pathlistnew[!index]
		}
		## Update pathlist and conpath
		pathlist <- c(pathlist[-conpath[[i]]], pathlistnew)
		dups <- duplicated(sapply(pathlist, function(x) paste(sort(unique(x)), collapse="_")))
		pathlist <- pathlist[!dups]
		names(pathlist) <- seq(along=pathlist)
		conpath <- t(sapply(pathlist, function(x) x[c(1, length(x))]))
		ends <- unique(as.vector(conpath))
		conpath <- lapply(ends, function(x) as.numeric(names(which(rowSums(conpath==x) > 0))))
		names(conpath) <- ends
		conpath <- conpath[sapply(conpath, length) > 1] # removes ends that occur only once
		pathlistnew <- list()
	}
	## Generate proper output format 
	rings <- unlist(rings, recursive=FALSE)
	dups <- duplicated(sapply(rings, function(x) paste(sort(unique(x)), collapse="_")))
	rings <- c(cyclist$rings, rings[!dups])
	l <- sapply(rings, length) 
	rings <- rings[order(l)]
	if(upper != Inf) { rings <- rings[l <= upper] }
	if(length(rings) > 0) { 
		names(rings) <- paste("ring", seq(along=rings), sep="") 
	} else {
		rings <- NULL
	}
	return(rings)
}

## (e) Identify inner rings
## Expected input x is list of rings from call: 
##      x <- rings(x=sdfset[[1]], upper=Inf, type="all", arom=FALSE) 
.is.inner <- function(x) {
        rnames <- rev(names(x)); names(rnames) <- rnames
        r <- x
        names(r) <- paste(names(r), "_", sep="")
        r <- unlist(r)
        for(i in rnames) {
                tmp <- x[[i]]
                r2 <- r[!gsub("_.*", "", names(r)) %in% i] 
                if(all(tmp %in% r2)) {
                        rnames[i] <- "redundant"
                        r <- r2
                }
        }
        return(x[rev(rnames!="redundant")])
}

## (f) Aromaticity assignment for rings
## Approach
##   (i) Identify rings where all atoms are sp2 hybridized. This means each atom has a 
##       double bond or at least one lone electron pair and is attached to a sp2 hybridized atom
##   (ii) Hueckel's rule needs to be true (4n+2=integer): 2, 6, 10, 14, 18, ... pi electrons 
##        per ring.
.is.arom <- function(sdf, rings) {
	if(length(rings)==0) { return(NULL) } 
	con <- conMA(sdf)
	b2 <- bonds(sdf)
	b <- b2[,"Nbondcount"] - b2[,"charge"]
	names(b) <- paste(b2[,"atom"], "_", seq(along=b2[,1]), sep="") 
	.neighborsFct <- function(x) { sapply(rownames(x), function(y) paste(sort(paste(x[y,][x[y,]!=0], sep="_")), collapse="_")) }
	## Determine aromaticity
	.arom <- function(con, b, r=rings[[1]]) {
		## Identify sp2 hybridized atoms
		sp <- .neighborsFct(x=con)[r]
		sp2 <- grepl("1_2$", sp)
		## Identify atoms with lone electron pairs
		el <- c(Al=3, As=5, At=7, B=3, Bi=5, Br=7, C=4, Cl=7, F=7, Ga=3, Ge=4, I=7, In=3, N=5, O=6, P=5, Pb=4, Po=6, S=6, Sb=5, Se=6, Si=4, Sn=4, Te=6, Tl=3)
		bsub <- b[names(sp)]
		lp <- el[gsub("_.*", "", names(bsub))] - bsub >= 2
		lp[is.na(lp)] <- 0 # If an element is not specified under el, then its NA value in lp is set to zero
		sp2lp <- all(sp2 | lp)
		## Hueckel's rule 
		d <- rowSums(con[r, r]==2)
		double <- sum(d) 
		if(length(d[lp]) != 0) {
			lpcount <- sum(lp[d==0]) * 2
		} else {
			lpcount <- 0
		}
		n <- ((double + lpcount) - 2) / 4
		is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
		hueckel <- is.wholenumber(n)
		return(sp2lp & hueckel)
	}
	return(sapply(names(rings), function(x) .arom(con, b, r=rings[[x]])))
}

## (g) Ring and aromaticity perception by running the functions (1)-(5) 
rings <- function(x, upper=Inf, type="all", arom=FALSE, inner=FALSE) {
        if(!any(c("SDF", "SDFset") %in% class(x))) stop("x needs to be of class SDF or SDFset")
	if(inner==TRUE & upper!=Inf) stop("Inner ring prediction requires upper=Inf")
        runAll <- function(x, upper) {
		con <- .cyclicCore(x)
		if(any(dim(con) == 0) | is.vector(con) | all(con == 0)) { # TRUE for linear compounds 
			if(arom==TRUE) {
				if(type=="all") return(list(RINGS=NULL, AROMATIC=NULL))
				if(type=="count") return(c(RINGS=0, AROMATIC=0))
			} else {
				if(type=="all") return(NULL)
				if(type=="count") return(0)
			}
		} else {
			path <- .linearCon(x=con)
			cyclist <- .update(con=con, path=path)
			myrings <- .rings(cyclist, upper+1) # Plus 'upper+1' is required because at this step all rings have duplicated atoms at ring closure
			myrings <- lapply(myrings, function(x) x[-1]) # Removes duplicated atom at ring closure 
			if(upper==Inf & inner==TRUE) myrings <- .is.inner(x=myrings) # Reduces myrings to inner rings only 
                        if(arom==TRUE) {
				myarom <- .is.arom(sdf=x, rings=myrings)
				if(type=="all") return(list(RINGS=myrings, AROMATIC=myarom))
				if(type=="arom") return(list(AROMATIC_RINGS=myrings[myarom]))
				if(type=="count") return(c(RINGS=length(myrings), AROMATIC=sum(myarom)))
			} else {
				if(type=="all") return(myrings)
				if(type=="count") return(length(myrings))
			}
		}
	}
	if(class(x)=="SDF") { 
		return(runAll(x, upper))
	}
	if(class(x)=="SDFset" & length(x) == 1) { 
		return(runAll(x[[1]], upper))
	}
	if(class(x)=="SDFset" &  length(x) > 1) {
		myrings <- lapply(1:length(x), function(y) runAll(x[[y]], upper))
		names(myrings) <- cid(x)
		if(type=="all" | type=="arom") return(myrings)
		if(type=="count") {
			mycol <- c("RINGS", "AROMATIC")
			return(matrix(unlist(myrings), ncol=length(myrings[[1]]), dimnames=list(names(myrings), mycol[1:length(myrings[[1]])]), byrow=T))
		}
	}
}
## Usage:
# rings(sdfset[1:4], upper=6, type="all", arom=TRUE)
# rings(sdfset[1:4], upper=6, type="arom", arom=TRUE)
# rings(sdfset[1:4], upper=6, type="count", arom=TRUE)
# plot(sdfset[1], print=F, atomnum=T, no_print_atoms="H") 
# plot(sdfset[1:4], print=F, atomnum=T, no_print_atoms="H") 

## (5.3.5) Enumerate Functional Groups
## (a) Generate neighbor information for each heavy atom in a molecule
.neighbors <- function(x, type="countMA") {
        ## Input checks        
        if(!any(c("SDF", "SDFset") %in% class(x))) stop("x needs to be of class SDF or SDFset")
        if(!any(c("all", "count", "countMA") %in% type)) stop("type can only be assigned: all, count or countMA") 
        ## Return neighbors
        .neighborsFct <- function(x, type=type) {
                colnames(x) <- paste(colnames(x), colSums(x>0), sep="_") # Adds number of bonded heavy atom neighbors (non-hydrogens)
                neighbors <- sapply(rownames(x), function(y) paste(sort(paste(gsub("_.*_", "_", colnames(x)[x[y,]!=0]), x[y,][x[y,]!=0], sep="_")), collapse="_"))
                if(type=="all") {
                        return(neighbors)
                }
                if(type=="count" | type=="countMA") {
                        return(table(paste(gsub("_.*", "", names(neighbors)), neighbors, sep=":")))
                }
        }
        ## Run on SDF object
        if(class(x)=="SDF") {
                x <- conMA(x, exclude="H")
                neighbors <- .neighborsFct(x, type)
                return(neighbors)
        }
        ## Run on SDFset objects containing one or many molecules
        if(class(x)=="SDFset") {
                cid <- cid(x)
                x <- conMA(x, exclude="H")
                neighbor_set <- lapply(seq(along=x), function(y) .neighborsFct(x[[y]], type))
                names(neighbor_set) <- cid
                if(type=="all" | type=="count") {
                        return(neighbor_set)
                }
                if(type=="countMA") {
                        columns <- unique(unlist(lapply(seq(along=neighbor_set), function(x) names(neighbor_set[[x]]))))
                        myMA <- matrix(NA, length(neighbor_set), length(columns), dimnames=list(NULL, columns))
                        for(i in seq(along=neighbor_set)) myMA[i, names(neighbor_set[[i]])] <- neighbor_set[[i]]
                        myMA[is.na(myMA)] <- 0
                        rownames(myMA) <- cid
                        return(myMA)
                }
        }
}
## Usage:
# .neighbors(sdfset[1:4], type="all")
# .neighbors(sdfset[1:4], type="countMA")

## (b) Count functional groups
groups <- function(x, groups="fctgroup", type) {
        ## Input checks        
        if(!any(c("SDF", "SDFset") %in% class(x))) stop("x needs to be of class SDF or SDFset")
        if(groups=="fctgroup" & (type=="count" | type=="all")) stop("when groups=\"fctgroup\", only type=\"countMA\" can be used")
        mylength <- length(x)
        ## Support for single molecule objects
	if(mylength==1) { x <- as(x, "SDFset"); y <- x; z <- x; cid(z) <- "dummy"; x <- c(y, z) }
        ## Generate neighbor counts
	neighbors <- .neighbors(x, type)
        ## Return neighbor counts if requested
	if(groups[1]=="neighbors") { 
		if(class(neighbors)=="matrix") {
			neighbors <- neighbors[!rownames(neighbors) %in% "dummy", ]
		} else { 
			neighbors <- neighbors[!names(neighbors) %in% "dummy"]
		}
		return(neighbors) 
	
	}
        ## Count functional groups based on data stored in neighbor matrix
        if(groups[1]=="fctgroup") { 
		groups <- c(RNH2="^C:.*N_1_1$", R2NH="^N:C_._1_C_._1$", R3N="^N:C_._1_C_._1_C_._1$", 
                            ROPO3="^P:O_._._O_._._O_._._O_._.$", ROH="(?=^C:.*O_1_1)(?=^(?:(?!(N|O|P|S)_._(2|3)).)*$)", 
		            RCHO="^O:C_2_2", RCOR="^C:C_._1_C_._1.*O_1_2", RCOOH="^C:.*O_1_1_O_1_2", 
			    RCOOR="^C:.*O_1_2_O_2_1", ROR="^O:.*C_._1_C_._1") 
	}
        groupMA <- sapply(names(groups), function(x) rowSums(neighbors[, rep(grep(groups[x], colnames(neighbors), perl=TRUE),2)]/2))
	## Fix counts for ambiguous functional groups
        if(c("ROR" %in% colnames(groupMA))) groupMA[, "ROR"] <- groupMA[, "ROR"] - groupMA[, "RCOOR"] 
	if(mylength>1) {
                return(groupMA)
        } else {
                return(groupMA[1,])
        }
}
## Usage:
# groups(sdfset[1:20], groups="fctgroup", type="countMA") 
# groups(sdfset[1:4], groups="neighbors", type="countMA")
# groups(sdfset[1:4], groups="neighbors", type="count")
# groups(sdfset[1:4], groups="neighbors", type="all")

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

##########################################################################
## (5.6) Subset SDF/SDFset Objects by Atom Index to Obtain Substructure ##
##########################################################################
## Function to obtain substructure from SDF/SDFset object by providing a row
## index for atom block. Both atom and bond blocks will be substetted accordingly. 
atomsubset <- function(x, atomrows, datablock=FALSE) {
        if(!any(c("SDF", "SDFset") %in% class(x))) stop("x needs to be of class SDF or SDFset")
        if(class(x)=="SDFset" & class(atomrows)!="list") stop("if x is of class SDFset, atomrows argument needs to be a list")
        if(class(x)=="SDFset") {
                if(!all(cid(x) == names(atomrows))) stop("x and atomrows need to have same length and identical component (molecule) names")
        }
        .atomsubset <- function(x, atomrows) {
                hb <- header(x)
                ab <- atomblock(x)[atomrows, ]
                bb <- bondblock(x)
                index <- rowSums(cbind(bb[,1] %in% atomrows, bb[,2] %in% atomrows)) == 2
                bb <- bb[index,]

                ## Update bb to positions in ab
		pos <- as.numeric(gsub(".*_", "", rownames(ab)))
		oripos <- 1:length(pos)
		names(oripos) <- pos
		tmp <- bb[,1:2]; tmp2 <- tmp 
		for(i in oripos) tmp2[tmp==as.numeric(names(oripos[i]))] <- i
		bb[,1:2] <- tmp2
		
		## Outputs
		if(is.vector(bb)) { bb <- t(as.matrix(bb)) }
                if(datablock==FALSE) {
                        sdf <- new("SDF", header=hb, atomblock=ab, bondblock=bb)
                        return(sdf)
                }
                if(datablock==TRUE) {
                        sdf <- new("SDF", header=hb, atomblock=ab, bondblock=as.matrix(bb), datablock=datablock(x))
                        return(sdf)
                }
        }
        if(class(x)=="SDF") {
                return(.atomsubset(x, atomrows))
        }
        if(class(x)=="SDFset") {
                ids <- cid(x)
                sdflist <- lapply(cid(x), function(y) atomsubset(sdfset[[y]], atomrows[[y]]))
                names(sdflist) <- ids
                sdfset <- as(sdflist, "SDFset")
                return(sdfset)
        }
}
## Usage: 
# atomsubset(sdfset[[1]], atomrows=1:18)
# atomsubset(sdfset[1:2], atomrows=list(CMP1=1:18, CMP2=1:12))

#################################
## (5.7) String Search Method ##
#################################
## String search function for SDFset
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
        toplot <- list(atomblock=cbind(atomblock(sdf)[,c(1:2)], as.matrix(bonds(sdf, type="bonds")[,-1])), bondblock=as.matrix(as.data.frame(bondblock(sdf))[,1:3]))
	## Create empty plot with proper dimensions
	plot(toplot[[1]], type="n", axes=F, xlab="", ylab="", ...)
	## Remove C-hydrogens including their bonds 
	if(noHbonds==TRUE) {
		nonbonded <- !1:length(toplot[[1]][,1]) %in% sort(unique(as.vector(toplot[[2]][,1:2])))
                nonbonded <- as.data.frame(toplot[[1]])[nonbonded,]
                CHbondindex <- sapply(seq(toplot[[2]][,1]), function(x) paste(sort(gsub("_.*", "", rownames(toplot[[1]]))[toplot[[2]][x,1:2]]), collapse="") == "CH")
		toplot[[1]] <- toplot[[1]][sort(unique(as.numeric(toplot[[2]][!CHbondindex,1:2]))), ]
		toplot[[2]] <- as.matrix(as.data.frame(toplot[[2]])[!CHbondindex,]) 
                toplot[[1]] <- as.matrix(rbind(toplot[[1]], nonbonded))
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
        ## Add charges 
	charge <- c("0"="", "3"="3+", "2"="2+", "1"="+", "-1"="-", "-2"="2-", "-3"="3-")
        charge <- charge[as.character(labelMA[,"charge"])]
        ## Add hydrogens to non-charged/non-C atoms according to valence rules (some SD files require this)
        Nhydrogens <- c("0"="", "1"="H", "2"="H2", "3"="H3", "4"="H4", "5"="H5", "6"="H6", "7"="H7", "8"="H8") 
        hydrogens <- (labelMA[, "Nbondrule"] + labelMA[, "charge"]) - labelMA[,"Nbondcount"]
        hydrogens[labelMA[,"charge"]!=0] <- 0; hydrogens[hydrogens < 0] <- 0
        hydrogens <- Nhydrogens[as.character(hydrogens)]
        ## Plot data
        if(is.vector(labelMA)) labelMA <- matrix(labelMA, 1, 2, byrow=TRUE, dimnames=list(rownames(toplot[[1]])[!grepl(exclude, rownames(toplot[[1]]))], c("C1", "C2")))
	if(is.matrix(labelMA) & length(labelMA[,1])>=1) {
		atomcol <- gsub("_.*", "", rownames(labelMA)); atomcol[!grepl("N|C|O|H", atomcol)] <- "any"; mycol <- c(C="black", H="black", N="blue", O="red", any="green"); atomcol <- mycol[atomcol]
		## Overplot nodes to display atom labels
		points(x=labelMA[,1], y=labelMA[,2], col="white", pch=16, cex=2.8)
		## Plot atom labels
		if(atomnum==TRUE) {
			text(x=labelMA[,1], y=labelMA[,2], paste(gsub("_", "", rownames(labelMA)), hydrogens, charge, sep=""), cex=atomcex, col=atomcol) 
		} else {
			text(x=labelMA[,1], y=labelMA[,2], paste(gsub("_.*", "", rownames(labelMA)), hydrogens, charge, sep=""), cex=atomcex, col=atomcol)
		}
	}
}
## Usage:
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


