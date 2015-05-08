cstrsplitSym=NA
.onLoad <- function(libname,pkgname) {

	#message("libname: ",libname)
	#message("pkgname: ",pkgname)
	cstrsplitSym <<- getNativeSymbolInfo("cstrsplit")
}

.db.header.size <- 16
# supported elements
.elements <- list(
    R=0,
    H=1,
    He=2,
    Li=3,
    Be=4,
    B=5,
    C=6,
    N=7,
    O=8,
    F=9,
    Ne=10,
    Na=11,
    Mg=12,
    Al=13,
    Si=14,
    P=15,
    S=16,
    Cl=17,
    Ar=18,
    K=19,
    Ca=20,
    Sc=21,
    Ti=22,
    V=23,
    Cr=24,
    Mn=25,
    Fe=26,
    Co=27,
    Ni=28,
    Cu=29,
    Zn=30,
    Ga=31,
    Ge=32,
    As=33,
    Se=34,
    Br=35,
    Kr=36,
    Rb=37,
    Sr=38,
    Y=39,
    Zr=40,
    Nb=41,
    Mo=42,
    Tc=43,
    Ru=44,
    Rh=45,
    Pd=46,
    Ag=47,
    Cd=48,
    In=49,
    Sn=50,
    Sb=51,
    Te=52,
    I=53,
    Xe=54,
    Cs=55,
    Ba=56,
    La=57,
    Ce=58,
    Pr=59,
    Nd=60,
    Pm=61,
    Sm=62,
    Eu=63,
    Gd=64,
    Tb=65,
    Dy=66,
    Ho=67,
    Er=68,
    Tm=69,
    Yb=70,
    Lu=71,
    Hf=72,
    Ta=73,
    W=74,
    Re=75,
    Os=76,
    Ir=77,
    Pt=78,
    Au=79,
    Hg=80,
    Tl=81,
    Pb=82,
    Bi=83,
    Po=84,
    At=85,
    Rn=86,
    Fr=87,
    Ra=88,
    Ac=89,
    Th=90,
    Pa=91,
    U=92,
    Np=93,
    Pu=94,
    Am=95,
    Cm=96,
    Bk=97,
    Cf=98,
    Es=99,
    Fm=100,
    Md=101,
    No=102,
    Lr=103,
    Rf=104,
    Db=105,
    Sg=106,
    Bh=107,
    Hs=108,
    Mt=109,
    Ds=110,
    Rg=111)



# similarity model
# input: descriptors for two compounds. both are vectors
# output: similarity
.cmp.similarity <- function(a, b, mode=1, worst=0,
    dbcon.a=NULL, db.intsize.a=4,
    dbcon.b=NULL, db.intsize.b=4)
{
    ## ThG: added for compatability with new S4 classes APset/AP
    if(class(a)=="APset") { a <- ap(a)[[1]] }
    if(class(b)=="APset") { b <- ap(b)[[1]] }
    if(class(a)=="AP") { a <- ap(a) }
    if(class(b)=="AP") { b <- ap(b) }
    ## ThG: end of lines
    # check argument, make sure they are vectors and not lists
    if (class(a) == 'character' && length(a) == 3 && a[[1]] == 'filedb:') 
        a <- .load.file.backed.desc(a, dbcon.a, db.intsize.a)
    if (class(b) == 'character' && length(b) == 3 && b[[1]] == 'filedb:') 
        b <- .load.file.backed.desc(b, dbcon.b, db.intsize.b)
    if (is.null(b) || is.null(a))
        return(0)
    if (class(b) != "numeric" || class(a) != "numeric") {
        stop("Both arguments must be AP/APset objects or descriptors in the form of vectors.\n",
        "Did you use \"[]\" instead of \"[[]]\" to index the descriptor ",
        "database?")
    }

    if (mode == 1 && worst != 0) {
        # estimate upper bound
        i <- length(a)
        j <- length(b)
        u_bond <- min(i, j) / max(i, j)
        if (u_bond < worst)
            return(0)
    }
    if (mode == 0) {
        # if no mode is given, guess a mode based on compound size
        # difference 
        if (length(a) < .25 * length(b) || length(a) > 4 * length(b))
            mode <- 3
        else
            mode <- 1
    }
    if (mode == 1) 
        return (length(intersect(a, b)) / length(union(a, b)))
    else if (mode == 2)
        return (length(intersect(a, b)) / min(length(a), length(b)))
    else if (mode == 3) {
        s <- length(intersect(a, b)) / min(length(a), length(b))
        return (s^3)
    }
}
cmp.similarity <- function(a, b, mode=1, worst=0)
{
    .cmp.similarity(a, b, mode=mode, worst=worst)
}
cmp.parse1 <- function(filename) sdf2ap(read.SDFset(filename))[[1]]
cmp.parse <- function(filename) apset2descdb(sdf2ap(read.SDFset(filename)))


# search db for similar compounds. `db' gives the database returned by
# `cmp.parse' or `sdfs_to_desc'. `query' gives the descriptor of one compound,
# either returned by `cmp.parse1/sdf_to_desc' or a reference to one instance in
# db such as db$descdb[[123]] two types of cutoff: score cutoff (<=1) or count
# cutoff (>1)
cmp.search <- function(db, query, type=1, cutoff=0.5, return.score=FALSE, quiet=FALSE,
	mode=1, visualize=FALSE, visualize.browse=TRUE, visualize.query=NULL) 
{
    ## ThG: added for compatability with new S4 classes APset/AP
    ## Note: type argument was also added (has no impact on old list object). 
    dbtype <- as.character(class(db))
    if(dbtype=="APset") { db <- apset2descdb(db) }
    if(class(query)=="APset") query <- ap(query[[1]]) 
    if(class(query)=="AP") query <- ap(query) 
    ## ThG: end of lines
    if (.is.file.backed.desc(query))
        query <- .load.file.backed.desc(query)
    dbcon <- NULL
    intsize <- 4
    if (db$type == 'file-backed') {
        for (i in 1:length(db$descdb)) {
            if (.is.file.backed.desc(db$descdb[[i]])) {
                dbcon <- file(paste(db$descdb[[i]][[2]], '.cdb', sep=''), 'rb')
                seek(dbcon, 16)
                intsize <- readBin(dbcon, integer(), n=1, size=1)
                break
            }
        }
    }
    scores <- array()
    ids <- array()
    pos <- 1    # tail of entries
    size <- length(db$cids)
    perstep <- round(size / 25)
    for (i in 1:size) {
        if (!quiet && i %% perstep == 0) {
            steps <- i / perstep
            .progress_bar(paste(min(steps*4, 100), "%", collapse=""))
        }
        cmp <- db$descdb[[i]]
        if (db$type == 'file-backed' && .is.file.backed.desc(db$descdb[[i]]))
            cmp <- .load.file.backed.desc(cmp, dbcon, intsize)
        if (cutoff <= 1)
            score <- .cmp.similarity(cmp, query, mode=mode, worst=cutoff)
        else if (pos - 1 < cutoff)
            score <- .cmp.similarity(cmp, query, mode=mode, worst=0)
        else
            score <- .cmp.similarity(cmp, query, mode=mode,
                    worst=scores[[pos - 1]])
        if (cutoff <= 1) {
            if (score >= cutoff) {
                scores[[pos]] <- score
                ids[[pos]] <- i
                pos <- pos + 1
            }
        } else {
            # using naive insertion sort
            if (pos - 1 < cutoff)
                pos <- pos + 1
            else if (scores[[pos - 1]] >= score)
                next

            pos_ <- pos - 2
            while (pos_ > 0 && scores[[pos_]] < score) {
                scores[[pos_ + 1]] <- scores[[pos_]]
                ids[[pos_ + 1]] <- ids[[pos_]]
                pos_ <- pos_ - 1
            }
            pos_ <- pos_ + 1
            scores[[pos_]] <- score
            ids[[pos_]] <- i
        }
    }
    if (!is.null(dbcon)) close(dbcon)
    if (!quiet) cat("\n")
    if (cutoff <= 1) {
        order <- sort(scores, index=TRUE)[["ix"]]
        order <- order[length(order):1]
    } else
        order <- 1:length(scores)

    ## ThG: modified/added for compatability with new S4 classes APset/AP
    if (return.score & dbtype=="list") {
        return(data.frame(ids=ids[order], scores=scores[order]))
     }
    if (!return.score & dbtype=="list") {
        return(ids[order])
    }
    if (dbtype=="APset") {
	if(type==1) {
		return(ids[order])
	}
	if(type==2) {
		index <- scores[order]
		names(index) <- db$cids[ids[order]]
		return(index)
	}
	if(type==3) {
		return(data.frame(index=ids[order], cid=db$cids[ids[order]], scores=scores[order]))
	}
     }
    ## ThG: end of lines
}

# segment SDFs. given indices of selected compounds, return the concatenation
# of SDFs for them
sdf.subset <- function(db, cmps)
{
    if (!grep('://', db$source) && !file.exists(db$source)) {
        stop("cannot find the source SDF file for this database\n")
    }
    output <- ""
    for (i in cmps)
    {
        sdf <- .extract_single_sdf(db, i)
        output <- paste(output, sdf, sep="")
    }

    return(output)
}

## ThG: function is obsolete for APset class which works with standard R subsetting utilities.
db.subset <- function(db, cmps)
{
    return(list(descdb=db$descdb[cmps], cids=db$cids[cmps],
                sdfsegs=db$sdfsegs[cmps], source=db$source))
}

# helper function for subset_sdf. This function takes an index, and return the
# SDF correponding compound at that index
.extract_single_sdf <- function(db, index)
{
    if (index > length(db$sdfsegs))
        stop("Index out of range in .extract_single_sdf")

    skip <- db$sdfsegs[[index]]
    con <- file(db$source, open="r")
    if (skip > 0) 
        #scan(con, "raw", skip=skip-1, nlines=1, quiet=TRUE)
        readLines(con, skip)
    sdf <- .read.one.sdf(con, collapse=FALSE)

    # replace line 1 with compound index if no name is available
    if (length(grep("[a-zA-z0-9]", sdf[[1]])) == 0)
        sdf[[1]] <- paste("ChemmineR_Unnamed_Compound_", index, sep="")

    # close connection
    close(con)
    return(paste(sdf, sep="", collapse="\n"))
}

# explain an atom pair descriptor
db.explain <- function(desc)
{
    ## ThG: added for compatability with new S4 classes APset/AP ##
    if(class(desc)=="APset") desc <- ap(desc[[1]])
    if(class(desc)=="AP") desc <- ap(desc) 
    ## ThG: end of lines ##
    if ('character' %in% class(desc)) {
        if (.is.file.backed.desc(desc))
            desc <- .load.file.backed.desc(desc)
        else
            stop("The descriptor(s) to be explained is not valid.")
    }
    if (length(desc) == 1)
        return(.explain(desc))
    else {
        result <- c()
        for (i in 1:length(desc))
            result <- c(result, .explain(desc[[i]]))
        return(result)
    }
}

.explain <- function(desc)
{
    # numbers of appearance
    occurence <- desc %% 2^7 + 1
    desc <- desc %/% 2^7
    # extract left atom, distance, and right atom
    left <- desc %/% 2^20
    desc <- desc %% 2^20
    dist <- desc %/% 2^13
    desc <- desc %% 2^13
    right <- desc
    # atom properties
    left <- .interpret_atom(left)
    right <- .interpret_atom(right)

    # output
    return (paste(left$symbol, " [", left$neighbors, " neighbor(s),",
                left$pi_electrons, " pi electrons]<---",
                dist,
                "--->", right$symbol, " [", right$neighbors,
                " neighbor(s),", right$pi_electrons, " pi electrons]",
                collapse="", sep=""))
}

# interpret an atom descriptor
# see also: .gen_atom_desc
.interpret_atom <- function(atom_desc)
{
    # atom
    atom_symbol <- names(.elements[.elements == atom_desc %/% 2^6])
    # neighbors
    neighbors <- atom_desc %% 2^6 %/% 2^3
    # pi electron
    pi_electrons <- atom_desc %% 2^3
    return(list(symbol=atom_symbol, neighbors=neighbors,
                pi_electrons=pi_electrons))
}

# convert a factor to a vector to simplify Union and Intersection operations.
# It basically eliminate repeated names by adding suffix. For example, if `a'
# appears twice in the factor, it will be converted to `a1' and `a2'. 
# this is an internal function and is not intended to be used by users
.factor_to_vector <- function(fac)
{
    vec <- vector(mode='numeric')
    index <- 1
    freq <- table(fac)
    keys <- names(freq)
    for (i in 1:length(keys)) {
        key <- keys[[i]]
        if (freq[[key]] == 1) {
            vec[[index]] <- as.numeric(key) * 2^7
            index <- index + 1
        } else {
            for (j in 1:freq[[key]]) {
                vec[[index]] <- as.numeric(key) * 2^7 + j - 1
                index <- index + 1
            }
        }
    }

    return(vec)
}

.uniquifyAtomPairs <- function(desc) {
	.Call("uniquifyAtomPairs",desc)
}

.progress_bar_int_cnt <- 0
.progress_bar <- function(label=NULL) {
    progress <- c("|", "/", "-", "\\")
    unlockBinding(".progress_bar_int_cnt", environment(.progress_bar))
    .progress_bar_int_cnt <<- .progress_bar_int_cnt %% 4 + 1
    if (is.null(label)) {
        cat("\r", progress[[.progress_bar_int_cnt]])
    } else {
        cat("\r                                           \r")
        cat(progress[[.progress_bar_int_cnt]], label)
    }
}

.postToHost <- function(host, path, data.to.send, port=80)
{
    if(missing(path)) path <- "/"
    ua <- "User-Agent: Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.8.1.1)"


    boundary <- paste(
            as.integer(runif(20, min=1, max=10)), collapse=""
    )

    header <- NULL
    header <- c(header, paste("POST ", path, " HTTP/1.1", sep=""))
    header <- c(header, paste("Host: ", host, sep=""))
    header <- c(header, "Connection: close")
    header <- c(header, ua)
    header <- c(header, "Accept: */*")
    header <- c(header,
            paste("Content-type: multipart/form-data; boundary=",
                boundary, sep="")
    )
    boundary <- paste('--', boundary, sep='')

    mcontent <- NULL # keeps the content.
    size <- 0

    for (i in 1:length(data.to.send)) {
        value <- data.to.send[[i]]
        key <- names(data.to.send)[i]

        content <- paste(
            boundary, "\n",
            "Content-Disposition: form-data; name=\"",
            as.character(key), "\"", "\n",
            "\n",
            as.character(value), "\n",
            sep=""
        )
        mcontent <- c(mcontent, charToRaw(content))
    }

    size <- length(mcontent) + length(strsplit(boundary, "")[[1]]) + 4
    header <- c(header, paste("Content-length: ", size, "\n\n", sep=""))

    postdata <- c(
        charToRaw(paste(header, collapse="\n")),
        mcontent,
        charToRaw(paste(boundary,"--\n\n",sep=""))
    )

    rm(header,mcontent)

    scon <- socketConnection(host=host,port=port,open="a+b",blocking=TRUE)
    writeBin(postdata, scon, size=1)

    output <- character(0)
    repeat{
        ss <- rawToChar(readBin(scon, "raw", 2048))
        output <- paste(output,ss,sep="")
        if(regexpr("\r\n0\r\n\r\n",ss)>-1) break()
        if(ss == "") break()
    }
    close(scon)
    return(output)
}

.read.one.sdf <- function(filename, collapse=TRUE)
{
    if (length(filename) != 1)
        stop('reference sdf must be a vector of only one entry!')

    if (!("connection" %in% class(filename)) &&
        !("character" %in% class(filename)))
        stop('reference sdf must be an SDF file or one character string of sdf')

    if ("connection" %in% class(filename))
        con <- filename
    # is filename really a filename a the SDF content? check existence of \n
    else if (length(grep("\n", filename)) == 0 &&
        length(grep("\r", filename)) == 0)
        con <- file(filename, open="r")
    else
        con <- textConnection(filename)

    result <- c()
    while (TRUE) {
        line <- readLines(con=con, n=1)
        if (length(line) == 0 || "$$$$" %in% line) {
            if (length(result) == 0)
                stop("empty SDF!")
            result <- c(result, '$$$$\n')
        if (collapse)
            return (paste(result, collapse="\n"))
        else
            return (result)
        } else
            result <- c(result, line)
    }
}

.data.frame.to.str <- function(x)
{
    con <- textConnection("string", "w")
    write.table(format(x), con, row.names=F, quote=F)
    close(con)
    return(paste(string, collapse='\n'))
} 

.is.file.backed.desc <- function(desc)
{
    "character" %in% class(desc) && length(desc) == 3 && desc[[1]] == 'filedb:'
}

.load.file.backed.desc <- function(desc, dbcon=NULL, intsize=4)
{
    dbname <- desc[[2]]
    pos <- desc[[3]]
    closeme <- FALSE
    if (is.null(dbcon)) {
        db_a <- paste(dbname, '.cdb', sep='')
        if (!file.exists(db_a)) stop("cannot find database file ", db_a)
        dbcon <- file(db_a, 'rb')
        closeme <- TRUE
        seek(dbcon, .db.header.size)
        intsize <- readBin(dbcon, integer(), n=1, size=1)
    }
    seek(dbcon, as.integer(pos))
    size <- readBin(dbcon, integer(), n=1, size=intsize)
    desc <- readBin(dbcon, integer(), n=size, size=intsize)
    if (closeme) close(dbcon)
    .factor_to_vector(as.factor(desc))
}

.haveOB <- function()
{
	if(suppressWarnings(require('ChemmineOB', quietly=T))) {
		return(TRUE)
   }else{
		return(FALSE)
	}

}
.ensureOB <- function(mesg = paste("ChemmineOB is required to make use of this function.",
										 "This package can be installed from BioConductor with the ",
										 "command 'biocLite(\"ChemmineOB\"). ",
										 "It is not currently available for windows however.",
										 "See http://bioconductor.org/packages/devel/bioc/html/ChemmineOB.html",
										 "for more information"))
{
	if(!.haveOB())
		stop(mesg)
}


##############################################
## PubChem Fingerprint Similarity Searching ##
##############################################

## Convert base 64 encoded PubChem fingerprints to binary matrix or string vector
## Definition of PubChem's substructure dictionary-based fingerprints:
## ftp://ftp.ncbi.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
## For debugging use command: as.integer(rev(intToBits("19")[1:6])) 

## Load PubChem's substructure dictionary from data dir of library 
data(pubchemFPencoding); pubchemFPencoding <- pubchemFPencoding 
fp2bit <- function(x, type=3, fptag="PUBCHEM_CACTVS_SUBSKEYS") {
	## Covert base 64 strings to matrix
	if(class(x)=="SDFset") {
		blockmatrix <- datablock2ma(datablocklist=datablock(x))
		fp <- blockmatrix[, fptag]
	}
	if(class(x)=="matrix") { 
		blockmatrix <- x
		fp <- blockmatrix[, fptag]
	}
	if(class(x)=="character") { 
		fp <- x
	}
	fpma <- unlist(strsplit(fp, ""))
	fpma <- matrix(fpma, length(fp), nchar(fp[1]), byrow=TRUE)
	fpma <- fpma[, 1:154] # remove padding signs '='

	## base 64 decoding (base 64 alphabet from http://www.faqs.org/rfcs/rfc3548.html)
	base64 <- c(A=0,B=1,C=2,D=3,E=4,F=5,G=6,H=7,I=8,J=9,K=10,L=11,M=12,N=13,O=14,P=15,
                    Q=16,R=17,S=18,T=19,U=20,V=21,W=22,X=23,Y=24,Z=25,a=26,b=27,c=28,d=29,
                    e=30,f=31,g=32,h=33,i=34,j=35,k=36,l=37,m=38,n=39,o=40,p=41,q=42,r=43,
                    s=44,t=45,u=46,v=47,w=48,x=49,y=50,z=51,"0"=52,"1"=53,"2"=54,"3"=55,
                    "4"=56,"5"=57,"6"=58,"7"=59,"8"=60,"9"=61,"+"=62,"/"=63)
	fpbitma <- as.integer(intToBits(base64[as.vector(t(fpma))]))
	fpbitma <- matrix(fpbitma, length(fpma[,1])*154, 32, byrow=TRUE)[,6:1]
	fpbitma <- matrix(t(fpbitma), length(fpma[,1]), 6*154, byrow=TRUE)[,33:913]
        pubchemFP <- pubchemFPencoding[,2]
	names(pubchemFP) <- pubchemFPencoding[,1]
	colnames(fpbitma) <- pubchemFP
	rownames(fpbitma) <- names(fp)
	if(type==1) {
                return(apply(fpbitma, 1, paste, collapse=""))  
        }
	if(type==2) {
                return(fpbitma)
        }
	if(type==3) {
                #return(as(fpbitma, "FPset"))
                return(new("FPset",fpma=fpbitma,type="pubchem"))
        }
}

## Fingerprint comparison and similarity search function 
fpSimOrig <- function(x, y, sorted=TRUE, method="Tanimoto", addone=1, cutoff=0, top="all", alpha=1, beta=1, ...) {
	## Predefined similarity methods
	if(class(method)=="character") {
	 	if(method=="Tanimoto" | method=="tanimoto") method <- function(a,b,c,d) (c+addone)/(a+b+c+addone)
	}
	if(class(method)=="character") {
	 	if(method=="Euclidean" | method=="euclidean") method <- function(a,b,c,d) sqrt((c+d+addone)/(a+b+c+d+addone))
	}
	if(class(method)=="character") {
	 	if(method=="Tversky" | method=="tversky") method <- function(a,b,c,d) (c+addone)/(alpha*a + beta*b+c+addone)
	}
	if(class(method)=="character") {
	 	if(method=="Dice" | method=="dice") method <- function(a,b,c,d) (2*c+addone)/(a+c+b+c+addone)
	}
	## Check for valid inputs
	if(!any(c(is.vector(x), class(x)=="FP", class(x)=="FPset" & length(x)==1))) stop("x needs to be object of class FP, FPset of length one, or vector")
        if(!any(c(is.vector(y), is.matrix(y), class(y)=="FP", class(y)=="FPset"))) stop("y needs to be object of class FP/FPset, vector or matrix")
	## Convert FP/FPset inputs into vector/matrix format
	if(class(x)=="FP") x <- as.numeric(x)
	if(class(x)=="FPset") x <- as.numeric(x[[1]])
	if(class(y)=="FP") y <- as.numeric(y)
	if(class(y)=="FPset") y <- as.matrix(y)
        ## Computate similarities
	if(class(y)=="matrix") {
		c <- colSums((t(y) + x) == 2)
		d <- colSums((t(y) + x) == 0)
		b <- rowSums(y) - c
	} else {
		c <- sum(x + y == 2)
		d <- sum(x + y == 0)
		b <- sum(y) - c
	}
	a <- sum(x) - c
	if(sorted==TRUE) {
                tmp <- rev(sort(method(a,b,c,d)))
                if(top!="all") tmp <- tmp[1:top]
		if(top!=1) { 
			return(tmp[tmp>=cutoff])
		} else {
			return(tmp) # returns always a value when top=1; required for compatibility with cmp.cluster
		}
	}
	if(sorted==FALSE) {
		tmp <- method(a,b,c,d)
		if(top!=1) {
                	return(tmp[tmp>=cutoff])
		} else {
			return(tmp) # returns always a value when top=1; required for compatibility with cmp.cluster
		}
	}
}

fpSim <- function(x, y, sorted=TRUE, method="Tanimoto", 
						addone=1, cutoff=0, top="all", alpha=1, beta=1,
						parameters=NULL,scoreType="similarity") {

	if(class(method)=="character") {
	 	if(method=="Tanimoto" | method=="tanimoto") 
			method <- 0
		else if(method=="Euclidean" | method=="euclidean") 
			method <- 1
		else if(method=="Tversky" | method=="tversky") 
			method <- 2
		else if(method=="Dice" | method=="dice") 
			method <- 3
		else 
			stop("invalid method found: ",method)
	}else
		return(fpSimOrig(x,y,sorted=sorted,method=method,addone=addone,cutoff=cutoff,top=top,
							  alpha=alpha, beta=beta))
		#stop("invalid method type found: ",class(method))


	if(!any(c(is.vector(x), class(x)=="FP", class(x)=="FPset" & length(x)==1))) 
		stop("x needs to be object of class FP, FPset of length one, or vector")
   if(!any(c(is.vector(y), is.matrix(y), class(y)=="FP", class(y)=="FPset"))) 
		stop("y needs to be object of class FP/FPset, vector or matrix")

	## Convert FP/FPset inputs into vector/matrix format
	if(class(x)=="FP") 
		x <- as.numeric(x)
	else if(class(x)=="FPset") 
		x <- as.numeric(x[[1]])

	if(class(y)=="FP") {
		dbSize = 1
		y <- as.numeric(y)
	}
	else if(class(y)=="FPset"){
		 dbSize=length(y)
		 y <- as.matrix(y)
	}
	else if(is.vector(y))
		dbSize=1
	else if(is.matrix(y))
		dbSize=nrow(y)

	if(is.vector(y)) y <- t(as.matrix(y))
   

	result=.Call("similarity",x,y,method,addone,alpha,beta)
	names(result) = rownames(y)

	if(!is.null(parameters)){
		numBitsSet = sum(x)+1

		#N = parameters$count[numBitsSet]
		#if(N==0){ # no stats collected for this number of bits so use global values
		if(parameters$count[numBitsSet]==0){ # no stats collected for this number of bits so use global values
			warning("no parameters avaliable for fingerprints with ",numBitsSet-1," bits set, using global parameters")
			numBitsSet=nrow(parameters) #global stats are last element of parameters
			#N = parameters$count[numBitsSet]
		}

		#message("using stats for ",numBitsSet-1," bits")
		#print(parameters[numBitsSet,])
		
		avg = parameters$avg[numBitsSet]
		varience = parameters$varience[numBitsSet]
		alpha = parameters$alpha[numBitsSet]
		beta = parameters$beta[numBitsSet]

		evalues = dbSize*(1-pbeta(result,alpha,beta))
		scores <<- data.frame(similarity=result,
									 zscore=(result - avg) /sqrt(varience),
									 evalue=evalues,
									 pvalue=1-exp(-evalues))
		
		titles = 1:4
		names(titles)=colnames(scores)
		if(sorted){
			if(titles[scoreType] <=2 ) #similarity or zscore, bigger is better
				scores = scores[order(-scores[,titles[scoreType]]),]
			else
				scores = scores[order(scores[,titles[scoreType]]),]
			if(top!="all")
				scores = scores[1:top,]
		}

		if(!missing(cutoff)){
			if(titles[scoreType] <=2 ) #similarity or zscore, bigger is better
				cutScores = scores[scores[,titles[scoreType]]>= cutoff,]
			else   # e or p values, lower is better
				cutScores = scores[scores[,titles[scoreType]]<= cutoff,]

			# make sure we don't lose all results do to cutoff
			if( nrow(cutScores) >= 1)
				scores = cutScores
		}
		scores
	}
	else{

		if(sorted) {
			result = sort(result, decreasing=TRUE)
			if(top!="all")
				result = result[1:top]
		}
		cutResult = result[result >= cutoff]
		if( length(cutResult) >= 1)
			cutResult
		else # make sure we don't lose all results do to cutoff
			result
	}
	

}
######################################
## Query ChemMine Web Tools Service ##
######################################
.serverURL <- "http://chemmine.ucr.edu/ChemmineR/"

# perform sdf to smiles conversion
sdf2smilesOB <- function(sdf) {
    if(! class(sdf) == "SDFset"){
        stop('reference compound must be a compound of class \"SDFset\"')
    } 

	 if(.haveOB()){
		 #sdfstrList=as(as(sdf,"SDFstr"),"list")
		 #defs = paste(Map(function(x) paste(x,collapse="\n"), sdfstrList),collapse="\n" )
		 defs = sdfSet2definition(sdf)
		 t=Reduce(rbind,strsplit(unlist(strsplit(convertFormat("SDF","SMI",defs),
															  "\n",fixed=TRUE)),
					 "\t",fixed=TRUE))
		 if(class(t)=="character"){ # R rearranged our matrix because there was only one result
			 smiles=t[1]
			 names(smiles)=t[2]
		 }else{
			 smiles = t[,1]
			 names(smiles)= t[,2]
		 }
		 as(smiles, "SMIset")
	 }else{
		 message("ChemmineOB not found, falling back to web service version. This will be much slower")
		 sdf2smilesWeb(sdf)
	 }
}
sdf2smiles <- sdf2smilesOB

sdf2smilesWeb <- function(sdfset,limit=100){
	#message("class of sdfset: ",class(sdfset))
	 if(length(sdfset) > limit)
		 sdfset = sdfset[1:limit]

	 smiles =c()
	 for(i in seq(along=sdfset)){

		#message("class of sdfset[[]]: ",class(sdfset[[i]]))
		 sdf <- sdf2str(sdfset[[i]])
		 sdf <- paste(sdf, collapse="\n")
		 response <- postForm(paste(.serverURL, "runapp?app=sdf2smiles", sep=""), sdf=sdf)[[1]]
		 if(grepl("^ERROR:", response)){
				  stop(response)
		 }
		 response <- sub("\n$", "", response) # remove trailing newline
		 id <- sub(".*\t(.*)$", "\\1", response) # get id
		 response <- sub("\t.*$", "", response) # get smiles
		 names(response) <- id
		 smiles = c(smiles,response)
	 }
	 as(smiles,"SMIset")
}

smiles2sdfOB <- function(smiles) {
    if(!class(smiles) %in% c("character", "SMIset")){
        stop('input must be SMILES strings stored as \"SMIset\" or \"character\" object')
    }
	 if(class(smiles)=="SMIset") smiles <- as.character(smiles)
	 if(.haveOB()) {
		 sdf = definition2SDFset(convertFormat("SMI","SDF",paste(paste(smiles,names(smiles),sep="\t"), collapse="\n")))
		 cid(sdf)=sdfid(sdf)
		 sdf
	 }else{
		 message("ChemmineOB not found, falling back to web service version. This will be much slower")
		 sdf =smiles2sdfWeb(smiles)
		 cid(sdf)=sdfid(sdf)
		 sdf
	 }
}
smiles2sdf <- smiles2sdfOB

regenCoordsOB <- function(sdf){
	applyOptions(sdf,data.frame(names="gen2D",args=""))
}
regenerateCoords <- regenCoordsOB

generate3DCoordsOB <- function(sdf){
	applyOptions(sdf,data.frame(names="gen3D",args=""))
}
generate3DCoords <- generate3DCoordsOB

canonicalizeOB <- function(sdf){
	applyOptions(sdf,data.frame(names="canonical",args=""))
}
canonicalize <- canonicalizeOB

canonicalNumberingOB <- function(sdf){
	.ensureOB()
	results=canonicalNumbering_OB(obmol(sdf))
	names(results) = cid(sdf)
	results
}
canonicalNumbering <- canonicalNumberingOB

applyOptions <- function(sdf,options){
	.ensureOB()

	if(class(sdf) == "SDFset" || class(sdf)=="SDF"){
		defs = sdfSet2definition(sdf)
		sdfNew = definition2SDFset(convertFormat("SDF","SDF",defs,options))
		cid(sdfNew)=sdfid(sdf)
		if(class(sdf)=="SDF")
			sdfNew[[1]]
		else
			sdfNew
	}else
		stop("input to applyOptions must be a SDFset or SDF object. Found ",class(sdf))
}

# perform smiles to sdf conversion through ChemMine Web Tools
smiles2sdfWeb <- function(smiles,limit=100) {
	if(length(smiles) > limit)
		smiles = smiles[1:limit]

	 smileStrings = if(! is.null(names(smiles)))
        paste(smiles, names(smiles), sep="\t")
	 else
		  smiles
    
	 sdfs = c()
	 for(smile  in smileStrings){
		 response <- postForm(paste(.serverURL, "runapp?app=smiles2sdf", sep=""), smiles=smile)[[1]]
		 if(grepl("^ERROR:", response))
			  stop(response)
		 
		 response <- strsplit(response, "\n")
		 response <- as(as(response, "SDFstr"), "SDFset")[[1]]
		 #response <- as(response, "SDFstr")
		 sdfs = c(sdfs,response)
	 }
	 names(sdfs)=names(smiles)
	 as(sdfs,"SDFset")
}

times = new.env()
times$descT = 0
times$uniqueT = 0
genAPDescriptors <- function(sdf,uniquePairs=TRUE){

	# t=Sys.time()
	d=.Call("genAPDescriptor",sdf)
	# times$descT <- times$descT + (Sys.time() - t)

	if(uniquePairs){
		# t=Sys.time()
		d= .uniquifyAtomPairs(d)
		# itimes$uniqueT <- times$uniqueT + (Sys.time() - t)
	}

	d
}

propOB <- function(sdfSet){
	.ensureOB()
	results = prop_OB(obmol(sdfSet))
	rownames(results) = cid(sdfSet)
	results
}

fingerprintOB <- function(sdfSet,fingerprintName){
	.ensureOB()
   x = fingerprint_OB(obmol(sdfSet),fingerprintName)
	if(is.vector(x)) x= t(as.matrix(x))

	fpset = new("FPset",fpma=x,
					type=fingerprintName)
	cid(fpset) = cid(sdfSet)
	fpset
}
smartsSearchOB <- function(sdfset,smartsPattern,uniqueMatches=TRUE){
	.ensureOB()
	smartsSearch_OB(obmol(sdfset),smartsPattern,uniqueMatches)
}
exactMassOB <- function(sdfset){
	.ensureOB()
	exactMass_OB(obmol(sdfset))
}


#compounds should be items that can be passed into similarity
maximallyDissimilar <- function(compounds,n,similarity = cmp.similarity) {
	dist = function(a,b) 1-similarity(a,b)
	if(length(compounds) <= 2)
		stop("maximallyDissimilar: compounds must have more than 2 elements")

	#pick first selection at random
	selected = sample(1:length(compounds),1)

	#compute dist between compounds[-selected] and compounds[selected]
	minDistsToSelected = sapply(1:length(compounds),
										 function(i) dist(compounds[[i]],compounds[[selected]]))

	for(k in 2:n){
		#find candidate farthest from those selected
		newPoint = which.max(minDistsToSelected)
		selected = c(selected,newPoint)
		minDistsToSelected[newPoint] = 0

		#update min distances
		minDistsToSelected = sapply(1:length(compounds),
											 function(i) 
												 if(minDistsToSelected[i]==0) 0 
												 else min(minDistsToSelected[i],dist(compounds[i],compounds[newPoint])))
	}
	
	compounds[selected]
}

sdf2OBMol<- function(sdfSet){
	.ensureOB()
	defs = paste(Map(function(x) paste(x,collapse="\n"),
						  as(as(sdfSet,"SDFstr"),"list")),"\n",
					 sep="",collapse="" )

	forEachMol("SDF",defs,identity,c)
}

