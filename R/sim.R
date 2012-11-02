.onLoad <- function(libname, pkgname) 
{
#pp#    if (!is.null(getOption('disable.chemminer.performance.pack'))
#pp#            && getOption('disable.chemminer.performance.pack') == 1) {
#pp#        cat("ChemmineR Performance Pack is explicitly disabled.\n")
#pp#        options(.use.chemminer.pp = 0)
#pp#    } else if (suppressWarnings(require('ChemmineRpp', quietly=T))) {
#pp#        cat("Using ChemmineR Performance Pack for calculation.",
#pp#        "Set `disable.chemminer.performance.pack' option to 1",
#pp#        "to disable the use of ChemmineR Performance Pack.\n")
#pp#        options(.use.chemminer.pp = 1)
#pp#    }
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

# read a sdf file and returns the atoms with bond list
# argument: the filename, which can also be a connection
# return: the list containing the vector of atoms and the vector of bonds
# this is an internal function and is not intended to be used by users
## ThG: function is not used by new S4 class system.
.parse <- function(file_name, skip=3)
{
    if ("connection" %in% class(file_name))
        con <- file_name
    else
        con <- file(file_name, open="r")
    num_atoms <- 0
    num_bonds <- 0
    # read header line, skip a few lines (given in `skip')
    if (skip > 1)
        readLines(con=con, n=skip)
    line <- readLines(con=con, n=1)
    num_atoms <- as.integer(substr(line, 1, 3))
    num_bonds <- as.integer(substr(line, 4, 6))

    if (.has.pp()) {
        buf <- paste(
            '', '', '',
            line, 
            paste(readLines(con=con, n=num_atoms + num_bonds), collapse='\n'),
            'M END\n$$$$\n',
            sep="\n")
        # close connection if it is opend here
        if (! "connection" %in% class(file_name))
            close(con)
        d <- Descriptors()
        if (Descriptors_parse_sdf(self=d, sdf=buf) == 0) {
            cat("SDF not well-formatted!")
            return(list(n_atoms=0, n_bonds=0, desc_obj=NULL))
        }
        return(list(n_atoms=num_atoms, n_bonds=num_bonds, desc_obj=d))
    }

    if (length(num_atoms) == 0 || num_atoms == 0)
        return(list(atoms=vector(), bonds=list(u=vector(), v=vector()),
            n_atoms=0, n_bonds=0))

    # parse atom block
    temp <- readLines(con=con, n=num_atoms)
    atoms <- gsub(' ', '', substr(temp, 32, 34))

    # parse bond block
    # bonds <- array(0, c(num_bonds, 3))
    temp <- readLines(con=con, n=num_bonds)
    u <- as.integer(substr(temp, 1, 3))
    v <- as.integer(substr(temp, 4, 6))
    t <- as.integer(substr(temp, 7, 9))
    bonds <- list(u=u, v=v, t=t)

# return: a compound object:
#   keys:
#       'atoms': a vector of atoms
#       'bonds': a list indexed by 'u', 'v', and 't'. Each entry will
#       be a vector 

    # close connection if it is opend here
    if (! "connection" %in% class(file_name))
        close(con)

    return(list(atoms=atoms, bonds=bonds, n_atoms=num_atoms, n_bonds=num_bonds))
}

# generates distance matrix for every pair of atoms in compound. Distance is
# the length of shortest path connecting them.
# Returns the distance matrix giving distance between any 2 atoms in the
# compound
# this is an internal function and is not intended to be used by users
.all_pairs_dist <- function(cmp)
{
    s <- length(cmp[['atoms']])
    dis <- array(1024, c(s, s))
    # for every edge, init the distance to 1
    for (i in 1:length(cmp[['bonds']][['u']])) {
        dis[cmp[['bonds']][['u']][i], cmp[['bonds']][['v']][i]] <- 1
        dis[cmp[['bonds']][['v']][i], cmp[['bonds']][['u']][i]] <- 1
    }
    for (k in 1:s) {
        for (i in 1:(s-1)) {
        if (i != k)
            for (j in (i + 1):s) {
            if (j != k)
                if (dis[i, k] + dis[k, j] < dis[i, j]) {
                    dis[i, j] <- dis[i, k] + dis[k, j]
                    dis[j, i] <- dis[i, j]
                }
            }
        }
    }

    return(dis)
}

# generates atom descriptor for atom #i. Descriptor is encoded as a number.
# this is an internal function and is not intended to be used by users
# see also: .interpret_atom
.gen_atom_desc <- function(cmp, i)
{
    neighbors <- 0 # neighbours (non hydrogens)
    pi_electrons <- 0 # number of pi electrons
    num_bonds <- length(cmp[['bonds']][['u']])

    # for every bond, update atom properties
    for (k in 1:num_bonds) {
        if (cmp[['bonds']][['u']][k] == i || cmp[['bonds']][['v']][k] == i) {
            if (cmp[['atoms']][ cmp[['bonds']][['u']][k] ] != 'H' &&
                    cmp[['atoms']][ cmp[['bonds']][['v']][k] ] != 'H') {
                neighbors <- neighbors + 1
                # single to triple bonds
                if (cmp[['bonds']][['t']][k] < 4)
                    pi_electrons <- pi_electrons + cmp[['bonds']][['t']][k] - 1
                # arromatic bonds
                else
                    pi_electrons <- pi_electrons + 0.5
            }
        }
    }

    pi_electrons <- floor(pi_electrons)
    atom_desc <- .elements[[cmp[['atoms']][i]]] * 2^6 + neighbors * 2^3 +
            pi_electrons
    return(atom_desc)
}

# generate atom descriptors for each atom of the compound. 
# this is an internal function and is not intended to be used by users
.gen_all_atom_desc <- function(cmp)
{
    num_atoms <- length(cmp[['atoms']])
    desc <- array()
    for (i in 1:num_atoms) {
        desc[[i]] <- .gen_atom_desc(cmp, i)
    }
    return(desc)
}

# generates atom pair descriptor for a compound. It will encode all atom pairs.
# descriptor = atom descriptors + atom pair information
# this is an internal function and is not intended to be used by users
.gen_atom_pair <- function(cmp)
{
    if (.has.pp() && class(cmp$desc_obj) == "_p_Descriptors") {
        if (is.null(cmp$desc_obj))
            return(vector())
        desc <- c()
        for (i in 1:Descriptors_get_len(self=cmp$desc_obj)) {
            desc <- c(desc,
                    Descriptors_get_descriptor(self=cmp$desc_obj, i=i-1))
        }
        if (length(desc) == 0) return(desc)
        return(.factor_to_vector(as.factor(desc)))
    }

    num_bonds <- length(cmp[['bonds']][['u']])
    if (num_bonds == 0)
        return(c())
    dis <- .all_pairs_dist(cmp)
    num_atoms <- length(cmp[['atoms']])
    desc <- c()
    tail <- 1

    a_desc <- .gen_all_atom_desc(cmp)
    # for every pair of atoms
    for (i in 1:num_atoms) {
        for (j in i:num_atoms) {
            if (i != j && dis[i,j] < 128 && cmp[['atoms']][i] != 'H' &&
                    cmp[['atoms']][j] != 'H') {
                # get the atom descriptor for the two atoms
                desc_i <- a_desc[[i]]
                desc_j <- a_desc[[j]]
                # construct descriptor for the atom pair
                if (desc_i >= desc_j)
                    desc[tail] <- desc_i * 2^20 + dis[i,j] * 2^13 + desc_j
                else
                    desc[tail] <- desc_j * 2^20 + dis[i,j] * 2^13 + desc_i
                tail <- tail + 1
            }
        }
    }
    if (length(desc) == 0)
        return(desc)
    else 
        # remove duplicated atom pair descriptors
        return(.factor_to_vector(as.factor(desc)))
}

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

# generate (atom pair) descriptor from sdf file encoding ONE compound. If the
# file encodes many compounds, only the first one will be parsed. Use
# cmp.parse/sdfs_to_desc to parse multiple compounds.
# returns descriptor for one compound as vector.
sdf_to_desc <- function(filename)
{
    cmp <- .parse(filename)
    desc <- .gen_atom_pair(cmp)
    if (.has.pp() && class(cmp$desc_obj) == "_p_Descriptors") 
        delete_Descriptors(self=cmp$desc_obj)
    desc
}
# alias
cmp.parse1 <- sdf_to_desc

# generate descriptors for all compounds in an sdf file
# returns a list of descriptors of compounds. The list contains two names:
# descdb: the descriptor db, itself being a vector of descriptors.
# names: the names of compounds as a vector.
sdfs_to_desc <- function(filename, quiet=FALSE, type="normal", dbname="")
{
#pp#    if (type == "file-backed") {
#pp#        if (! .has.pp())
#pp#            stop("file-backed parsing is only available with ChemmineR",
#pp#            " Performance Pack package\n")
#pp#        if (! suppressWarnings(require('ChemmineRpp', quietly=T))) 
#pp#            stop('Error: cannot load ChemmineR Performance Pack package\n')
#pp#        if ("connection" %in% class(filename) ||
#pp#            length(grep('(ht|f)tp[s]://', filename)))
#pp#            stop('file-backed parsing can only work with path to a local',
#pp#                ' file.\n')
#pp#        if (dbname == "")
#pp#            stop("file-backed parsing requies `dbname' to be set to an",
#pp#            " output filename without extension.\n")
#pp#    }
    # count compounds first
    con <- file(filename, open="r")
    cat("counting number of compounds in sdf...\n")
    n_compounds <- 0
    sdf_seg <- array()
    cur_line <- 0
    sdf_seg[[1]] <- 0
    read_cid <- TRUE
    cid_buf <- ''
    cids <- array()
    while (TRUE) {
        line <- readLines(con, n=1)
        cur_line <- cur_line + 1
        if (length(line) == 0)
            break;
        if (line == "$$$$") {
            n_compounds <- n_compounds + 1
            sdf_seg[n_compounds + 1] <- cur_line
            read_cid <- TRUE
            cids[n_compounds] <- cid_buf
        } else if (read_cid) {
            cid_buf <- line
            read_cid <- FALSE
        }
    }
    cat(n_compounds, " compounds found\n")
    close(con)

    # real parsing
    descdb <- list()
    if (type != 'file-backed') {
        con <- file(filename, open="r")
        id <- 1
        cur_line <- 0

        for (id in 1:n_compounds) {
            # parse
            cmp <- .parse(con, skip=3)
            # update db
            descdb[[id]] <- .gen_atom_pair(cmp)
            cur_line <- cur_line + 4 + cmp$n_atoms + cmp$n_bonds
            if (! quiet) {
                cid_p <- substr(cids[[id]], 1, 20)
                if (cid_p == cids[[id]]) 
                    msg <- paste("\r", id, "/", n_compounds, " parsed (", 
                        cid_p, ")",
                        " now at line ", cur_line,
                        '                                                     ',
                        sep='')
                else
                    msg <- paste("\r", id, "/", n_compounds, " parsed (", 
                        cid_p, "...)",
                        " just passed line ", cur_line,
                        '                                                     ',
                        sep='')
                cat(substr(msg, 1, 79))
            }
            # proceed to next compound, which might not exist
            temp <- readLines(con, n=sdf_seg[[id + 1]] - cur_line)
            cur_line <- sdf_seg[[id + 1]]
        }
        cat("\nYou can use save(..., file='...', compress=TRUE)",
            "to save the database\n")
        close(con)
    } else {
        dbfile <- paste(dbname, '.cdb', sep='')
        if (batch_parse(filename, dbfile) == 0) {
            stop('Error in parsing using ChemmineR Performance Pack',
                ' package. Check your input file.\n')
        }
        # indexing binary db
        dbf <- file(dbfile, 'rb')
        # skip header
        seek(dbf, 16)
        # read int size
        intsize <- readBin(dbf, integer(), size=1)
        for (id in 1:n_compounds) {
            descdb[[id]] <- c("filedb:", dbname, seek(dbf))
            d_size <- readBin(dbf, integer(), size=intsize)
            seek(dbf, d_size * intsize, origin='current')
        }
        close(dbf)
        cat("\nYour database has been generated and is backed by file",
            dbfile, '\n')
        cat("You can use save(..., file='...', compress=TRUE)", 
            "to save the database\n")
        cat("Also make sure", dbfile, "is not deleted or overwritten", 
            "unless you do not need the database any more.\n")
    }
    return(list(descdb=descdb, cids=cids, sdfsegs=sdf_seg,
                source=filename, type=type))
}
# alias
cmp.parse <- sdfs_to_desc

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

    if (visualize) {
        notes <- as.character(scores[order])
        names(notes) <- rep("similarity", length(notes))
        if (is.null(visualize.query))
            url <- sdf.visualize(db, ids[order], extra=notes,
                                    browse=visualize.browse)
        else
            url <- sdf.visualize(db, ids[order], reference.sdf=visualize.query,
                extra=notes, browse=visualize.browse)
        if (!visualize.browse)
            print(url)
    }
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

# view sdfs in ChemMine
.sdf.visualize.max.files = 100
sdf.visualize <- function(db, cmps, extra=NULL, reference.sdf=NULL,
reference.note=NULL, browse=TRUE, quiet=TRUE)
{
    ## ThG: added for compatability with new S4 classes APset/AP ##
    if(class(db)=="SDF") db <- as(db, "SDFset")
    if(missing(cmps) & class(db)=="SDFset") cmps <- cid(db) # turns cmps into optional argument for SDFset class
    ## ThG: end of lines ##
    if (length(cmps) > 100)
        stop(paste('\n\tYou cannot visualize more than 100 compounds.\n',
            '\tYou supplied ', length(cmps), ' compounds.\n',
            '\tSending too many compounds will take a long time and too much\n',
            '\tresources on the server, or it could crash your browser.\n',
            sep=''))

    # read the reference file if there is any
    if (!is.null(reference.sdf))
        reference.sdf <- .read.one.sdf(reference.sdf)

    if (!is.null(extra) && length(extra) != length(cmps))
        stop(paste('\n\tthe indices and the extra information have",
            " different lengths.\n',
            '\tYou supplied ', length(cmps), ' compounds.\n',
            '\t`extra\' has a length of ', length(extra),
            sep=''))
    
    ## ThG: added for compatability with new S4 classes APset/AP ##
    if(class(db)=="SDFset") {
    	sdfstr <- as(db, "SDFstr")
	sdfs <- paste(paste(unlist(as(sdfstr, "list")), collapse="\n"), "\n", sep="")
	cids <- cmps
    } else {
    	sdfs <- sdf.subset(db, cmps)
    	cids <- db$cids[cmps]
    }
    ## ThG: end of lines ##
    
    # build query
    query <- list(sdf=sdfs, cids=paste(cids, collapse="\1"))
    if (! is.null(extra)) {
        # if an entry in extra is a data frame, format it
        .extra <- c()
        for (i in extra) {
            if (class(i) == 'data.frame' || class(i) == 'matrix')
                .extra <- c(.extra, .data.frame.to.str(i))
            else
                .extra <- c(.extra, paste(as.character(i)))
        }
        query <- c(query, list(extra=paste(.extra, collapse="\1")))
        if (! is.null(names(extra)))
            query <- c(query, list(names=paste(names(extra), collapse="\1")))
    }
    if (! is.null(reference.sdf)) {
        query <- c(query, list(referencesdf=reference.sdf))
        if (! is.null(reference.note)) {
            if (class(reference.note) == 'list' &&
                length(reference.note) == 1) {
                    if (! is.null(names(reference.note)))
                        query <- c(query,
                            list(referencenotename=names(reference.note)))
                    reference.note <- reference.note[[1]]
            }
            if (class(reference.note) == 'data.frame' ||
                class(reference.note) == 'matrix')
                query <- c(query, 
                    list(referencenote=.data.frame.to.str(reference.note)))
            else
                query <- c(query, 
                    list(referencenote=paste(as.character(reference.note))))
        }
    }

    # sending
    if (! quiet) cat('sending SDF to ChemMine\n')
    response <- .postToHost("bioweb.ucr.edu",
        "/ChemMineV2/chemminer/postsdfs", query)
    # reading response
    ref <- gsub(
    '.*\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.(.*)\\.\\.\\.\\.\\.\\.\\.\\.\\.\\..*',
    '\\1', response)
    if (! quiet)
        cat("starting your browser (use `options(browser=\"...\")' to", 
            "customize browser)\n")
    url <- paste(
        c("http://bioweb.ucr.edu/ChemMineV2/chemminer/viewsdfs?", ref),
        collapse='')
    if (browse) browseURL(url)
    return(url)
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

.has.pp <- function()
{
    !is.null(getOption('.use.chemminer.pp')) &&
        getOption('.use.chemminer.pp') != 0
}

########### New Functions ###########

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
                return(as(fpbitma, "FPset"))
        }
}

## Fingerprint comparison and similarity search function 
fpSim <- function(x, y, sorted=TRUE, method="Tanimoto", cutoff=0, top="all", alpha=1, beta=1, ...) {
	## Predefined similarity methods
	if(class(method)=="character") {
	 	if(method=="Tanimoto" | method=="tanimoto") method <- function(a,b,c,d) c/(a+b+c)
	}
	if(class(method)=="character") {
	 	if(method=="Euclidean" | method=="euclidean") method <- function(a,b,c,d) sqrt((c+d)/(a+b+c+d))
	}
	if(class(method)=="character") {
	 	if(method=="Tversky" | method=="tversky") method <- function(a,b,c,d) c/(alpha*a + beta*b + c)
	}
	if(class(method)=="character") {
	 	if(method=="Dice" | method=="dice") method <- function(a,b,c,d) 2*c/(a + c + b + c)
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

######################################
## Query ChemMine Web Tools Service ##
######################################
.serverURL <- "http://chemmine.ucr.edu/ChemmineR/"

# get CIDs from PubChem through ChemMine Web Tools
getIds <- function(cids) {
    if(! class(cids) == "numeric"){
        stop('reference compound ids must be of class \"numeric\"')
    }
	cids <- paste(cids, collapse=",")
	# query server
	response <- postForm(paste(.serverURL, "runapp?app=getIds", sep=""), cids=cids)[[1]]
	if(grepl("^ERROR:", response)){
        stop(response)
    }
    if(grepl("linux", sessionInfo()$platform)) {
        # temporary workaround for linux: save to file with curl and re-open
        # related to R bug 14533 https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14533
        temp <- tempfile()
        command <- paste("curl", response, ">", temp, "2> /dev/null")
        system(command)
        z <- gzcon(file(temp, "rb"))
        sdf <- read.SDFset(read.SDFstr(z))
        close(z)
        unlink(temp)
    } else {
        z <- gzcon(url(response))
        sdf <- read.SDFset(read.SDFstr(z))
        close(z)
    }
	return(sdf)
}

# search PubChem through ChemMine Web Tools with smiles query
searchString <- function(smiles) {
    if(! class(smiles) == "character"){
        stop('reference compound must be a smiles string of class \"character\"')
    } 
    response <- postForm(paste(.serverURL, "runapp?app=searchString", sep=""), smiles=smiles)[[1]]
    if(grepl("^ERROR:", response)){
        stop(response)
    }
	response <- as.numeric(strsplit(response, ",")[[1]])
	return(getIds(response))
}

# search PubChem through ChemMine Web Tools with sdf query
searchSim <- function(sdf) {
    if(! class(sdf) == "SDFset"){
        stop('reference compound must be a compound of class \"SDFset\"')
    } 
    smiles <- sdf2smiles(sdf)
    return(searchString(smiles))
}

# perform sdf to smiles conversion through ChemMine Web Tools
sdf2smiles <- function(sdf) {
    if(! class(sdf) == "SDFset"){
        stop('reference compound must be a compound of class \"SDFset\"')
    } 
	sdf <- sdf2str(sdf[[1]])
	sdf <- paste(sdf, collapse="\n")
	response <- postForm(paste(.serverURL, "runapp?app=sdf2smiles", sep=""), sdf=sdf)[[1]]
	if(grepl("^ERROR:", response)){
        stop(response)
    }
	response <- sub("\n$", "", response) # remove trailing newline
	id <- sub(".*\t(.*)$", "\\1", response) # get id
	response <- sub("\t.*$", "", response) # get smiles
	names(response) <- id
	return(response)
}

# perform smiles to sdf conversion through ChemMine Web Tools
smiles2sdf <- function(smiles) {
    if(! class(smiles) == "character"){
        stop('reference compound must be a smiles string of class \"character\"')
    }
    if(! is.null(names(smiles))){
        smiles <- paste(smiles, names(smiles)[1], sep="\t")
    }
	response <- postForm(paste(.serverURL, "runapp?app=smiles2sdf", sep=""), smiles=smiles)[[1]]
	if(grepl("^ERROR:", response)){
        stop(response)
    }
	response <- strsplit(response, "\n")
	response <- as(as(response, "SDFstr"), "SDFset")
	return(response)
}
