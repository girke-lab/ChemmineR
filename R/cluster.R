# db: the compound descriptor db. generated with cmp.parse(..)
# cutoff: cutoff used to build the clusters. Cutoff is the distance cutoff, not
# similarity cutoff. It can be a vector
# is.similarity: set when the cutoff is a similarity rather than a distance
# cutoff
# save.distance: whether to compute and save the distance matrix before
# clustering. It will prompt for a filename to save the distance matrix for
# future use.
# use.distance: the pre-computed distance matrix to use, if you have any.
# RETURN: a matrix giving the cluster ID for each compound in the db under
# different cutoffs. Each row corresponds to one compound, and each column
# corresponds to a cutoff

cmp.cluster <- function(db, cutoff, is.similarity=TRUE, save.distances=FALSE,
        use.distances=NULL, quiet=FALSE, ...)
{
    ## ThG: added for compatibility with new S4 classes APset/AP
    dbtype <- as.character(class(db))
    if(dbtype=="APset") { db <- apset2descdb(db) }
    if(dbtype=="FPset") { db <- list(descdb=db, cids=cid(db), sdfsegs=NULL, source="FPset", type="FPset") }
    ## ThG: end of lines
    
    # see if db if file-backed
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
    # prepare column names from cutoffs
    colname_f <- function(x) {paste("CLID_", x, collapse="", sep="")}
    colname <- apply(as.array(cutoff), 1, colname_f)
    o_colname <- colname

    # prepare cutoffs for real computation
    if (is.similarity)
        cutoff <- 1 - cutoff

    if (save.distances != FALSE) {
        distmat <- .calc.distmat(db$descdb, quiet=quiet, 
            dbcon.a=dbcon, dbcon.b=dbcon,
            db.intsize.a=intsize, db.intsize.b=intsize,
        ...)
        if ('character' %in% class(save.distances)) {
            save(distmat, file=save.distances)
        }
        distf <- function(i, j) {
            return(distmat[i, j])
        }
    }
    else if (! is.null(use.distances)) {
        distf <- function(i, j) {
            return(use.distances[i, j])
        }
    } 
    ## ThG: added to make function easier to use with new S4 classes APset/AP
    else if (class(db$descdb) == "FPset") {
        distf <- function(i, j) {
            return(1-fpSim(db$descdb[[i]], db$descdb[[j]], top=1, ...))
        }
    ## ThG: end of lines
    } else {
        distf <- function(i,j) {
            args <- list(...)
            if (! is.null(args$worst) || ! is.null(args$mode))
                return(1 - .cmp.similarity(db$descdb[[i]], db$descdb[[j]], 
                        dbcon.a=dbcon, dbcon.b=dbcon,
                        db.intsize.a=intsize, db.intsize.b=intsize,
                        ...))
            else
                return(1 - .cmp.similarity(db$descdb[[i]], db$descdb[[j]], 
                        dbcon.a=dbcon, dbcon.b=dbcon,
                        db.intsize.a=intsize, db.intsize.b=intsize,
                        mode=1, worst=1 - cutoff[[1]], ...))
        }
    }
    if (!is.null(dbcon)) close(dbcon)

    # sort cutoffs decreasingly
    o <- order(cutoff, decreasing=TRUE)
    cutoff <- cutoff[o]
    colname <- colname[o]
    n_cutoff <- length(cutoff)
    # record of current clustering result
    cluster_id <- matrix(0, ncol=n_cutoff, nrow=length(db$descdb))
    colnames(cluster_id) <- colname

    # all nodes to label/cluster
    all <- 1:length(db$descdb)
    # for progress bar
    perstep <- round(length(db$descdb) / 100) + 1
    for (i in all) {
        # has i already been clustered?
        if (cluster_id[[i, 1]] != 0)
            next
        # else, use i as the lead to get all neighbors
        cluster_id <- .cluster_g(i, all, distf, cutoff, cluster_id, quiet=quiet)
        # progress bar
        ratio <- floor(sum(cluster_id[, 1] != 0) / perstep)
        if (! quiet) .progress_bar(paste(min(ratio, 100), "%", collapse=""))
    }
    cat("\n")

    cluster_id <- as.data.frame(cluster_id)
    cluster_id$ids <- 1:length(cluster_id[,1])
    # convert each column to factor; compute cluster size
    for(i in 1:n_cutoff) {
        col_label <- colname[[i]]
        # size
        clsize <- as.data.frame(table(cluster_id[,col_label]))
        # add column of cluster size
        cluster_id <- merge(cluster_id, clsize, by.x=col_label, by.y=1, all.x=T)
        names(cluster_id)[length(cluster_id)] <- gsub('CLID', 'CLSZ', col_label)
    }

    cat("sorting result...")
    # reorder columns based on user-supplied order
    all_col_labels = array("ids")
    for (i in 1:n_cutoff) {
        all_col_labels[[2 * i]] <- gsub('CLID', 'CLSZ', o_colname[[i]])
        all_col_labels[[2 * i + 1]] <- o_colname[[i]]
    }
    cluster_id <- cluster_id[,all_col_labels]
    # sort by cluster size corresponding to the first clustering cutoff
    if (n_cutoff == 1)
        cluster_id <- cluster_id[order(-cluster_id[,2], cluster_id[,3],
                cluster_id[,1]),]
    else
        cluster_id <- cluster_id[order(-cluster_id[,2], cluster_id[,3],
                -cluster_id[,4], cluster_id[,5], cluster_id[,1]),]

    cat("\n")
    rownames(cluster_id) <- cluster_id[,1]
    ## ThG: added to make function easier to use with new S4 classes APset/AP
    if(dbtype=="APset") { cluster_id[,"ids"] <- db$cids[cluster_id[,"ids"]] }
    if(dbtype=="FPset") { cluster_id[,"ids"] <- db$cids[cluster_id[,"ids"]] }
    ## ThG: end of lines
    return(cluster_id)
}

.calc.distmat <- function(descdb, quiet=FALSE, ...) {
	if(class(descdb)=="FPset") {
    		return(1 - sapply(cid(descdb), function(x) fpSim(descdb[x], descdb, sorted=FALSE)))
	} else { 
  		if (!quiet) cat("calculating distance matrix\n")
		len = length(descdb)
		distmat <- matrix(1, ncol=len, nrow=len)
		for (i in 1:(len-1)) {
			distmat[i, i] <- 0
			for (j in (i+1):len) {
				d <- 1 - .cmp.similarity(descdb[[i]], descdb[[j]], ...)
				distmat[i, j] <- d
				distmat[j, i] <- d
        		}
        		prog_ratio <- i / (len - 1)
        		prog_ratio <- prog_ratio * 2 - prog_ratio * prog_ratio
        		if (! quiet)
            		.progress_bar(paste(min(prog_ratio * 100, 100), "%", collapse=""))
		}
	if (!quiet) cat("distance matrix is successfully generated\n")
	return(distmat)
	}
}

# will only consider the first cluster cutoff, if multiple cutoffs are used
# generating the cluster
cluster.sizestat <- function(cls, cluster.result=1)
{
   st <- data.frame(table(factor(cls[,cluster.result * 2])))
   # count clusters of each size
   st[,2] <- st[,2] / as.numeric(as.vector(st[,1]))
   names(st) <- c("cluster size", "count")
   return(st)
}

# visualize a cluster. db is the database generated by cmp.parse. cls is the 
# cluster generated by cluster, and size_cutoff the cutoff size for the
# clusters to be considered - clusters with smaller size will be ignored in the
# visualization. The distance matrix will be calculated on fly for compounds in
# the cluster under concern, but you can provide the (full) distance matrix in
# distmat, which must cover all compounds in db. colors of points will be
# selected randomly, but you can also provide color.vector.
# non.interactive can be set to a filename
cluster.visualize <- function(db, cls, size.cutoff, distmat=NULL,
        color.vector=NULL, non.interactive="", cluster.result=1,
        dimensions=2, quiet=FALSE, highlight.compounds=NULL, 
        highlight.color=NULL, ...)
{
    ## ThG: added for compatibility with new S4 classes APset/AP
    dbtype <- as.character(class(db))
    if(dbtype=="APset") { db <- apset2descdb(db) }
    ## ThG: end of lines
    cluster_col <- cluster.result * 2 + 1
    cls_is_large <- cls[,2] >= size.cutoff
    cls_sel <- cls[cls_is_large,]
    cluster_ids <- levels(as.factor(cls_sel[,cluster_col]))

    if (is.null(color.vector)) {
        .colors <- rainbow(length(cluster_ids))
        color.vector <- sample(.colors, length(.colors))
    }

    if (is.null(distmat)) {
        descdb_sel <- db$descdb[cls_sel[,1]]
        distmat_sel <- .calc.distmat(descdb_sel, quiet=quiet, ...)
    } else
        distmat_sel <- distmat[cls_sel[,1], cls_sel[,1]]
    mds_points <- cmdscale(distmat_sel, k=dimensions)

    if (dimensions != 2) {
        coord <- cbind(mds_points, cls_sel[,cluster_col])
        rownames(coord) <- cls_sel$ids 
        colnames(coord) <-
            c(paste('V', 1:dimensions, sep=''), colnames(cls_sel)[cluster_col])
        return(coord)
    }
 
    if (non.interactive != "") {
        if (length(grep('\\.pdf$', non.interactive)) != 0)
            pdf(non.interactive)
    else if (length(grep('\\.ps$', non.interactive)) != 0 ||
            length(grep('\\.eps$', non.interactive)) != 0)
            postscript(non.interactive, onefile=TRUE, horizontal=FALSE)
    else {
            warning("The filename you supplied has an unsupported extension.",
                " We will add .pdf to the name.")
            pdf(paste(non.interactive, 'pdf', sep='.'))
    }
    }
    
    # set up the plot
    plot(range(mds_points[,1]), range(mds_points[,2]), type="n", xlab="",
        ylab="", main="Clustering Result")

    # check highlight
    if (!is.null(highlight.compounds)) {
        cls_is_highlight <- cls[,1] %in% highlight.compounds
        cls_sel_is_highlight <- cls_is_highlight[cls_is_large]
    }
    
    for (i in 1:length(cluster_ids)) {
        col <- color.vector[i%%length(color.vector) + 1]
        cid <- cluster_ids[[i]]
        if (! quiet) {
            cat(paste("cluster", cid, "colored", col))
            cat("\n")  
        }
        # indices of points in this cluster
        in_cluster <- seq(1, length(cls_sel[,1]))[cls_sel[,cluster_col] == cid]
        points(mds_points[in_cluster,1], mds_points[in_cluster,2], col=col)
        # check highlight
        if (!is.null(highlight.compounds)) {
            highlight_in_cluster <- cls_sel_is_highlight & cls_sel[,cluster_col] == cid
            if (is.null(highlight.color)) .col <- col
            else .col <- highlight.color
            points(mds_points[highlight_in_cluster,1], mds_points[highlight_in_cluster,2], col=.col, pch=19)
        }
    }

    # interactive
    all.clicked <- rep(FALSE, length(cls_sel[,1]))
    if (non.interactive == "") {
        cat("=============================================================\n")
        cat("| Click points in a plot to get information on compounds    |\n")
        cat("| they represent.                                           |\n")
        cat("|                                                           |\n")
        g_dev <- .Platform$GUI
        if (is.character(g_dev) && g_dev == 'X11') 
            cat(
            "| right click on the plot to stop.                          |\n")
        else
            cat(
            "| press ESC key in the plot window to stop.                 |\n")
        cat("=============================================================\n\n")
        while(TRUE) {
            clicked <- identify(mds_points[,1], y=mds_points[,2],
                            labels=cls_sel$ids, n=1)
            if (length(clicked) == 0) break
            index <- cls_sel$ids[clicked]
            all.clicked[clicked] <- TRUE
            output <- data.frame(c(db$cids[index], cls[cls$ids==index,]))
            names(output) <- c('Compound ID', colnames(cls))
            print(output)
        }
    } else {
        dev.off()
    }

    coord <- cbind(mds_points, cls_sel[,cluster_col], all.clicked)
    rownames(coord) <- cls_sel$ids 
    colnames(coord) <- c(paste('V', 1:dimensions, sep=''),
                        colnames(cls_sel)[cluster_col], 'Clicked')
    return(coord)
}

# fast detection of duplicated compounds
# first all descriptors are concatenated as character strigns. then `duplicated'
# are called to detect (potentially) duplicated compounds
cmp.duplicated <- function(db, sort=FALSE, type=1)
{
    ## ThG: added for compatability with new S4 classes APset/AP
    dbtype <- as.character(class(db))
    if(dbtype=="APset") { db <- apset2descdb(db) }
    ## ThG: end of lines
    if (!sort)
        f <- function(x) {paste(x, collapse=" ")}
    else
        f <- function(x) {paste(sort(x), collapse=" ")}

    db_str <- lapply(db$descdb, f)
    ## ThG: added to also return duplicates in a cluster data frame
    if(type==1) { 
	dup <- duplicated(db_str) 
    }
    if(type==2) {
    	names(db_str) <- db$cids; 
    	db_str <- tapply(names(db_str), as.character(db_str), paste)
    	names(db_str) <- 1:length(db_str)
    	clsz <- sapply(db_str, length)
    	dup <- data.frame(ids=unlist(db_str), CLSZ_100=rep(clsz, clsz), CLID_100=rep(1:length(db_str), clsz))
    }
    ## ThG: end of lines
    return(dup)
}

# helper function used by cluster.visualization
# given a set of points and a point, find the nearest match of point in the set
.search.nearest <- function(mds_points, point)
{
    best <- 1
    dist <- (mds_points[1,1] - point$x)^2  + (mds_points[1,2] - point$y)^2
    for (i in 2:length(mds_points[,1])) {
        d <- (mds_points[i,1] - point$x)^2 + (mds_points[i,2] - point$y)^2
        if (d < dist) {
           best <- i
           dist <- d
        }
    }
    return(best)
}

# intertanl cluster procedure. It takes a leader compounds, and use graph
# traversal to find all compounds that should be grouped in the cluster leaded
# by the leader. If multiple cutoffs are given, it will use the loosest one;
# but at the same time, clusters for stricter cutoffs are also built as a
# desired side-effect, without need to recompute the distance.
# RETURN: the new cluster_id matrix
.cluster_g <- function(lead, all, dist_func, cutoff, cluster_id,
symmetric=TRUE, quiet=FALSE)
{
    n_cutoff <- length(cutoff)
    # where to start in finding the elements in cluster
    start_from <- 1
    if (symmetric)
        start_from <- lead + 1
    # queue to store the visited and to-visit nodes
    # queue will be sorted by number of "unknowns"
    queue <- c(lead)
    if (cluster_id[lead, 1] == 0)
        cluster_id[lead,] <- lead
    current_id <- cluster_id[lead, ]
    # head points to node currently being visited
    head <- 1
    # tail points to tail of queue, ie where to add new nodes
    tail <- 2
    # do we really need to search?
    if (start_from > length(all))
        return(cluster_id)
    while (head != tail) {
        if (! quiet) .progress_bar()
        # pop from head
        elem <- queue[[head]]
        # init cluster IDs for elem: clear the unknowns
        cluster_id[elem,cluster_id[elem,] == 0] <- elem
        # find neighbors
        for (i in all[start_from:length(all)]) {
            # does i contain any unknown? only need to check
            # strictest cutoff
            if (cluster_id[i,n_cutoff] != 0)
                next
            # else, update
            d = dist_func(elem, i)
            indices <- (d <= cutoff)
            cluster_id[i, indices] <- cluster_id[elem, indices]
            # satisfying the loosest cutoff?
            if (indices[[1]]) {
                # add i to queue
                num_unknowns <- sum(cluster_id[i,] == 0)
                # is it already in the queue?
                if (length(queue[queue == i]))
                    pos <- (1:length(queue))[queue == i]
                else {
                    pos <- tail
                    tail <- tail + 1
                }
                while (pos > head && num_unknowns <
                        sum(cluster_id[queue[[pos - 1]],] == 0)) {
                    queue[[pos]] <- queue[[pos - 1]]
                    pos <- pos - 1
                }
                queue[[pos]] <- i
            }
        }
        head <- head + 1
    }
    return(cluster_id)
}

###############################
## Jarvis-Patrick Clustering ##
###############################
## Added by ThG on 28-Oct-12
## Function to perform Jarvis-Patrick clustering. The algorithm requires a
## nearest neighbor table, which consists of 'j' nearest neighbors for each item
## in the dataset. This information is then used to join items into clusters
## that share at least 'k' nearest neighbors. The values for 'j' and 'k' are
## user-defined parameters. The jarvisPatrick() function can generate the
## nearest neighbor table for APset and FPset objects and then perform Jarvis-
## Patrick clustering on that table. It also accepts a precomputed nearest 
## neighbor table in form of an object of class matrix. The output is a cluster
## vector with the item labels in the name slot and the cluster IDs in the data
## slot. Alternatively, the function can return the nearest neighbor matrix.
## As third parameter the user can set a minimum similarity value for generating 
## the nearest neighbor table. The latter is an optional setting that is not part 
## of the original Jarvis-Patrick algorithm. It allows to generate more tight 
## clusters and minimizes some limitations of this method, such as joining unrelated
## items when clustering small datasets.  
jarvisPatrick <- function(x, j, k, cutoff=NA, type="cluster", ...) {      
        ## Check inputs
        if(!any(c("APset", "FPset", "matrix") %in% class(x))) stop("class(x) needs to be APset, FPset or matrix")
        ## If class(x) is APset or FPset, generate nearest neighbor matrix (nnm)
        if(any(c("APset", "FPset") %in% class(x))) {
                if(is.na(cutoff)) { # Standard Jarvis-Patrick clustering without cutoff
                        if(class(x)=="FPset") {
                                nnm <- t(sapply(seq(along=x), function(y) names(fpSim(x[y], x, top=j, ...))))
                        } 
                        if(class(x)=="APset") {
                                nnm <- t(sapply(seq(along=x), function(y) names(cmp.search(x, x[y], type=2, cutoff=j, quiet = TRUE, ...))))
                        }
                } 
                if(is.numeric(cutoff) & cutoff <= 1) { # Non-standard Jarvis-Patrick clustering with cutoff
                        nnm <- matrix(NA, length(x), j)
                        if(class(x)=="FPset") {
                                for(i in seq(along=nnm[,1])) {
                                        tmp <- names(fpSim(x[i], x, cutoff=cutoff, top=j, ...))
                                        nnm[i,1:length(tmp)] <- tmp
                                }
                        }
                        if(class(x)=="APset") {
                                for(i in seq(along=nnm[,1])) {
                                        tmp <- names(cmp.search(x, x[i], type=2, cutoff=cutoff, quiet = TRUE, ...)[1:j])
                                        nnm[i,1:length(tmp)] <- tmp
                                }
                        }
                }
                rownames(nnm) <- cid(x)
                colnames(nnm) <- seq(along=nnm[1,])
        }
        if(type=="matrix") {
                return(nnm) 
        }
        if(type=="cluster") {
                ## Run Jarvis-Patrick clustering on nearest neighbor matrix (nnm)
                ## (i) Initialize algorithm by generating numeric vector with increasing cluster numbers
                clusters <- seq(along=nnm[,1]); names(clusters) <- rownames(nnm) 
                "%in%" <- function(x, table) match(x, table, nomatch = 0, incomparables=NA) > 0 # Sets NAs of %in% function as incompatible entries which is required when non-standard cutoff is used.
                for(i in seq(along=nnm[,1])) {
                        ## (ii) Identify for each nnm row the remaining rows with at least k common neighbors
                        nncount <- rowSums(matrix(as.character(t(nnm[i:length(nnm[,1]),])) %in% nnm[i,], length(nnm[,1])-i+1, length(nnm[1,]), byrow=TRUE)) 
                        above <- nncount >= k
                        ## (iii) Assign to new members always smallest cluster number to maintain assignments of previous 
                        ##       iterations. This allows single pass through nnm and results in a single linkage like behavior.
                        if(length(clusters[i:length(clusters)][above])!=0) {
                                clusters[i:length(clusters)][above] <- min(clusters[i:length(clusters)][above])
                        }
                }
                ## Assign continuous numbers as cluster names
                clusterstmp <- sort(clusters)
                tmp <- 1:length(unique(clusterstmp)); names(tmp) <- unique(clusterstmp)
                tmp <- tmp[as.character(clusterstmp)]; names(tmp) <- names(clusterstmp)
                clusters <- tmp[names(clusters)]
                return(clusters)
        }
}
## Usage:
# library(ChemmineR)
# data(apset)
# fpset <- desc2fp(apset)
# jarvisPatrick(x=apset, j=6, k=5)
# jarvisPatrick(x=fpset, j=6, k=2, cutoff=0.4)
# jarvisPatrick(x=fpset, j=2, k=2, type="matrix")







