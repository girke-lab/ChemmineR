
ErrorClass <- "try-error" 


postgresqlWriteTable <- function(con, name, value, field.types, row.names = TRUE,
                                 overwrite = FALSE, append = FALSE, ..., allow.keywords = FALSE) {
	 if(!require(RPostgreSQL))
		 stop("The RPostgreSQL package is required to use this function (postgresqlWriteTable)")
    if(overwrite && append)
        stop("overwrite and append cannot both be TRUE")
    if(!is.data.frame(value))
        value <- as.data.frame(value)
    if(row.names){
        value <- cbind(row.names(value), value)  ## can't use row.names= here
        names(value)[1] <- "row.names"
    }
    if(missing(field.types) || is.null(field.types)){
        ## the following mapping should be coming from some kind of table
        ## also, need to use converter functions (for dates, etc.)
        field.types <- sapply(value, dbDataType, dbObj = con)
    }

    i <- match("row.names", names(field.types), nomatch=0)
    if(i>0) ## did we add a row.names value?  If so, it's a text field.
        ## MODIFIED -- Sameer
        field.types[i] <- dbDataType(dbObj=con, field.types[row.names])
    new.con <- con

    if(dbExistsTable(con,name)){
        if(overwrite){
            if(!dbRemoveTable(con, name)){
                warning(paste("table", name, "couldn't be overwritten"))
                return(FALSE)
            }
        }
        else if(!append){
            warning(paste("table",name,"exists in database: aborting assignTable"))
            return(FALSE)
        }
    }
    if(!dbExistsTable(con,name)){      ## need to re-test table for existance
        ## need to create a new (empty) table
        sql1 <- paste("create table ", postgresqlTableRef(name), "\n(\n\t", sep="")
        sql2 <- paste(paste(postgresqlQuoteId(names(field.types)), field.types), collapse=",\n\t",
                      sep="")
        sql3 <- "\n)\n"
        sql <- paste(sql1, sql2, sql3, sep="")
        rs <- try(dbSendQuery(new.con, sql))
        if(inherits(rs, ErrorClass)){
            warning("could not create table: aborting assignTable")
            return(FALSE)
        } else {
            dbClearResult(rs)
        }
    }

    ## convert columns we can't handle in C code
    value[] <- lapply(value, function(z) {
        if(is.object(z) && !is.factor(z)) as.character(z) else z
    })
    oldenc <- dbGetQuery(new.con, "SHOW client_encoding")
    postgresqlpqExec(new.con, "SET CLIENT_ENCODING TO 'UTF8'")
    sql4 <- paste("COPY  ", postgresqlTableRef(name),"(",paste(names(value),collapse=","),") FROM STDIN")
    postgresqlpqExec(new.con, sql4)
    postgresqlCopyInDataframe(new.con, value)
    rs<-postgresqlgetResult(new.con)

    retv <- TRUE
    if (inherits(rs, ErrorClass)) {
        warning("could not load data into table")
        retv <- FALSE
    }

    dbClearResult(rs)
    sql5 <- paste("SET CLIENT_ENCODING TO '", oldenc, "'", sep="")
    dbGetQuery(new.con, sql5)

    retv
}


