
debug = FALSE

dbOp<-function(dbExpr){
	#print(as.character(substitute(dbExpr)))
	#print(system.time(dbExpr))
	dbExpr
}




initDb <- function(handle){

	if(is.character(handle)){ #assume sqlite filename
		require(RSQLite)
		driver = dbDriver("SQLite")
		conn = dbConnect(driver,dbname=handle,cache_size=100000)
	}else if(inherits(handle,"DBIConnection")){ #the user passed in a connection object
		conn=handle
	}else{
		stop("handle must be a SQLite database name or a DBIConnection")
	}

	tableList=dbListTables(conn)

	if( ! all(c("compounds","descriptor_types","descriptors") %in% tableList)) {
		print("createing db")

		statements = unlist(strsplit(paste(
							  readLines(system.file("schema/compounds.SQLite",package="ChemmineR",mustWork=TRUE)),
							  collapse=""),";",fixed=TRUE))
		#print(statements)

		Map(function(sql) dbOp(dbGetQuery(conn,sql)),statements)
	}
	conn
}

dbTransaction <- function(conn,expr){
	tryCatch({
		dbGetQuery(conn,"BEGIN TRANSACTION")
		ret=expr
		dbCommit(conn)
		ret
	},error=function(e){
		dbRollback(conn)
		stop(paste("db error inside transaction: ",e$message))
	})
}
loadDb <- function(conn,data,featureGenerator){

	if(debug) print(paste("loading ",paste(dim(data),collapse=" "),"compounds"))

	data=cbind(data,definition_checksum=sapply(as.vector(data[,"definition"]),
													  function(def) digest(def,serialize=FALSE) ))

	numCompounds = dim(data)[1]
	dups=0

	addNeededFeatures(conn,data,featureGenerator)

	if(all(c("name","definition","format") %in% colnames(data))){
		insertNamedDef(conn,data)
	}else if(all(c("definition","format") %in% colnames(data))){
		insertDef(conn,data)
	}else {
		stop("given data must have columns either (definition,format), or (name,definiton,format)")
	}

	insertUserFeatures(conn,data)
	as.matrix(data["definition_checksum"],rownames.force=FALSE)
}
loadDescriptors <- function(conn,data){
	#expects a data frame with "definition_checksum" and "descriptor"

	req_columns=c("definition_checksum","descriptor","descriptor_type")
	if(!all(req_columns %in% colnames(data)))
		stop(paste("missing some names, found",paste(colnames(data),collapse=","),"need",paste(req_columns,collapse=",")))

	#ensure the needed descriptor types are available for the insertDescriptor function to use
	unique_types = unique(data[["descriptor_type"]])
	all_descriptor_types=dbGetQuery(conn,"SELECT distinct descriptor_type FROM descriptor_types")[[1]]
	newTypes = setdiff(unique_types,all_descriptor_types)
	if(length(newTypes) > 0)
		insertDescriptorType(conn,data.frame(descriptor_type=newTypes))

	insertDescriptor(conn,data)
}
addNeededFeatures <- function(conn,data, featureGenerator){
	if(dim(data)[1] == 0){ # no data given
		warning("no data given to addNeededFeatures, doing nothing")
		return()
	}

	features = featureDiff(conn,data)

	if(length(features$missing) != 0)
		stop(paste("missing features:",paste(features$missing,collapse=",")))


	#will need to query all existing compounds and add these features
	if(length(features$new) != 0){
		message("adding new features to existing compounds. This could take a while")
		lapply(features$new,function(name) createFeature(conn,name,is.numeric(data[,name])))
		indexExistingCompounds(conn,features$new,featureGenerator)
	}
}
featureDiff <- function(conn,data) {

	if( nrow(data) == 0){ # no data given
		list(new=c(),missing=c())
	}else{

		compoundFields = dbListFields(conn,"compounds")
		userFieldNames = setdiff(colnames(data),compoundFields)

		tableList=dbListTables(conn)
		existingFeatures = sub("^feature_","",tableList[grep("^feature_.*",tableList)])

		missingFeatures = setdiff(existingFeatures,userFieldNames)
		newFeatures = setdiff(userFieldNames,existingFeatures)

		list(new=newFeatures, missing=missingFeatures)
	}
}
insertUserFeatures<- function(conn,data){
	if(debug) print("inserting user features")

	if(dim(data)[1] == 0){ # no data given
		warning("no data given to insertUserFeatures, doing nothing")
		return()
	}


	compoundFields = dbListFields(conn,"compounds")
	userFieldNames = setdiff(colnames(data),compoundFields)
	if(length(userFieldNames)==0) return()

	tableList=dbListTables(conn)
	existingFeatures = sub("^feature_","",tableList[grep("^feature_.*",tableList)])
	if(length(existingFeatures)==0) return()


		
	#index all new compounds for all features

	#fetch newly inserted compounds
	
	df = dbGetQuery(conn,paste("SELECT compound_id,definition_checksum FROM compounds LEFT JOIN feature_",
								 userFieldNames[1]," as f USING(compound_id)  WHERE f.compound_id IS
								 NULL",sep=""))
	data= merge(data,df)
	sapply(userFieldNames,function(name) insertFeature(conn,name,data))

}


########### Large query/file utilities ################
bufferLines <- function(fh,batchSize,lineProcessor){
	while(TRUE){
		lines = readLines(fh,n=batchSize)
		if(length(lines) > 0)
			lineProcessor(lines)
		else
			break;
	}
}
bufferResultSet <- function(rs,rsProcessor,batchSize=1000,closeRS=FALSE){
	while(TRUE){
		chunk = fetch(rs,n=batchSize)
		if(dim(chunk)[1]==0) # 0 rows in data frame
			break;
		rsProcessor(chunk)
	}
	if(closeRS) dbClearResult(rs)
}
batchByIndex <- function(allIndices,indexProcessor, batchSize=100000){

	numIndices=length(allIndices)

	if(numIndices==0)
		return()

	for(start in seq(1,numIndices,by=batchSize)){
		end = min(start+batchSize-1,numIndices)
		if(debug) print(paste(start,end))
		indexSet = allIndices[start:end]
		indexProcessor(indexSet)

	}
}
#this does not guarentee a consistant ordering of the result
selectInBatches <- function(conn, allIndices,genQuery,batchSize=100000){
	#print(paste("all indexes: ",paste(allIndices,collapse=", ")))

	indices = unique(allIndices)
	#TODO: possibly pre-allocate result here, if performance is a problem
	result=NA
	batchByIndex(indices, function(indexBatch){
			#print(paste("query:",genQuery(indexBatch)))
			df = dbGetQuery(conn,genQuery(indexBatch))
			#print(paste("got",paste(dim(df),collapse=" "),"results"))
			result <<- if(is.na(result)) df else  rbind(result,df)
			#print(paste("total results so far: ",length(result)))
	},batchSize)
	#print(paste("final count: ",length(result)))
	result

}
###############################################################


definition2SDFset <- function(defs){
	read.SDFset(unlist(strsplit(defs,"\n",fixed=TRUE)))
}

loadSdf <- function(conn,sdfFile,fct=function(x) data.frame(),
						  descriptors=function(x) data.frame(descriptor=c(),descriptor_type=c()), 
						  Nlines=10000, startline=1, restartNlines=100000){

	if(inherits(sdfFile,"SDFset")){
		if(debug) print("loading SDFset")
		sdfset=sdfFile

		sdfstrList=as(as(sdfset,"SDFstr"),"list")
		names = unlist(Map(function(x) x[1],sdfstrList))
		defs = unlist(Map(function(x) paste(x,collapse="\n"), sdfstrList) )

		systemFields=data.frame(name=names,definition=defs,format="sdf")

		userFeatures <- fct(sdfset)
		allFields = if(length(userFeatures)!=0) cbind(systemFields,userFeatures) else systemFields
		cmdIds=dbTransaction(conn,{
			checksums=loadDb(conn,allFields,fct)

			#We assume descriptors are in the same order as compounds
			descriptor_data = descriptors(sdfset)
			if(length(descriptor_data) != 0)
				loadDescriptors(conn,cbind(checksums,descriptor_data))

			findCompoundsByChecksum(conn,checksums)
	   })
		return(cmdIds)
	}
	compoundIds=c()
	## Define loop parameters 
	stop <- FALSE 
	f <- file(sdfFile, "r")
	n <- Nlines
	offset <- 0
	## For restarting sdfStream at specific line assigned to startline argument. If assigned
        ## startline value does not match the first line of a molecule in the SD file then it 
        ## will be reset to the start position of the next molecule in the SD file.
	if(startline!=1) { 
		fmap <- file(sdfFile, "r")
		shiftback <- 2
		chunkmap <- scan(fmap, skip=startline-shiftback, nlines=restartNlines, what="a", blank.lines.skip=FALSE, quiet=TRUE, sep ="\n")
		startline <- startline + (which(grepl("^\\${4,4}", chunkmap, perl=TRUE))[1] + 1 - shiftback)
		if(is.na(startline)) stop("Invalid value assigned to startline.")
		dummy <- scan(f, skip=startline-2, nlines=1, what="a", blank.lines.skip=FALSE, quiet=TRUE, sep ="\n")
		close(fmap)
		offset <- startline - 1 # Maintains abolut line positions in index
	}
	counter <- 0
	cmpid <- 1
	partial <- NULL

	dbTransaction(conn,{
		while(!stop) {
			counter <- counter + 1
			chunk <- scan(f, n=n, what="a", blank.lines.skip=FALSE, quiet=TRUE, sep ="\n") # scan has more flexibilities for reading specific line ranges in files.
			if(length(chunk) > 0) {
				if(length(partial) > 0) {
					chunk <- c(partial, chunk)
				}
				## Assure that lines of least 2 complete molecules are stored in chunk if available
				inner <- sum(grepl("^\\${4,4}", chunk, perl=TRUE)) < 2
				while(inner) {
					chunklength <- length(chunk)
					chunk <- c(chunk, readLines(f, n = n))
					if(chunklength == length(chunk)) { 
						inner <- FALSE 
					} else {
						#inner <- sum(grepl("^\\${4,4}", chunk, perl=TRUE)) < 2
						inner <- sum(grepl("$$$$", chunk, fixed=TRUE)) < 2
					}
				}
				y <- regexpr("^\\${4,4}", chunk, perl=TRUE) # identifies all fields that start with a '$$$$' sign
				index <- which(y!=-1)
				indexDF <- data.frame(start=c(1, index[-length(index)]+1), end=index)
				complete <- chunk[1:index[length(index)]]
				if((index[length(index)]+1) <= length(chunk)) {
					partial <- chunk[(index[length(index)]+1):length(chunk)]
				} else {
					partial <- NULL
				}

				sdfset <- read.SDFset(read.SDFstr(complete))
				valid = validSDF(sdfset)
				sdfset=sdfset[valid]


				userFeatures <- fct(sdfset)

		#		sdfstrList=as(as(sdfset,"SDFstr"),"list")
		#		names = unlist(Map(function(x) x[1],sdfstrList))
		#		defs = unlist(Map(function(x) paste(x,collapse="\n"), sdfstrList) )

				defs=apply(indexDF,1,function(row) paste(complete[row[1]:row[2]],collapse="\n"))
				names = complete[indexDF[,1]]
				systemFields=data.frame(name=names[valid],definition=defs[valid],format="sdf")

				allFields = if(length(userFeatures)!=0) cbind(systemFields,userFeatures) else systemFields
				checksums=loadDb(conn,allFields,fct)

				descriptor_data = descriptors(sdfset)
				if(length(descriptor_data) != 0)
					loadDescriptors(conn,cbind(checksums,descriptor_data))

				compoundIds = c(compoundIds,findCompoundsByChecksum(conn,checksums))
				
			}
			if(length(chunk) == 0) {
				stop <- TRUE
				close(f)
			}
		}

		compoundIds
	})
}



loadSmiles <- function(conn, smileFile,batchSize=10000){

	f = file(smileFile,"r")
	compoundIds = c()
	bufferLines(f,batchSize=batchSize,function(lines) 
					compoundIds <<- c(compoundIds,
											findCompoundsByChecksum(conn,loadDb(conn,data.frame(definition=lines,format="smile")))))
	close(f)
	compoundIds
}

findCompounds <- function(conn,featureNames,tests){

	# SELECT compound_id FROM compounds join f1 using(compound_id) join f2 using(compound_id)
	# ... where test1 AND test2 AND ...
	featureTables = paste("feature_",featureNames,sep="")
	tryCatch({
		sql = paste("SELECT compound_id FROM compounds JOIN ",paste(featureTables,collapse=" USING(compound_id) JOIN "),
					" USING(compound_id) WHERE ",paste("(",paste(tests,collapse=") AND ("),")") ) 
	
		#print(paste("query sql:",sql))
		result = dbGetQuery(conn,sql)
		result[1][[1]]
	},error=function(e){
		if(length(grep("no such column",e$message))!=0){
			stop("featureNames must contain every feature used in a test")
		}else{
			stop(paste("error in findCompounds:",e$message))
		}
	})

}

#undefined ordering by default
findCompoundsByChecksum <- function(conn,checksums,keepOrder=FALSE)
	findCompoundsByX(conn,"definition_checksum",checksums,keepOrder)
findCompoundsByName<- function(conn,names,keepOrder=FALSE)
	findCompoundsByX(conn,"name",names,keepOrder)

findCompoundsByX<- function(conn,fieldName,data,keepOrder=FALSE){
	result = selectInBatches(conn,data,function(batch) 
			paste("SELECT compound_id,",fieldName
					," FROM compounds WHERE ",fieldName," IN
					('",paste(batch,collapse="','"),"')",sep=""),1000)
	ids = result$compound_id
	if(length(ids)!=length(data))
		stop(paste("found only",length(ids),"out of",length(data),
					  "queries given"))
	if(keepOrder){
		names(ids)=result[[fieldName]]
		ids[data]
	}else{
		ids
	}
}
getCompounds <- function(conn,compoundIds,filename=NA){
	
	processedCount=0
	if(!is.na(filename)){
		f=file(filename,"w")

		resultProcessor = function(rows){
			lapply(rows[2][[1]],function(def) cat(def,"\n",sep="",file=f))
			processedCount <<- processedCount + dim(rows)[1]
		}
	}else{
		sdfset = rep(NA,length(compoundIds))
		count=1
		resultProcessor = function(rows){
					l=length(rows[2][[1]])
					sdfset[count:(count+l-1)]<<-as(definition2SDFset(rows[2][[1]]),"SDF")
					names(sdfset)[count:(count+l-1)]<<-as.character(rows[1][[1]])
					count<<-count+l
					processedCount <<- processedCount + dim(rows)[1]
			}
	}

	batchByIndex(compoundIds,function(compoundIdSet){
		rs = dbSendQuery(conn,
							  paste("SELECT compound_id,definition FROM compounds where compound_id in (",
									  paste(compoundIdSet,collapse=","),")"))
		bufferResultSet(rs, resultProcessor,1000)
		dbClearResult(rs)
	})

	if(length(compoundIds) != processedCount) {
		if(debug) print(str(sdfset))
		warning(paste("not all compounds found,",length(compoundIds),"given but",processedCount,"found"))
	}

	if(!is.na(filename)){
		close(f)
	}else{
		return(as(sdfset,"SDFset"))
	}
}
getCompoundNames <- function(conn, compoundIds){

	result = selectInBatches(conn,compoundIds,function(ids)
					  paste("SELECT compound_id, name FROM compounds where compound_id in (",
									  paste(ids,collapse=","),")"))

	n=result$name
	names(n)=result$compound_id
	#print(n[as.character(compoundIds)])
	n[as.character(compoundIds)]
	#as.matrix(merge(data.frame(compound_id=compoundIds),result,sort=FALSE)[[2]])
}
indexExistingCompounds <- function(conn,newFeatures,featureGenerator){

	# we have to batch by index because we need to execute an insert statment along
	# the way and R DBI does not allow you to do two things at once.
	batchByIndex(dbGetQuery(conn,"SELECT compound_id FROM
									compounds WHERE format!='junk' ")[1][[1]],function(compoundIdSet){
		tryCatch({
				rows=dbGetQuery(conn,paste("SELECT compound_id,definition FROM compounds WHERE compound_id in (",
								  paste(compoundIdSet,collapse=","),")"))

				sdfset = definition2SDFset(rows[2][[1]])
				userFeatures  = featureGenerator(sdfset)
				lapply(newFeatures, function(name){
						 insertFeature(conn,name,cbind(compound_id=rows[1][[1]],userFeatures[name]))
				})
			},error=function(e) stop(paste("error in indexExistingCompounds:",e$message))
		)

	},1000)
}
addNewFeatures <- function(conn,featureGenerator){

	firstBatch = TRUE
	newFeatures = c()
	batchByIndex(dbGetQuery(conn,"SELECT compound_id FROM
									compounds WHERE format!='junk' ")[1][[1]],function(compoundIdSet){
		tryCatch({
				rows=dbGetQuery(conn,paste("SELECT compound_id,definition FROM compounds WHERE compound_id in (",
								  paste(compoundIdSet,collapse=","),")"))

				sdfset = definition2SDFset(rows[2][[1]])

				data = featureGenerator(sdfset)

				if(firstBatch){
					firstBatch<<-FALSE
					features = featureDiff(conn,data)
					newFeatures <<- features$new
					for(name in features$new)
						createFeature(conn,name,is.numeric(data[[name]]))
				}
				if(debug) print(paste("new features ",paste(newFeatures,collapse=",")))
				lapply(newFeatures, function(name){
						 insertFeature(conn,name,cbind(compound_id=rows[1][[1]],data[name]))
				})
			},error=function(e) stop(paste("error in indexExistingCompounds:",e$message))
		)

	},1000)


}

createFeature <- function(conn,name, isNumeric){


	sqlType = if(isNumeric) "NUMERIC" else "TEXT"
	if(debug) print(paste("adding",name,", sql type: ",sqlType))
	dbGetQuery(conn,
		paste("CREATE TABLE feature_",name," (
			compound_id INTEGER PRIMARY KEY REFERENCES compound(compound_id) ON DELETE CASCADE, ",
			name," ",sqlType," )",sep=""))
	#print("made table")
	dbGetQuery(conn,paste("CREATE INDEX feature_",name,"_index ON
								 feature_",name,"(",name,")",sep=""))
	#print("made index")

}

rmDups <- function(data,columns) data[!duplicated(data[,columns]),]
insertDef <- function(conn,data)  {
	data = rmDups(data,"definition_checksum")
	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn,paste("INSERT INTO compounds(definition,definition_checksum,format) ",
								 "VALUES(:definition,:definition_checksum,:format)",sep=""), bind.data=data)
	}else{
		apply(data,1,function(row) dbOp(dbGetQuery(conn, 
						 paste("INSERT INTO compounds(definition,definition_checksum,format)
																			  VALUES('",row[1],"','",row[2],"','",row[3],"')", sep=""))))
	}
}

insertNamedDef <- function(conn,data) {
	data = rmDups(data,"definition_checksum")
	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn,paste("INSERT INTO compounds(name,definition,definition_checksum,format) ",
								 "VALUES(:name,:definition,:definition_checksum,:format)",sep=""), bind.data=data)
	}else{
		apply(data,1,function(row) dbOp(dbGetQuery(conn, 
						 paste("INSERT INTO compounds(name,definition,definition_checksum,format) VALUES('",
																	  row[1],"','",row[2],"','",row[3],"','",row[4],"')", sep=""))))
	}
}

insertFeature <- function(conn,name,values){
	if(debug) print(paste("name:",name))
	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn, paste("INSERT INTO feature_",name,"(compound_id,",name,") ",
												 "VALUES(:compound_id,:",name,")",sep=""), bind.data=values)
	}else{
		apply(data,1,function(row) 
				dbGetQuery(conn,paste("INSERT INTO feature_",name,"(compound_id,",name,")
											 VALUES(",row[1],",", if(!is.numeric(row[2]))
													  {paste("'",row[2],"'",sep="")} else {row[2]},")")))
	}
}
insertDescriptor <- function(conn,data){
	data = rmDups(data,c("definition_checksum","descriptor_type"))
	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn, paste("INSERT INTO descriptors(compound_id, descriptor_type_id,descriptor) ",
				"VALUES( (SELECT compound_id FROM compounds WHERE definition_checksum = :definition_checksum),
							(SELECT descriptor_type_id FROM descriptor_types WHERE descriptor_type = :descriptor_type), 
							:descriptor )") ,bind.data=data)
	}else{
		apply(data,1,function(row) 
			dbGetQuery(conn,paste("INSERT INTO descriptors(compound_id, descriptor_type_id,descriptor) ",
				"VALUES( (SELECT compound_id FROM compounds WHERE definition_checksum = '",row["definition_checksum"] ,"'),
					(SELECT descriptor_type_id FROM descriptor_types WHERE descriptor_type = '",row["descriptor_type"],"'), 
						'",row["descriptor"],"' )" ) ))
	}
}
insertDescriptorType <- function(conn,data){
	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn,"INSERT INTO descriptor_types(descriptor_type) VALUES(:descriptor_type)",
								 bind.data=data)
	}else{
		apply(data,1,function(row) 
				dbGetQuery(conn,paste("INSERT INTO descriptor_types(descriptor_type) VALUES('",row[1],"')")))
	}
}
