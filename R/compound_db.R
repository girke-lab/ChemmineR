
debug = FALSE
#debug = TRUE

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
	
	enableForeignKeys(conn)

	tableList=dbListTables(conn)

	if( ! all(c("compounds","descriptor_types","descriptors") %in% tableList)) {
		print("createing db")

		sqlFile = file.path("schema",if(inherits(conn,"SQLiteConnection")) "compounds.SQLite" 
								  else if(inherits(conn,"PostgreSQLConnection")) "compounds.RPostgreSQL")
																	
		statements = unlist(strsplit(paste(
							  readLines(system.file(sqlFile,package="ChemmineR",mustWork=TRUE)),
							  collapse=""),";",fixed=TRUE))
		#print(statements)

		Map(function(sql) dbOp(dbGetQuery(conn,sql)),statements)
	}
	conn
}
enableForeignKeys <- function(conn){
	#SQLite needs this to enable foreign key constrains, which are off by default
	if(inherits(conn,"SQLiteConnection"))
		dbSendQuery(conn,"PRAGMA foreign_keys = ON")
}

dbTransaction <- function(conn,expr){
	tryCatch({

		# be paranoid about setting this as bad things will happen if its not set
		enableForeignKeys(conn)

		dbGetQuery(conn,"BEGIN TRANSACTION")
		ret=expr
		dbCommit(conn)
		ret
	},error=function(e){
		dbRollback(conn)
#		print(sys.calls())
		stop(paste("db error inside transaction: ",e$message))
	})
}
dbGetQueryChecked <- function(conn,statement,...){
	ret=dbGetQuery(conn,statement)
	err=dbGetException(conn)
	if(err$errorMsg[1] != "OK")
		stop("error in dbGetQuery: ",err$errorMsg,"  ",traceback())
	ret
}
loadDb <- function(conn,data,featureGenerator){

	names(data)=tolower(names(data))
	if(debug) print(paste("loading ",paste(dim(data),collapse=" "),"compounds"))

	if( ! ("definition_checksum" %in% colnames(data) ))
		data=cbind(data,definition_checksum=definitionChecksums(data[,"definition"]) )

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
definitionChecksums <- function(defs) {
	sapply(as.vector(defs),function(def) digest(def,serialize=FALSE),USE.NAMES=FALSE)
}
loadDescriptors <- function(conn,data){
	#expects a data frame with  req_columns

	req_columns=c("definition_checksum","descriptor","descriptor_type")
	if(!all(req_columns %in% colnames(data)))
		stop(paste("descriptor function is missing some fields, found",paste(colnames(data),collapse=","),
					  "need",paste(req_columns,collapse=",")))

	#ensure the needed descriptor types are available for the insertDescriptor function to use
	unique_types = unique(data[["descriptor_type"]])
	all_descriptor_types=dbGetQuery(conn,"SELECT distinct descriptor_type FROM descriptor_types")

	newTypes = if(nrow(all_descriptor_types)==0) unique_types 
				  else setdiff(unique_types,all_descriptor_types$descriptor_type)

	if(debug)
		print(paste("existing desc types:",paste(all_descriptor_types,collapse=","),"needed right now: ",
				paste(unique_types,collapse=",")," to be added: ",paste(newTypes,collapse=",")))
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

		existingFeatures = listFeatures(conn)

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
								 tolower(userFieldNames[1])," as f USING(compound_id)  WHERE f.compound_id IS
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
parBatchByIndex <- function(allIndices,indexProcessor,reduce,cl,batchSize=100000){

	numIndices=length(allIndices)

	if(numIndices==0)
		return()

	starts = seq(1,numIndices,by=batchSize)
	f = function(jobId){
		tryCatch({
				start = starts[jobId]
				end = min(start+batchSize-1,numIndices)
				ids=allIndices[start:end]
				indexProcessor(ids,jobId)
			},
			error=function(e){
				#write the error to a file to ensure it doesn't get lost
				cat(as.character(e),"\n",file=paste("error-",jobId,".out",sep=""))
				stop(e)
			}
		)
	}

	require(snow)
	#copy this to all nodes once so it is not copied for each iteration
	#of clusterApply
	clusterExport(cl,"allIndices",envir=environment())

	#we explicitly create an environment for f with just what it needs
	# so that serialization does not pull in un-nessacary things
	fEnv = new.env(parent=globalenv())
	fEnv$numIndices = numIndices
	fEnv$starts = starts
	fEnv$indexProcessor = indexProcessor
	fEnv$batchSize = batchSize

	environment(f) <- fEnv

	reduce(clusterApplyLB(cl,1:length(starts),f ))
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

loadSmiles <- function(conn, smileFile,...){
	loadSdf(conn,smile2sdfFile(smileFile),...)
}
loadSdf <- function(conn,sdfFile,fct=function(x) data.frame(),
						  descriptors=function(x) data.frame(descriptor=c(),descriptor_type=c()), 
						  Nlines=50000, startline=1, restartNlines=100000,updateByName=FALSE){

	if(inherits(sdfFile,"SDFset")){
		if(debug) print("loading SDFset")
		sdfset=sdfFile


		sdfstrList=as(as(sdfset,"SDFstr"),"list")
		names = unlist(Map(function(x) x[1],sdfstrList))
		defs = unlist(Map(function(x) paste(x,collapse="\n"), sdfstrList) )
		processAndLoad(conn,names,defs,sdfset,fct,descriptors,updateByName)
	}else{
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

					defs=apply(indexDF,1,function(row) paste(complete[row[1]:row[2]],collapse="\n"))[valid]
					names = complete[indexDF[,1]][valid]
					cmdIds = processAndLoad(conn,names,defs,sdfset,fct,descriptors,updateByName,inTransaction=TRUE)

					compoundIds = c(compoundIds,cmdIds)
					
				}
				if(length(chunk) == 0) {
					stop <- TRUE
					close(f)
				}
			}

			compoundIds
		})
	}
}

processAndLoad <- function(conn,names,defs,sdfset,featureFn,descriptors,updateByName,inTransaction=FALSE ) {
	# - compute checksums on defs
	# - if(updateByName)
	#		query checksums for each name
	#		exclude those that match our checksum
	#		updates <- those with different checksum
	#		new compounds <- everything else
	#	else
	#		exclude existing checksums
	#		insert the rest	
	# - include checksums in systemFields
	#		- have loadDb check for the existance of checksums and don't re-compute
	# - compute and insert/update descriptors for new/updated compounds

	checksums = definitionChecksums(defs)
	index=1:length(names)
	deleteCompIds=c()


	if(updateByName){
		#we assume compounds are unique by name and so changes in checksum
		# indicate updates to the same compounds

		if(debug) message("given names: ",paste(names,collapse=","))
		names(index)=names
		index=rev(index) # so last matches win


		existingByName = findCompoundsByX(conn,"name",names,allowMissing=TRUE,
														 extraFields = c("definition_checksum","name"))
		if(debug) {message("existing names: "); print(existingByName); }

#		idsToLoad = Filter(function(i) {
#			#either the name is completely new, or it exists, but the checksum is different, so this is 
#			# an update
#			(! names[i] %in% existingByName$name) || (! checksums[i] %in% existingByName$definition_checksum)
#
#		},index)


		namesToLoad = unique(Filter(function(name) {
			#either the name is completely new, or it exists, but the checksum is different, so this is 
			# an update
			(! name %in% existingByName$name) || (! checksums[index[name]] %in% existingByName$definition_checksum)
		},names))
		if(debug) message("names to load: ",paste(namesToLoad,collapse=","))

		namesToDelete = unique(Filter(function(name) {
			#delete those needing an update, so name exists and checksum is different
			( name %in% existingByName$name) && (! checksums[index[name]] %in% existingByName$definition_checksum)
		},names))
		if(debug) message("names to delete: ",paste(namesToDelete,collapse=","))

		#delete modified and missing compounds, will cascade to all descriptors

		deleteCompIds = c()
		for(name in namesToDelete){
			i = Position(function(x) x==name,existingByName$name,right=TRUE)
			deleteCompIds[checksums[index[name]]] = existingByName$compound_id[i]
		}


		#existingNameToCompId = existingByName$compound_id
		#names(existingNameToCompId) = existingByName$name
		#existingNameToCompId = rev(existingNameToCompId)
		#deleteCompIds = existingNameToCompId[namesToDelete]
		#rownames(existingByName)=existingByName$compound_id

		#rownames(existingByName)=existingByName$name
		#deleteCompIds = existingByName[namesToDelete,]$compound_id

		#if(length(namesToDelete) != 0)
			#names(deleteCompIds) = checksums[index[namesToDelete]]
		ids = index[namesToLoad]
		message("loading ",length(ids)," new compounds, updating ",length(deleteCompIds)," compounds")
	}else{
		#we do not assume names are unique, therefore if a checksum does not
		#exist, then it is added as if it where a new compound
		names(index)=checksums
		existingByChecksum = findCompoundsByX(conn,"definition_checksum",checksums,allowMissing=TRUE,
														 extraFields = c("definition_checksum","name"))

		if(debug){ message("existing checksums: "); print(existingByChecksum); }
		#select only those whose checksum does not already exist
		#setdiff also makes its result unique
		checksumsToLoad = setdiff(checksums,existingByChecksum$definition_checksum)
		if(debug) message("checksumsToLoad: ",paste(checksumsToLoad,collapse=","))

		ids = index[checksumsToLoad]
		message("loading ",length(ids)," new compounds")
	}
	if(debug) message("updating ids: ",paste(ids,collapse=","))

	tx = if(inTransaction) function(a,x) x  else dbTransaction
	if(length(ids)==0){
		tx(conn,deleteCompounds(conn,deleteCompIds))
		c() # no compound ids to return
	}else{

		names = names[ids]
		defs = defs[ids]
		checksums = checksums[ids]
		sdfset = sdfset[ids]
		userFeatures = featureFn(sdfset)

		if(debug) message("Features: ",colnames(userFeatures))
		if(debug) message("names: ",length(names)," defs: ",length(defs),", cksm: ",length(checksums),", features: ",length(userFeatures))

		systemFields=data.frame(name=names,definition=defs,format=rep("sdf",length(names)),definition_checksum=checksums)

		allFields = if(length(userFeatures)!=0) cbind(systemFields,userFeatures) else systemFields
		tx(conn,{
			deleteCompounds(conn,deleteCompIds)
			loadedChecksums=loadDb(conn,allFields,featureFn)

			#We assume descriptors are in the same order as compounds
			descriptor_data = descriptors(sdfset)
			if(length(descriptor_data) != 0)
				loadDescriptors(conn,cbind(loadedChecksums,descriptor_data))


			#print("original deleted comp ids")
			#print(deleteCompIds)
			#cmdIds=findCompoundsByChecksum(conn,loadedChecksums)
			#print("new comp ids: ")
			#print(cmdIds)

			# to keep stable compound id numbers, we reset the compound id of those
			# compounds that were merely modified.
			for(i in seq(along=deleteCompIds)){
				cksum=names(deleteCompIds)[i]
				id = deleteCompIds[i]
				dbSendQuery(conn,paste("UPDATE compounds SET compound_id = ",id,
											  " WHERE definition_checksum = '",cksum,"'",sep=""))
			}


			cmdIds=findCompoundsByChecksum(conn,loadedChecksums)
			if(nrow(loadedChecksums) != length(cmdIds)){
				stop("failed to insert all compounds. recieved ",nrow(loadedChecksums), 
					  " but only inserted ",length(cmdIds))
			}
			#print("updated comp ids:")
			#print(cmdIds)
		
			#those in deleteCompIds were not really added, just modified. so remove them here.
			#return newly added compound ids.
			setdiff(cmdIds,deleteCompIds)
		})
	}
}


smile2sdfFile <- function(smileFile,sdfFile=tempfile()){
	.ensureOB("smile format only suppported with ChemmineOB package")
	convertFormatFile("SMI","SDF",smileFile,sdfFile)
	sdfFile
}


deleteCompounds <- function(conn,compoundIds) {
	if(length(compoundIds)!=0)
		dbSendQuery(conn,paste("DELETE FROM compounds WHERE compound_id IN (",
							paste(compoundIds,collapse=","),")"))
}

#TODO: document feature names must be lower case in the tests clause
findCompounds <- function(conn,featureNames,tests){

	# SELECT compound_id FROM compounds join f1 using(compound_id) join f2 using(compound_id)
	# ... where test1 AND test2 AND ...
	featureNames = tolower(featureNames)

	featureTables = paste("feature_",featureNames,sep="")
	tryCatch({
		sql = paste("SELECT compound_id FROM compounds JOIN ",paste(featureTables,collapse=" USING(compound_id) JOIN "),
					" USING(compound_id) WHERE ",paste("(",paste(tests,collapse=") AND ("),")") ) 
	
		if(debug) print(paste("query sql:",sql))
		result = dbGetQueryChecked(conn,sql)
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
findCompoundsByChecksum <- function(conn,checksums,keepOrder=FALSE,allowMissing=FALSE)
	findCompoundsByX(conn,"definition_checksum",checksums,keepOrder,allowMissing)
findCompoundsByName<- function(conn,names,keepOrder=FALSE,allowMissing=FALSE)
	findCompoundsByX(conn,"name",names,keepOrder,allowMissing)

findCompoundsByX<- function(conn,fieldName,data,keepOrder=FALSE,allowMissing=FALSE,extraFields=c()){


	xf = if(length(extraFields)!=0) paste(",",paste(extraFields,collapse=",")) else ""
	result = selectInBatches(conn,data,function(batch) 
			paste("SELECT compound_id,",fieldName," ",xf,
					" FROM compounds WHERE ",fieldName," IN
					('",paste(batch,collapse="','"),"')",sep=""),1000)

	#no column names preserved for empty dataframes, so we can't just
	# handle it the same way, we need a special case :(
	if(length(result)==0)
		return(result)

	ids = result$compound_id

	if(!allowMissing && length(ids)!=length(data))
		stop(paste("found only",length(ids),"out of",length(data),
					  "queries given"))
	if(keepOrder){
		if(length(extraFields)!=0){
			rownames(result)=result[[fieldName]]
			result[as.character(data),c("compound_id",extraFields)]
		}
		else{
			names(ids)=result[[fieldName]]
			ids[as.character(data)]
		}
	}else{
		if(length(extraFields)!=0)
			result[,c("compound_id",extraFields)]
		else
			ids
	}
}
getCompounds <- function(conn,compoundIds,filename=NA,keepOrder=FALSE,allowMissing=FALSE){
	
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
		f=if(keepOrder)
			function(rows){
				rownames(rows) = rows$compound_id
				orderedIds = intersect(compoundIdSet,rows$compound_id)
				resultProcessor(rows[as.character(orderedIds),])
			}
		else
			f=resultProcessor

		bufferResultSet(rs,f,1000)
		dbClearResult(rs)
	})

	if(!allowMissing && length(compoundIds) != processedCount) {
		if(debug) print(str(sdfset))
		stop(paste("not all compounds found,",length(compoundIds),"given but",processedCount,"found"))
	}

	if(!is.na(filename)){
		close(f)
	}else{
		return(as(sdfset,"SDFset"))
	}
}
getCompoundNames <- function(conn, compoundIds,keepOrder=FALSE,allowMissing=FALSE){

	result = selectInBatches(conn,compoundIds,function(ids)
					  paste("SELECT compound_id, name FROM compounds where compound_id in (",
									  paste(ids,collapse=","),")"))

	if(!allowMissing)
		if(nrow(result) != length(compoundIds))
			stop(paste("found only",nrow(result),"out of",length(compoundIds), "ids given"))

	n=result$name
	if(keepOrder){
		names(n)=result$compound_id
		n[as.character(compoundIds)]
	}else
		n
}
indexExistingCompounds <- function(conn,newFeatures,featureGenerator){

	# we have to batch by index because we need to execute an insert statment along
	# the way and R DBI does not allow you to do two things at once.
	ids = dbGetQuery(conn,"SELECT compound_id FROM compounds WHERE format!='junk' ")
	if(length(ids) > 0)
		batchByIndex(ids[1][[1]],function(compoundIdSet){
			tryCatch({
					rows=dbGetQuery(conn,paste("SELECT compound_id,definition FROM compounds WHERE compound_id in (",
									  paste(compoundIdSet,collapse=","),")"))

					sdfset = definition2SDFset(rows[2][[1]])
					userFeatures  = featureGenerator(sdfset)
					names(userFeatures)=tolower(names(userFeatures))
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
	ids = dbGetQuery(conn,"SELECT compound_id FROM compounds WHERE format!='junk' ")
	if(length(ids) > 0)
		batchByIndex(ids[1][[1]],function(compoundIdSet){
			tryCatch({
					rows=dbGetQuery(conn,paste("SELECT compound_id,definition FROM compounds WHERE compound_id in (",
									  paste(compoundIdSet,collapse=","),")"))

					sdfset = definition2SDFset(rows[2][[1]])

					data = featureGenerator(sdfset)
					names(data)=tolower(names(data))

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
				},error=function(e) stop(paste("error in addNewFeature:",e$message))
			)

		},1000)


}

listFeatures <- function(conn){
	
	 tableList=dbListTables(conn)
	 sub("^feature_","",tableList[grep("^feature_.*",tableList)])
}

createFeature <- function(conn,name, isNumeric){

	sqlType = if(isNumeric) "NUMERIC" else "TEXT"
	if(debug) print(paste("adding",name,", sql type: ",sqlType))

	dbGetQueryChecked(conn,
		paste("CREATE TABLE feature_",name," (
			compound_id INTEGER PRIMARY KEY REFERENCES compounds(compound_id) ON DELETE CASCADE ON UPDATE CASCADE, ",
			"",name," ",sqlType," )",sep=""))

	#print("made table")
	dbGetQuery(conn,paste("CREATE INDEX feature_",name,"_index ON
								 feature_",name,"(\"",name,"\")",sep=""))
	#print("made index")

}

rmDups <- function(data,columns) data[!duplicated(data[,columns]),]
insertDef <- function(conn,data)  {
	data = rmDups(data,"definition_checksum")
	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn,paste("INSERT INTO compounds(definition,definition_checksum,format) ",
								 "VALUES(:definition,:definition_checksum,:format)",sep=""), bind.data=data)
	}else if(inherits(conn,"PostgreSQLConnection")){
		if(debug) print(data[,"definition_checksum"])
		fields = c("definition","definition_checksum","format")
		apply(data[,fields],1,function(row) dbOp(dbGetQuery(conn, 
						 "INSERT INTO compounds(definition,definition_checksum,format) VALUES($1,$2,$3)",
						 row)))
	}else{
		stop("database ",class(conn)," unsupported")
	}
}

insertNamedDef <- function(conn,data) {
	data = rmDups(data,"definition_checksum")
	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn,paste("INSERT INTO compounds(name,definition,definition_checksum,format) ",
								 "VALUES(:name,:definition,:definition_checksum,:format)",sep=""), bind.data=data)
	}else if(inherits(conn,"PostgreSQLConnection")){
		fields = c("name","definition","definition_checksum","format")
		postgresqlWriteTable(conn,"compounds",data[,fields],append=TRUE,row.names=FALSE)


		#apply(data[,fields],1,function(row) dbOp(dbGetQuery(conn, 
						 #"INSERT INTO compounds(name,definition,definition_checksum,format) VALUES($1,$2,$3,$4)",
						 #row)))
	}else{
		stop("database ",class(conn)," unsupported")
	}
}

insertFeature <- function(conn,name,values){
	if(debug) print(paste("name:",name))
	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn, paste("INSERT INTO feature_",name,"(compound_id,",name,") ",
												 "VALUES(:compound_id,:",name,")",sep=""), bind.data=values)
	}else if(inherits(conn,"PostgreSQLConnection")){
		fields = c("compound_id",name)

		postgresqlWriteTable(conn,paste("feature_",name,sep=""),values[,fields],append=TRUE,row.names=FALSE)


	#	apply(values[,fields],1,function(row) 
	#			dbGetQuery(conn,paste("INSERT INTO feature_",name,"(compound_id,\"",name,"\") ",
	#										 "VALUES($1,$2)",sep=""),row))
	}else{
		stop("database ",class(conn)," unsupported")
	}
}
insertDescriptor <- function(conn,data){
	data = rmDups(data,c("definition_checksum","descriptor_type"))
	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn, paste("INSERT INTO descriptors(compound_id, descriptor_type_id,descriptor) ",
				"VALUES( (SELECT compound_id FROM compounds WHERE definition_checksum = :definition_checksum),
							(SELECT descriptor_type_id FROM descriptor_types WHERE descriptor_type = :descriptor_type), 
							:descriptor )") ,bind.data=data)
	}else if(inherits(conn,"PostgreSQLConnection")){
		fields = c("definition_checksum","descriptor_type","descriptor")
		apply(data[,fields],1,function(row) 
			dbGetQuery(conn,paste("INSERT INTO descriptors(compound_id, descriptor_type_id,descriptor) ",
					"VALUES( (SELECT compound_id FROM compounds WHERE definition_checksum =$1 ),
								(SELECT descriptor_type_id FROM descriptor_types WHERE descriptor_type = $2 ), $3)"),
							row))
	}else{
		stop("database ",class(conn)," unsupported")
	}
}
insertDescriptorType <- function(conn,data){
	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn,"INSERT INTO descriptor_types(descriptor_type) VALUES(:descriptor_type)",
								 bind.data=data)
	}else if(inherits(conn,"PostgreSQLConnection")){
		apply(data,1,function(row) 
				dbGetQuery(conn,paste("INSERT INTO descriptor_types(descriptor_type) VALUES($1)"),row))
	}else{
		stop("database ",class(conn)," unsupported")
	}
}
