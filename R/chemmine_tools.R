# Purpose: R interface to ChemMine Tools
# by Tyler William H Backman

.serverURL <- "http://chemmine.ucr.edu/ChemmineR/"

# check job status
status <- function(object){
    if(class(object) != "jobToken"){
        stop("input not of class jobToken")
    }
    response <- postForm(paste(.serverURL, "jobStatus", sep=""), task_id=slot(object, "jobId"))[[1]]
    if(grepl("^ERROR:", response)){
        stop(response)
    }
    return(response)
}

# browse job online (works only once, saves to user account)
browseJob <- function(object){
    if(class(object) != "jobToken"){
        stop("input not of class jobToken")
    }
    url <- paste(.serverURL, "showJob", "/", slot(object, "jobId"), sep="")
    browseURL(url)
    return(url)
}

# get result 
result <- function(object){
    if(class(object) != "jobToken"){
        stop("input not of class jobToken")
    }
    response <- "RUNNING"
    while(response == "RUNNING"){
        Sys.sleep(2)
        response <- postForm(paste(.serverURL, "jobResult", sep=""), task_id=slot(object, "jobId"))[[1]]
    }
    if(grepl("^ERROR:", response)){
        stop(response)
    }
    if(response == "FAILED"){
        stop("Job Failed")
    }
    response <- .convertOutput(response, slot(object, "tool_name"))
    return(response)
}

# Purpose: retrieve list of all tools from server
listCMTools <- function(){
    response <- postForm(paste(.serverURL, "listCMTools", sep=""), category="all")[[1]]
    if(grepl("^ERROR:", response)){
        stop(response)
    }
    read.table(text=response, sep="\t", header=T)
}

# Purpose: retrieve details on a tool 
toolDetails <- function(tool_name){
    response <- postForm(paste(.serverURL, "toolDetails", sep=""), tool_name=tool_name)[[1]]
    if(grepl("^ERROR:", response)){
        stop(response)
    }
    cat(response)
}

# Purpose: launch a ChemMine Tools job on server
launchCMTool <- function(tool_name, input = "", ...){
    toolList <- listCMTools()
    if(! tool_name %in% toolList$Name){
        stop("invalid tool name")
    }
    input <- .convertInput(input, tool_name)
    response <- postForm(paste(.serverURL, "launchCMTool", sep=""), tool_name = tool_name, input = input, ...)[[1]]
    if(grepl("^ERROR:", response)){
        stop(response)
    }
    new("jobToken",
        tool_name = tool_name,
        jobId = response
    )
}

# Purpose: convert input to correct format
.convertInput <- function(input, toolName){
    response <- postForm(paste(.serverURL, "getConverter", sep=""), converterType="input", toolName=toolName)[[1]]
    objectClass <- gsub("\n.*", "", response)
    converter <- gsub("^.*?\n", "", response)

    if(objectClass == "data.frame"){
        try(input <- as.data.frame(input), silent = TRUE)
    } else {
        try(input <- as(input, objectClass), silent = TRUE)
    }

    if(class(input) != objectClass){
        stop(paste("input not of class", objectClass))
    }

    eval(parse(text = converter))
}

# Purpose: convert output to correct format 
.convertOutput <- function(output, toolName){
    response <- postForm(paste(.serverURL, "getConverter", sep=""), converterType="output", toolName=toolName)[[1]]
    converter <- gsub("^.*?\n", "", response)

    eval(parse(text = converter))
}

##################################
# Wrappers for old web functions #
##################################

# view sdfs in ChemMine Tools
sdf.visualize <- function(sdf){
    if(! class(sdf) == "SDFset"){
        stop('input not of class \"SDFset\"')
    } 
    job <- launchCMTool("sdf.visualize", sdf)
    browseJob(job)
}

# get CIDs from PubChem through ChemMine Web Tools
getIds <- function(cids) {
    if(! class(cids) == "numeric"){
        stop('reference compound ids must be of class \"numeric\"')
    }
    if(length(cids) == 0){
        stop('no compounds to retrieve- input must contain at least one cid')
    }
    jobToken <- launchCMTool("pubchemID2SDF", cids)
    result(jobToken)
}

# search PubChem through ChemMine Web Tools with smiles query
searchString <- function(smiles) {
    if(class(smiles) == "SMIset")
        smiles = as.character(smiles)
    if(! class(smiles) == "character"){
        stop('reference compound must be a smiles string of class \"character\"')
    } 	
    sdfquery <- smiles2sdf(smiles)
    searchSim(sdfquery)
}

# search PubChem through ChemMine Web Tools with sdf query
searchSim <- function(sdf) {
    if(! class(sdf) == "SDFset"){
        stop('reference compound must be a compound of class \"SDFset\"')
    } 
    jobToken <- launchCMTool('Fingerprint Search', sdf, 'Similarity Cutoff'=0.9, 'Max Compounds Returned'=200)
    cids <- as.numeric(result(jobToken))
    if(length(cids) == 0) return(SDFset())
    getIds(cids)
}
