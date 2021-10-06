
pubchemServerURL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

pubchemCidToSDF = function(cids){

   if(! class(cids) == "numeric")
       stop('reference compound ids must be of class \"numeric\"')
   
   if(length(cids) == 0)
       stop('no compounds to retrieve- input must contain at least one cid')

	finalSDF=NA
   
	batchByIndex(cids,function(indexBatch){
		url = paste(pubchemServerURL,"compound","cid",paste(indexBatch,collapse=","),"SDF",sep="/")
		suppressWarnings(
			if(is.na(finalSDF)){
				finalSDF<<-read.SDFset(readLines(url(url)))
			}else{
				finalSDF<<-c(finalSDF,read.SDFset(readLines(url(url))))
			})
	},100)
	cid(finalSDF) = makeUnique(cid(finalSDF),silent=TRUE)
	finalSDF

}

# PubChem only allow you to request by InChI one at a time, so this is possible 
# but should NOT be recommended by the package. Remeber sleep for some time like half 
# or one second if this is in a for loop. To be polite with API calls. 

# The only useful thing you can get through API is CID., so the return is a vector 
# of CIDs. You are not allowed to request other info by InChI, like SDF, PNG, etc. 

#' Query pubchem by InChI sttrings and return CIDs
#' @description Use PubChem API to get CIDs by InChI sttrings. This function sends 
#' one request per InChI. For courtesy, it is not recommended to parellelize this function. 
#' @param inchis character vector, InChI strings
#' @param verbose logical, show verbose information? 
#' @return a numeric vector of CIDs with names. Successful requests will have empty 
#' names, requests with invalid InChI strings will have name "invalid" and requests 
#' with valid InChI but not found in PubChem will have name "not_found"
#' @export
#' @importFrom jsonlite fromJSON
#' @importFrom stringi stri_replace_all_fixed
#' @examples 
#' # Imagine we have an InChIs as following
#' # first two are valid, thrid has no result, last is invalid 
#' if(interactive()) {
#'     inchis <- c(
#'         "InChI=1S/C15H26O/c1-9(2)11-6-5-10(3)15-8-7-14(4,16)13(15)12(11)15/h9-13,16H,5-8H2,1-4H3/t10-,11+,12-,13+,14+,15-/m1/s1", 
#'         "InChI=1S/C3H8/c1-3-2/h3H2,1-2H3", 
#'         "InChI=1S/C15H20Br2O2/c1-2-12(17)13-7-3-4-8-14-15(19-13)10-11(18-14)6-5-9-16/h3-4,6,9,11-15H,2,7-8,10H2,1H3/t5-,11-,12+,13+,14-,15-/m1/s1",
#'         "InChI=abc"
#'     )
#'     pubchemInchi2cid(inchis)
#' }
pubchemInchi2cid <- function(inchis, verbose = TRUE) {
    stopifnot(is.character(inchis) && length(inchis) > 0)
    if(length(inchis) == 1) stopifnot(nchar(inchis) > 0)
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    
    # InChI strings need to be escaped, all possible special characters I am aware 
    # are listed below, add more if I missed some one. minus/dash does not need to be 
    # escaped -LZ
    inchi_escaped <- stringi::stri_replace_all_fixed(
        str = inchi,
        pattern = c("/", "[", "]", "&", ",", "(", ")", "=", "+"),
        replacement = c("%2F", "%5B", "%5D", "%26", "%2C", "%28", "%29", "%3D", "%2B"),
        vectorize_all = FALSE
    )
    
    inchi_request_base_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON?inchi="
    
    cids <- c()
    inchi_total <- length(inchi_escaped)
    for(i in seq_along(inchi_escaped)) {
        Sys.sleep(0.5)
        cat("requesting for InChI", i, "/", inchi_total, "\n")
        res <- tryCatch(
            {
                req_url <- paste0(inchi_request_base_url, inchi_escaped[i])
                jsonlite::fromJSON(url(req_url))$IdentifierList$CID
            },
            error = function(e) {
                return(c("invalid" = 0))
            },
            warning = function(w) {
                if(grep("400 Bad Request", w$message) == 1) cat("Invalid InChI, skip current one, No.", i, "\n", sep = "")
                return(c("invalid" = 0))
            }
        )
        if(res == 0 && is.null(names(res))) {
            res <- c("not_found" = res)
            cat("InChI not found in PubChem, No.", i, "\n", sep = "")
        }
        cids <- c(cids, res)
    }
    cids
}

#' Query pubchem by InChIKeys sttrings and get SDF back
#' @description Use PubChem API to get CIDs by InChIKeys
#' @param inchikeys character vector, InChIKey strings
#' @return a list of 2 items. the first item "sdf_set" is a `SDFset` object. It 
#' contains all queried and successful SDF infomation. The second item "sdf_index"
#' is a named numeric vector. It records whether all input InChIKeys have successful 
#' returns in the  `SDFset` object. If so, a non-zero value is returned as the index of
#' where it exists in the `SDFset` object, if not, 0 is returned.
#' @export
#'
#' @examples
#' inchikeys <- c(
#' "ZFUYDSOHVJVQNB-FZERPYLPSA-N", 
#' "KONGRWVLXLWGDV-BYGOPZEFSA-N", 
#' "AANKDJLVHZQCFG-WLIQWNBFSA-N",
#' "SNFRINMTRPQQLE-JQWAAABSSA-N"
#' )
#' pubchemInchikey2sdf(inchikeys)
pubchemInchikey2sdf <- function(inchikeys){
    stopifnot(is.character(inchikeys) && length(inchikeys) > 0)
    if(length(inchikeys) == 1) stopifnot(nchar(inchikeys) > 0)
    if(any(duplicated(inchikeys))) stop("Make sure all InChI strings are unique.")
    
    finalSDF <- NA
    batchByIndex(inchikeys,function(indexBatch){
        url <- paste(pubchemServerURL,"compound","inchikey",paste(indexBatch,collapse=","),"SDF",sep="/")
        suppressWarnings(
            if(is.na(finalSDF)){
                finalSDF<<-read.SDFset(readLines(url(url)))
            }else{
                finalSDF<<-c(finalSDF,read.SDFset(readLines(url(url))))
            })
    },100)
    # handle missing not found inchikeys
    finalSDF_list <- SDFset2list(finalSDF)
    # get inchikey and cid from SDF as dataframe
    sdf_ma <- as.data.frame(t(data.frame(lapply(finalSDF_list, function(x) x$datablock[c("PUBCHEM_COMPOUND_CID", "PUBCHEM_IUPAC_INCHIKEY")]))))
    # match queried inchikeys to returned inchikeys and record position
    sdf_index <- unlist(lapply(inchikeys, function(x){
        index <- which(sdf_ma$PUBCHEM_IUPAC_INCHIKEY %in% x)
        if(length(index) == 0) 0 else index
    }))
    names(sdf_index) <- inchikeys
    # rename SDF by CIDs
    cid(finalSDF) <- sdf_ma$PUBCHEM_COMPOUND_CID
    
    list(
        sdf_set = finalSDF,
        sdf_index = sdf_index
    )
}
						
pubchemSmilesSearch = function(smiles){
	if(class(smiles) == "SMIset")
		smiles = as.character(smiles)
   if(! class(smiles) == "character"){
      stop('reference compound must be a smiles string of class \"character\"')
   } 
	url = paste(pubchemServerURL,"compound","fastsimilarity_2d","smiles",smiles,"SDF",sep="/")
	read.SDFset(readLines(suppressWarnings(url(url))))
}

pubchemSDFSearch = function(sdf){
   if(! class(sdf) == "SDFset"){
        stop('reference compound must be a compound of class \"SDFset\"')
   } 
   
	url = paste(pubchemServerURL,"compound","fastsimilarity_2d","sdf","SDF",sep="/")
	#message("url: ",url)

	sdfStr = paste(as(sdf[[1]],"character"),collapse="\r\n")
	read.SDFset(readLines(textConnection(rawToChar(postForm(url,sdf = sdfStr)))))
}
pubchemSDF2PNG = function(sdf,outputFile){
	if(! class(sdf) == "SDFset"){
        stop('reference compound must be a compound of class \"SDFset\"')
   } 
   
	url = paste(pubchemServerURL,"compound","fastsimilarity_2d","sdf","PNG",sep="/")
	#message("url: ",url)

	sdfStr = paste(as(sdf[[1]],"character"),collapse="\r\n")
	pngData = postForm(url,sdf=sdfStr)
	writeBin(pngData[1:length(pngData)],outputFile)
}
pubchemName2CID = function(name){
	
	url = paste(pubchemServerURL,"compound","name",name,"cids","txt",sep="/")
	#message("url: ",url)

	tryCatch(
		readLines(suppressWarnings(url(url))),
		error=function(e) NA,
		warning=function(e) NA)
}
