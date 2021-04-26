


test.formatConversions <- function() {
	if(! ChemmineR:::.haveOB()) 
		DEACTIVATED("this test requires ChemmineOB, but not available")
	message("test.formatConversions")

	data(sdfsample)

	smiles = sdf2smiles(sdfsample[1])
	sdfs = smiles2sdf(smiles)
	smiles2 = sdf2smiles(sdfs)
	smiles =gsub("[\r\n]", "", smiles) #remove newlines and convert to simple string
  smiles2 =gsub("[\r\n]", "", smiles2)

	checkEquals(smiles,smiles2)

}

test.genAPDescriptors <- function(){

	DEACTIVATED("removed old version of function")

	data(sdfsample)
	
	for(i in 1:100){
		sdf = sdfsample[[i]]
		desc = genAPDescriptors(sdf) 
		#print(head(desc));

		oldDesc=ChemmineR:::.gen_atom_pair(ChemmineR:::SDF2apcmp(sdf))
                 
		#print(oldDesc);
		compResult = desc==oldDesc
		if(!all(compResult)){
			message("descriptor mismatch")
			print(oldDesc[!compResult])
			print(desc[!compResult])
			message("----------")
			firstFalse = match(FALSE,compResult)
			print(oldDesc[(firstFalse-5):(firstFalse+5)])
			print(desc[(firstFalse-5):(firstFalse+5)])

		}
		checkTrue(all(compResult))
		#checkEqualsNumeric(desc,oldDesc)
	}
}
test.propOB <- function() {
	DEACTIVATED("fails on ubuntu 16.04")
	data(sdfsample)
	p = propOB(sdfsample[1:5])
	#print(p)
	#checkEquals(ncol(p),15)
	checkEquals(nrow(p),5)
   checkEquals(p$MW[2],MW(sdfsample[2])[[1]])

}
test.fingerprintOB <- function(){

	if(! ChemmineR:::.haveOB()) 
		DEACTIVATED("this test requires ChemmineOB, but not available")

	data(sdfsample)
	fp = fingerprintOB(sdfsample[1:5],"FP2")
	checkEquals(fptype(fp),"FP2")
	fpSingle = fingerprintOB(sdfsample[1],"FP2")
	checkEquals(as.character(class(fpSingle)),"FPset")
	checkEqualsNumeric(as.matrix(fpSingle[1]), as.matrix(fp[1]))
}
test.obmolRefs <- function() {
	if(! ChemmineR:::.haveOB()) 
		DEACTIVATED("this test requires ChemmineOB, but not available")

	data(sdfsample)
	obmolRef = obmol(sdfsample[[1]])

	checkEquals(class(obmolRef)[[1]],"_p_OpenBabel__OBMol")

	obmolRefs = obmol(sdfsample)
	checkEquals(class(obmolRefs),"list")
	checkEquals(class(obmolRefs[[2]])[[1]],"_p_OpenBabel__OBMol")
	checkEquals(length(sdfsample),length(obmolRefs))
}
test.smartsSearchOB <- function(){
	if(! ChemmineR:::.haveOB()) 
		DEACTIVATED("this test requires ChemmineOB, but not available")
	data(sdfsample)
	rotableBonds = smartsSearchOB(sdfsample[1:5],"[!$(*#*)&!D1]-!@[!$(*#*)&!D1]",uniqueMatches=FALSE)
	print("rotable bonds: ")
	print(rotableBonds)
	print(sdfid(sdfsample[1:5]))
	checkEquals(as.vector(rotableBonds[1:5]),c(24,20,14,30,10))

}
test.fpSim <- function(){
	data(apset)
	fpset = desc2fp(apset)
	dists = fpSim(fpset[[1]],fpset,top=6)
	checkEqualsNumeric(dists, 
							 c(1.0000000,0.4719101,0.4288499,0.4275229,0.4247423,0.4187380),
							 tolerance = 0.0001)

	for(m in c("tanimoto","euclidean","tversky","dice")){
		sim = ChemmineR:::fpSimOrig(fpset[[1]],fpset,
					method=m,cutoff=0.4,top=6)
		simFast= fpSim(fpset[[1]],fpset,
					method=m,cutoff=0.4,top=6)
		#message("method: ",m)
		#print(sim)
		#print(simFast)
		checkEqualsNumeric(sim,simFast,tolerance=0.00001)
	}

	
}
test.fpSimParameters<- function(){
	data(apset)
	fpset = desc2fp(apset)

	params = genParameters(fpset)

	similarities = fpSim(fpset[[1]],fpset,top=6,parameters=params)
	#print(similarities)

	checkEqualsNumeric(similarities$similarity, 
							 c(1,0.471910112359551,0.428849902534113,0.427522935779817,0.424742268041237,0.418738049713193),
							 tolerance = 0.0001)
	checkEqualsNumeric(similarities$zscore, 
							 c(6.41202501933441,1.6379277432492,1.24865001753133,1.23665382424036,1.21151572069559,1.15723571451518),
							 tolerance = 0.0001)
	checkEqualsNumeric(similarities$evalue, 
							 c(0.00000,6.36243,11.64270,11.84515,12.27769,13.25050),
							 tolerance = 0.0001)
	checkEqualsNumeric(similarities$pvalue, 
							 c(0,0.998274830604939,0.999991217102066,0.999992826756578,0.999995345559029,0.999998240535916),
							 tolerance = 0.0001)

}
test.exactMassOB <- function(){
	if(! ChemmineR:::.haveOB()) 
		DEACTIVATED("this test requires ChemmineOB, but not available")
	data(sdfsample)
	mass = exactMassOB(sdfsample[1:5])
	checkEqualsNumeric(mass,c(456.2009,357.1801,
									  370.1100,461.1733,
									  318.1943),tolerance=0.00001)
}
test.3dCoords <-function(){
	DEACTIVATED("causing timeout on bioc, disabling for now")
	data(sdfsample)
	sdf3d = generate3DCoords(sdfsample[1])

	checkTrue(!any(atomblock(sdf3d)[[1]][,3]==0))
	
}
test.canonicalize <- function(){
	if(! ChemmineR:::.haveOB()) 
		DEACTIVATED("this test requires ChemmineOB, but not available")
	data(sdfsample)
	cansdf = canonicalize(sdfsample[1])

	bb=bondblock(cansdf)[[1]]

	checkEqualsNumeric(bb[1,1:3],c(2,3,1))
	checkEqualsNumeric(bb[2,1:3],c(2,4,1))

}
test.parseV3000 <- function() {

	DEACTIVATED("requires local files")
	sdfset2 = read.SDFset("~/runs/v3000/DrugLike-0_2-3K3K_1.v2k.sdf")
	sdfset3 = read.SDFset("~/runs/v3000/DrugLike-0_2-3K3K_1.sdf")  

	compareSdfVersions = function(v2k,v3k){
		checkEquals(sdfid(v2k),sdfid(v3k))

		#message("v2k: ",nrow(atomblock(v2k)),"x",ncol(atomblock(v2k)))
		#message("v3k: ",nrow(atomblock(v3k)),"x",ncol(atomblock(v3k)))
		#print(head(atomblock(v2k)))
		#print(head(atomblock(v3k)))
		toCompare = c(1:5,7:10) #exclude colum 6
		checkTrue( all(atomblock(v2k)[,toCompare] ==
							atomblock(v3k)[,toCompare]))

#		cmp = bondblock(v2k)[,1:3] == bondblock(v3k)
#		if(! all(cmp)){
#			mismatched = which(cmp==FALSE)
#			print(cmp)
#			print("mismatched: ")
#			print(mismatched)
#			print("data:")
#			print(bondblock(v2k)[mismatched,1:3])
#			print(bondblock(v3k)[mismatched,])
#		}
#		checkTrue( all(bondblock(v2k)[,1:3] == bondblock(v3k)))

		checkTrue( all(datablock(v2k) == datablock(v3k)))
	}
	for(i in seq(along=sdfset2)){
		#if(!(i  %in% c(38,39,89))){ # this differ in acceptable ways
		#	message("testing ",i, " id: ",sdfid(sdfset2[i]))
			compareSdfVersions(sdfset2[[i]],sdfset3[[i]])
		#}
	}


}


test.pubchemPUG <- function(){

	DEACTIVATED("Causes SSL protocol version error on BioC")
	sdf = ChemmineR:::pubchemCidToSDF(c(434,435))
	message("cid to sdf, length: ",length(sdf))	
	checkEquals(length(sdf),2)


	sdf = ChemmineR:::pubchemSmilesSearch("C1CCCCCC1")
	message("smiles search, length: ",length(sdf))
	checkTrue(length(sdf) >= 1)

	data(sdfsample)
	sdf = ChemmineR:::pubchemSDFSearch(sdfsample[1])
	message("sdf search, length: ",length(sdf))
	checkTrue(length(sdf) >= 1)

	# just check that no exception is thrown
	ChemmineR:::pubchemSDF2PNG(sdfsample[1],"test-sample.png")


}

test.largestComponent <- function(){
	DEACTIVATED("just for manual testing")
	print("----------------")
	testSdf = smiles2sdf(c(
								  "O=C(NC1CCCC1)CN(c1cc2OCCOc2cc1)C(=O)CCC(=O)Nc1noc(c1)C	TEST1",
								  "CC.CC(=O)C1=CC2=C(C=C1)SC3=CC=CC=C3N2CCCN(C)C.C(=CC(=O)O)C(=O)O	TEST2",
								  "Cl.CCC1C2CC3C4C5(CC(C2C5O)N3C1O)C6=CC=CC=C6N4C	TEST3"))
	answer =  (c( "O=C(NC1CCCC1)CN(c1cc2OCCOc2cc1)C(=O)CCC(=O)Nc1noc(c1)C",
					 "CC(=O)C1=CC2=C(C=C1)SC3=CC=CC=C3N2CCCN(C)C",
					 "CCC1C2CC3C4C5(CC(C2C5O)N3C1O)C6=CC=CC=C6N4C"))
	names(answer) = c("TEST1","TEST2","TEST3")

	largestComps = largestComponent(testSdf)

	message("result: ")
	print(as.character(sdf2smiles(largestComps)))
	print(answer)

	print(largestComps==answer)

	checkTrue(all(largestComps==answer))
	print("----------------")
}

