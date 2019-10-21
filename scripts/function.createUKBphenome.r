# Rscript written by Lars Fritsche to extract / reformat UKB data
options(stringsAsFactors=F)
library(data.table)
library(tidyr)
library(parallel)

## list all TAB-delimited baskets in text file here (one basket per line): e.g. ukb#####.tab
baskets <- readLines("./data/baskets.txt")

# speed up things for data.table by using multiple threads
setDTthreads(detectCores()/2)

today <- format(Sys.Date(), format =  "%Y%m%d")

# function to get ranges from list of integers
source("./scripts/function.getRanges.r")
source("./scripts/function.expandPhecodes.r")
source("./scripts/function.simpleCap.r")
source("./scripts/function.harmonizeICD9.r")
# source("./scripts/function.harmonizeICD10.r")

dir.create("./results")

# read PheWAS map (downloaded from https://phewascatalog.org/phecodes; selected "export all" top right corner)
icd9map <- fread("./data/phecode_icd9_rolled.csv",colClasses="character")

# read PheWAS map (downloaded from https://phewascatalog.org/phecodes_icd10; selected "export all" top right corner)
icd10map <- fread("./data/phecode_icd10.csv",colClasses="character")

# collected phecode information by combining tables from https://github.com/PheWAS/ and the mapping tables
if(!file.exists("./PheWAS")) system("git clone git@github.com:PheWAS/PheWAS.git")

# Phecode information
load(file="./PheWAS/data/pheinfo.rda")
pheinfo <- data.table(pheinfo)
pheinfo <- pheinfo[,.(phecode,description,groupnum,group,color)]

# Gender rules
load(file="./PheWAS/data/gender_restriction.rda")
gender_restriction <- data.table(gender_restriction)
gender_restriction <- rbind(
	data.table('phecode'=gender_restriction[male_only == F & female_only == F,phecode],'sex'="Both"),
	data.table('phecode'=gender_restriction[male_only == T,phecode],'sex'="Male"),
	data.table('phecode'=gender_restriction[female_only == T,phecode],'sex'="Female"))

# Exclusion / roll up rules
pheinfo2 <- unique(fread("./data/phecode_icd9_rolled.csv",colClasses="character",
	select=c("PheCode","Excl. Phecodes","Excl. Phenotypes","Rollup","Leaf"),
	col.names=c("phecode","phecode_exclude_range","phecode_exclude_phenotypes","rollup","leaf")))

# create one data.table with all criteria
pheinfo <- merge(pheinfo,pheinfo2,by="phecode")
pheinfo <- merge(pheinfo,gender_restriction,by="phecode")
pheinfo <- pheinfo[,c("phecode","description","sex","rollup","leaf","groupnum","group","color","phecode_exclude_range","phecode_exclude_phenotypes")]

# add phecode information that's missing (collected form earlier versions)
pheinfoOLD <- fread("./data/Phecode_Definitions_FullTable_Modified.txt",colClasses="character")
pheinfo <- rbind(pheinfo,pheinfoOLD[!phecode %in% pheinfo$phecode,])

# Phecode that should not be rolled up
norollup <- pheinfo$phecode[which(pheinfo$rollup == 0)]

# map ICD9 codes to PheCodes
# read UKB coding
ICD9codes <- fread("./data/coding87.tsv")
ICD9codes[!grepl("Block",coding),ICD9:=sapply(coding,harmonizeICD9)]
codeICD9 <- ICD9codes[,sort(unique(ICD9))]

mappedICD9Codes <- NULL
icd9map_new <- list()
for(i in 1:nrow(icd9map)){
	mapped9 <- grep(icd9map$ICD9[i],codeICD9)
	if(length(mapped9) > 0) {
		mappedICD9Codes <- unique(c(mappedICD9Codes,codeICD9[mapped9]))
		icd9map_new[[icd9map$ICD9[i]]] <- data.table(
			'phecode'=icd9map$PheCode[i],
			'ICD9'=codeICD9[mapped9])
	}
}
icd9key <- rbindlist(icd9map_new)
icd9unmapped <- codeICD9[!codeICD9 %in% mappedICD9Codes]

# roll up PheWAS codes
icd9key$added <- "original"
pcodes <- unique(c(icd9key$phecode,gsub("\\..+","",icd9key$phecode),gsub("(\\..).","\\1",icd9key$phecode)))
pcodes <- sort(pcodes)

for(p in 1:length(pcodes)){
    if(grepl("\\.",pcodes[p])) {
        iSub <- which(icd9key$phecode %in% pcodes[which(grepl(paste0("^",pcodes[p]),pcodes)
        				& nchar(pcodes) > nchar(pcodes[p]))]
        				& !icd9key$phecode %in% norollup)
    } else {
        iSub <- which(icd9key$phecode %in% pcodes[grep(paste0("^",pcodes[p],"\\."),pcodes)]
        				& !icd9key$phecode %in% norollup)
    }
    if(length(iSub) == 0) next
    iTop <- icd9key$ICD9[which(icd9key$phecode == pcodes[p])]
    addTop <- which(icd9key$ICD9 %in% unique(icd9key$ICD9[iSub]) & ! icd9key$ICD9 %in% iTop)
    if(length(addTop) == 0) next

    addKey <- icd9key[addTop,]
    addKey$phecode <- pcodes[p]
    addKey$added <- "rolled up PheWAS code"
    icd9key <- rbind(icd9key,addKey)
}

## Add ICD code description and phecode description
icd9key <- unique(icd9key)
print("Rollup of phewas codes (ICD9 code)")
print(table(icd9key$added))

icd9key <- merge(icd9key,ICD9codes[,.(ICD9,meaning,node_id,parent_id,selectable)],by="ICD9")
icd9key <- merge(icd9key,pheinfo,by="phecode")
icd9key <- icd9key[,c("ICD9","meaning","node_id","parent_id","selectable","phecode","description",
	"group","groupnum","added","sex","rollup","leaf")]


# map ICD10 codes to PheCodes
# read UKB coding
ICD10codes <- fread("./data/coding19.tsv")
ICD10codes <- ICD10codes[!grepl("Block",coding),]
ICD10codes[,ICD10category:=gsub("([A-Z][0-9]{2}).+","\\1",coding)]
ICD10codes[,ICD10suffix:=gsub("[A-Z].+","",gsub("^[A-Z][0-9]{2}","",coding))]
ICD10codes[,ICD10:=paste0(ICD10category,ifelse(ICD10suffix == "","","."),ICD10suffix)]
ICD10codes <- ICD10codes[,c("ICD10category","ICD10suffix"):=NULL]

codeICD10 <- ICD10codes[,sort(unique(ICD10))]
mappedICD10Codes <- NULL
icd10map_new <- list()
for(i in 1:nrow(icd10map)){
	mapped10 <- grep(icd10map$ICD10[i],codeICD10)
	if(length(mapped10) > 0) {
		mappedICD10Codes <- unique(c(mappedICD10Codes,codeICD10[mapped10]))
		icd10map_new[[icd10map$ICD10[i]]] <- data.table(
			'phecode'=icd10map$PheCode[i],
			'ICD10'=codeICD10[mapped10])
	}
}
icd10key <- rbindlist(icd10map_new)
icd10unmapped <- codeICD10[!codeICD10 %in% mappedICD10Codes]

#### roll up phewas codes
icd10key$added <- "original"
pcodes <- unique(c(icd10key$phecode,gsub("\\..+","",icd10key$phecode),gsub("(\\..).","\\1",icd10key$phecode)))
pcodes <- sort(pcodes)
for(p in 1:length(pcodes)){
    if(grepl("\\.",pcodes[p])) {
        iSub <- which(icd10key$phecode %in% pcodes[which(grepl(paste0("^",pcodes[p]),pcodes)
        				& nchar(pcodes) > nchar(pcodes[p]))]
        				& !icd10key$phecode %in% norollup)
    } else {
        iSub <- which(icd10key$phecode %in% pcodes[grep(paste0("^",pcodes[p],"\\."),pcodes)]
        				& !icd10key$phecode %in% norollup)
    }
    if(length(iSub) == 0) next
    iTop <- icd10key$ICD10[which(icd10key$phecode == pcodes[p])]
    addTop <- which(icd10key$ICD10 %in% unique(icd10key$ICD10[iSub]) & ! icd10key$ICD10 %in% iTop)
    if(length(addTop) == 0) next

    addKey <- icd10key[addTop,]
    addKey$phecode <- pcodes[p]
    addKey$added <- "rolled up PheWAS code"
    icd10key <- rbind(icd10key,addKey)
}

icd10key <- unique(icd10key)

print("Rollup of phewas codes (ICD10 code)")
print(table(icd10key$added))

icd10key <- merge(icd10key,ICD10codes[,.(ICD10,meaning,node_id,parent_id,selectable)],by="ICD10")
icd10key <- merge(icd10key,pheinfo,by="phecode")
icd10key <- icd10key[,c("ICD10","meaning","node_id","parent_id","selectable","phecode","description",
	"group","groupnum","added","sex","rollup","leaf")]

print(paste0("Mapping tables created and stored in ./results/UKB_PHENOME_ICD*_PHECODE_MAP_",today,".txt"))

fwrite(icd10key,paste0("./results/UKB_PHENOME_ICD10_PHECODE_MAP_",today,".txt"),sep="\t",quote=T)
fwrite(icd9key,paste0("./results/UKB_PHENOME_ICD9_PHECODE_MAP_",today,".txt"),sep="\t",quote=T)

### Create Phenome

# get Field Listings and store in file `./data/FieldListing.txt`
source("./scripts/function.collectFields.r")
field.file <- "./data/FieldListing.txt"

# load all baskets:
names(baskets) <- basename(baskets)

# get all field ids
fields <- fread(field.file)

# ICD fields without addendums or dates
fieldids9 <- fields[grepl("ICD9",Description) & !grepl("addendum|date",Description,ignore.case=T),`Field ID`]
fieldids10 <- fields[grepl("ICD10",Description) & !grepl("addendum|date",Description,ignore.case=T),`Field ID`]

# Collapse ICD9 and ICD10 separately
# Also collect genetic sex and self-reported gender:
ICD9 <- ICD10 <- SEXGENDER <-  list()

for(basket in baskets){
	ukb_data_columns <- scan(basket,what=character(0),sep="\t",nlines=1,quiet=T)

	extractICD9columns <- grep(paste("^f\\.",fieldids9,"\\.",sep="",collapse="|"),ukb_data_columns)
	if(length(extractICD9columns) > 0) {
		ICD9columns <- getRanges(c(1,extractICD9columns))
		basketData <- fread(cmd=paste("cut -f",ICD9columns,basket),header=T)
		ICD9[[basket]] <- na.omit(gather(basketData, "field", "ICD9", -f.eid))
	}

	extractICD10columns <- grep(paste("^f\\.",fieldids10,"\\.",sep="",collapse="|"),ukb_data_columns)
	if(length(extractICD10columns) > 0) {
		ICD10columns <- getRanges(c(1,extractICD10columns))
		basketData <- fread(cmd=paste("cut -f",ICD10columns,basket),header=T)
		ICD10[[basket]] <- na.omit(gather(basketData, "field", "ICD10", -f.eid))
	}

	extractSexGenderColumns <- grep(paste("^f\\.",c("31","22001"),"\\.",sep="",collapse="|"),ukb_data_columns)
	if(length(extractSexGenderColumns) == 2) {
		SexGenderColumns <- getRanges(c(1,extractSexGenderColumns))
		SEXGENDER[[basket]]  <- fread(cmd=paste("cut -f",SexGenderColumns,basket),header=T)
	}
}

# Process sex/gender information: Female = 0; Male = 1
# Only keep samples where sex == gender; unclear why sex might differ from gender (gender identity, bone marrow transplan, sample swap)
SEXGENDER <- rbindlist(SEXGENDER)
setnames(SEXGENDER,c("f.31.0.0","f.22001.0.0"),c("Sex","GeneticSex"))
SEXGENDER[,Sex:=ifelse(Sex != GeneticSex,NA,Sex)]
SEXGENDER <- na.omit(unique(SEXGENDER))
females <- SEXGENDER[Sex==0,f.eid]
males <- SEXGENDER[Sex==1,f.eid]

ICD9 <- unique(rbindlist(ICD9)[,field:=NULL])
ICD10 <- unique(rbindlist(ICD10)[,field:=NULL])

# Cleanup ICD9 codes; i.e. add dot and remove suffixes
ICD9[,ICD9:=sapply(ICD9,harmonizeICD9)]

# Cleanup ICD10 codes; i.e. add dot and remove suffixes
ICD10[,ICD10category:=gsub("([A-Z][0-9]{2}).+","\\1",ICD10)]
ICD10[,ICD10suffix:=gsub("[A-Z].+","",gsub("^[A-Z][0-9]{2}","",ICD10))]
ICD10[,ICD10:=paste0(ICD10category,ifelse(ICD10suffix == "","","."),ICD10suffix)]
ICD10 <- ICD10[,c("ICD10category","ICD10suffix"):=NULL]

# Merge with phewas map then remove ICD9 codes
phecode1 <- merge(ICD9,icd9key[phecode != "",.(phecode,ICD9)],by="ICD9",allow.cartesian=T)
phecode1[,ICD9:=NULL]

# Merge with phewas map then remove ICD10 codes
phecode2 <- merge(ICD10,icd10key[phecode != "",.(phecode,ICD10)],by="ICD10",allow.cartesian=T)
phecode2[,ICD10:=NULL]

# combine ICD9 and ICD10-based phecodes
phenotypeData <- unique(rbind(phecode1,phecode2))
phecodes <- unique(phenotypeData$phecode)
phecodes <- phecodes[order(as.numeric(phecodes))]

# only keep samples with valid sex (see above)
sampleNames <- SEXGENDER$f.eid

# create two empty data.frames (one with exclusion criteria applied, one without apply exclusion criteria); add X to phecodes
phenoOut <- phenoOut0 <- data.frame('IID'=sampleNames,matrix(NA,ncol=length(phecodes),nrow=length(sampleNames)))
colnames(phenoOut) <- c("IID",paste0("X",phecodes))

# collect inclusions
inclusions <- phenotypeData[,c("f.eid","phecode")]
inclusions <- split(inclusions$f.eid,inclusions$phecode)

# collect exclusions (at this point identical to inclusions)
exclusions <- phenotypeData[,c("f.eid","phecode")]
exclusions <- split(exclusions$f.eid,exclusions$phecode)

# prepare summary table
pheinfo2 <- data.frame(pheinfo[phecode %in% phecodes,c("phecode","description","group","groupnum","color","sex","phecode_exclude_range")],'ncases'=0,'ncontrols'=0)

# collect cases with potentially wrong sex
notFemale <- character(0)
notMale <- character(0)

# create case control studies, one phecode at a time
for(p in 1:nrow(pheinfo2)){
    phecode_remove <- ""

	# collect phecodes to include from controls
    exclude_phecodes <- pheinfo2$phecode[p]
    if(pheinfo2$phecode_exclude_range[p] != ""){
        phecode_remove <- unlist(strsplit(gsub(" ","",pheinfo2$phecode_exclude_range[p]),",")[[1]])
        exclude_phecodes <- c(exclude_phecodes,unlist(sapply(phecode_remove,function(x) expandPhecodes(x,T),USE.NAMES=F)))
    }
    exclude_phecodes <- unique(exclude_phecodes[which(exclude_phecodes %in% pheinfo2$phecode)])

    # non-sex specific traits
    if(pheinfo2$sex[p] == "Both"){
        phenoOut[which(!sampleNames %in% unique(unlist(exclusions[exclude_phecodes]))),p+1] <- 0
        phenoOut[which(sampleNames %in% unlist(inclusions[pheinfo2$phecode[p]])),p+1] <- 1
        phenoOut0[which(!sampleNames %in% unlist(inclusions[pheinfo2$phecode[p]])),p+1] <- 0
        phenoOut0[which(sampleNames %in% unlist(inclusions[pheinfo2$phecode[p]])),p+1] <- 1
    } else if (pheinfo2$sex[p] == "Female"){
	    # female-specific traits
        phenoOut[which(!sampleNames %in% unique(unlist(exclusions[exclude_phecodes])) & sampleNames %in% females),p+1] <- 0
        phenoOut[which(sampleNames %in% unlist(inclusions[pheinfo2$phecode[p]]) & sampleNames %in% females),p+1] <- 1
        phenoOut0[which(!sampleNames %in% unlist(inclusions[pheinfo2$phecode[p]]) & sampleNames %in% females),p+1] <- 0
        phenoOut0[which(sampleNames %in% unlist(inclusions[pheinfo2$phecode[p]]) & sampleNames %in% females),p+1] <- 1
        notFemale <- unique(c(notFemale,sampleNames[which(sampleNames %in% unlist(inclusions[pheinfo2$phecode[p]]) & !sampleNames %in% females)]))
    } else if (pheinfo2$sex[p] == "Male"){
	    # male-specific traits
        phenoOut[which(!sampleNames %in% unique(unlist(exclusions[exclude_phecodes])) & sampleNames %in% males),p+1] <- 0
        phenoOut[which(sampleNames %in% unlist(inclusions[pheinfo2$phecode[p]]) & sampleNames %in% males),p+1] <- 1
        phenoOut0[which(!sampleNames %in% unlist(inclusions[pheinfo2$phecode[p]]) & sampleNames %in% males),p+1] <- 0
        phenoOut0[which(sampleNames %in% unlist(inclusions[pheinfo2$phecode[p]]) & sampleNames %in% males),p+1] <- 1
        notMale <- unique(c(notMale,sampleNames[which(sampleNames %in% unlist(inclusions[pheinfo2$phecode[p]]) & !sampleNames %in% males)]))
    }

    # get sample sizes
    pheinfo2$ncontrols[p] <- length(which(phenoOut[,p+1] == 0))
    pheinfo2$ncases[p] <- length(which(phenoOut[,p+1] == 1))
    print(p)
}

print(paste0("Phenome created and stored in ./results/UKB_PHENOME_*",today,".txt"))

# summary
fwrite(pheinfo2,paste0("./results/UKB_PHENOME_DESCRIPTION_",today,".txt"),sep="\t",row.names=F,col.names=T,quote=T)

# phenome
fwrite(phenoOut,paste0("./results/UKB_PHENOME_",today,".txt"),sep="\t",row.names=F,col.names=T,quote=T)

# phenome without exclusion criteria (sex filter was applied)
fwrite(phenoOut0,paste0("./results/UKB_PHENOME_NO_EXCLUSIONS_",today,".txt"),sep="\t",row.names=F,col.names=T,quote=T)
