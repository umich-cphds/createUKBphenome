options(stringsAsFactors=F)
library("data.table",quietly=T)
library("XML")
library("RCurl")
library("rlist")

datatype <- c('Integer',
	'Categorical (single)',
	'Categorical (multiple)',
	'Continuous',
	'Text',
	'Date',
	'Time',
	'Compound',
	'unknown',
	'unknown',	
	'unknown')

vt <- c("11","21","22","31","41","51","61","101")

all.tables <- list()

for(i in 1:length(vt)){
	tables <- readHTMLTable(paste0("http://biobank.ndph.ox.ac.uk/showcase/list.cgi?it=0&vt=",vt[i]))
	tables <- list.clean(tables, fun = is.null, recursive = FALSE)
	tableX <- tables[[1]]
	if(ncol(tableX) == 4){
		tableX <- tableX[,-4]
		colnames(tableX)[1:3] <- c("Field ID","Description","Category")
	}
	for(j in 1:3){
		tableX[[j]] <- gsub("[  ]+$|^[  ]+","",tableX[[j]])
	}
	tableX$datatype <- datatype[i]
	all.tables[[i]] <- tableX
}

all.tables <- rbindlist(all.tables)

fwrite(all.tables,"./data/FieldListing.txt",sep="\t")
