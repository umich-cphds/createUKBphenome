options(stringsAsFactors=F)
library("data.table",quietly=T)
library("htmltab")

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
	table_url <- paste0("https://biobank.ndph.ox.ac.uk/showcase/list.cgi?it=0&vt=",vt[i])
	tableX <- data.table(htmltab(doc = table_url,which=1))
	tableX$datatype <- datatype[i]
	all.tables[[i]] <- tableX
}

all.tables <- rbindlist(all.tables)

fwrite(all.tables,"./data/FieldListing.txt",sep="\t")
