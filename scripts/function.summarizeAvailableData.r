library(data.table)

withdrawn <- list.files("./data","w.+.csv",full.name=T)

all.tables <- fread("./data/FieldListing.txt",sep="\t")
all.tables[,`:=`(subfield_id=as.character(NA),basket=as.character(NA),basketfile=as.character(NA))]

baskets <- readLines("./data/baskets.txt")

names(baskets) <- gsub(".tab","",basename(baskets))

baskets <- baskets[order(names(baskets))]

counts <- NULL

for(b in 1:length(baskets)){
	btemp <- names(fread(baskets[b],nrow=1,header=T))
	bcols <- data.table('field'=gsub("[f\\.]{0,2}([0-9]+)[\\.\\-].+","\\1",btemp),'subfield'=btemp)

	bcols <- split(bcols,bcols$field)
	
	for(bcol in names(bcols)){
		available <- which(all.tables$'Field ID' == bcol)
		counts <- c(counts,length(available))
		if(length(available)>0) {
			all.tables$subfield_id[available] <- paste(bcols[[bcol]]$subfield,collapse=",")
			all.tables$basket[available] <- names(baskets)[b]
			all.tables$basketfile[available] <- baskets[b]
		}
	}
}

all.tables <- all.tables[order(as.integer(all.tables$'Field ID')),]

fwrite(all.tables,"./data/Fields_in_Available_Data.txt",
	sep="\t",quote=T,row.names=F,col.names=T,na="n/a")
	

# merge baskets in one huge data file
# Overwrite the older data with the new data

for(b in 1:length(baskets)){

	btemp <- fread(baskets[b],header=T,quote="")

	if(grepl("^f",names(btemp)[1])){
		eid <- "f.eid"
		pattern <- "^f\\.([0-9]+)\\.([0-9]+)\\.([0-9]+)"
	} else {
		eid <- "eid"
		pattern <- "(^[0-9]+)\\-([0-9]+)\\.([0-9]+)"
	}

	newnames <- c("id",
		gsub(pattern,"\\1-\\2-\\3",names(btemp)[-1]))
	setnames(btemp,names(btemp),newnames)

	if(b == 1){
		ukb_data <- btemp
	} else {
		# removeCols <- names(btemp)[which(names(btemp) %in% names(ukb_data) & names(btemp) != "id"),]
		keepCols <- names(ukb_data)[which(!names(ukb_data) %in% names(btemp) | names(ukb_data) == "id")]
		
		ukb_data <- ukb_data[,keepCols,with=F]
		ukb_data <- merge(ukb_data,btemp,by="id")
	}
	print(dim(ukb_data))
	rm("btemp")
	gc()
}

# Remove individuals who withdrew
wids <- readLines(sort(withdrawn,decreasing=T)[1])

ukb_data <- ukb_data[!id %in% wids,]

fwrite(ukb_data,"./data/merged_baskets.txt",sep="\t",quote=F)
