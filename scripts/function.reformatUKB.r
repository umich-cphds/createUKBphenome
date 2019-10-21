reformatUKB <- function(fields,dataCoding=F,data.file,field.file){
	require("data.table")

	# print info about fields
	fieldinfo <- fread(field.file,header=T)
	fieldinfo <- fieldinfo[which('Field ID' %in% fields & !is.na(basket)),]

	if(dataCoding){
		require("RCurl")
		require("htmltidy")
		require("XML")
		fieldinfo$URL <- URL <- paste0("https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=",fieldinfo$'Field ID')
		dataCodes <- NULL
		for(i in 1:nrow(fieldinfo)){
			u <- URL[i]
			doc.raw <- getURL(u)
			doc <- tidy_html(doc.raw)
			html <- htmlTreeParse(doc, useInternal = TRUE)
			txt <- xpathApply(html, "//body//text()[not(ancestor::script)][not(ancestor::style)][not(ancestor::noscript)]", xmlValue)
			txt <- unlist(unlist(txt))
			codedData <- which(grepl("Data-Coding",txt))
			if(length(codedData)==1){
				dataCodes <- c(dataCodes,txt[which(grepl("Data-Coding",txt))+1])
			} else {
				dataCodes <- c(dataCodes,NA)
			}
		}
		fieldinfo$'Data Coding' <- paste0("https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=",dataCodes)
		print(data.table(fieldinfo))
	} else {
		print(fieldinfo[,1:4])
	}
	
	# Extract all fields 
	allfields <- fieldinfo$subfield_id[which(fieldinfo$subfield_id != "n/a")]
	
	selected_columns <- c("id",unique(unlist(strsplit(allfields,","))))	
 	print(paste("Reading all entries across",length(selected_columns),"columns"))

	reformatted <- fread(data.file,select=selected_columns)
	
	fieldpattern <- "(.+)\\-.+\\-.+"
	entrypattern <- ".+\\-(.+\\-.+)"

	keep <- which(rowSums(is.na(reformatted[,-1]) | reformatted[,-1] == "") != ncol(reformatted[,-1]))
	reformatted <- reformatted[keep,]

	entries <- names(reformatted)[-1]
	entryNumbers <-	gsub(entrypattern,"\\1",entries)
	fieldCols <- unique(gsub(fieldpattern,"\\1",entries))
	uentries <- unique(entryNumbers)

	# reformat cancer registry (allow multiple entries by person)
	reformatted2 <- list()
	for(e in uentries){
		newEntries <- entries[which(entryNumbers == e)]
		newRows <- reformatted[,c("id",newEntries),with=F]
		setnames(newRows,newEntries,gsub(fieldpattern,"\\1",newEntries))
		missingFields <- fields[which(!fields %in% names(newRows))]
		for(missingField in missingFields){
			newRows[[as.character(missingField)]] <- NA
		}
		reformatted2[[e]] <- newRows[,c("id",fields),with=F]
	}

	reformatted2 <- rbindlist(reformatted2)

	# remove empty entries	
	keep <- which(rowSums(is.na(reformatted2[,-1]) | reformatted2[,-1] == "") != ncol(reformatted2[,-1]))
	reformatted2	<- reformatted2[keep,]
	print(paste(nrow(reformatted2),"lines after filtering empty lines"))
	if(!dataCoding){
		return(reformatted2)
	} else {
		return(list('data'=reformatted2,'info'=data.table(fieldinfo,dataCodes)))
	}
}
