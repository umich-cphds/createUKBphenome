# Rscript written by Lars Fritsche to extract / reformat UKB data

options(stringsAsFactors=F)
library("data.table")
library("optparse")
library("intervals")

# data.table with the fields and their descriptions
all.tables <- fread("/net/junglebook/home/larsf/Projects/UKB/data/Fields_in_Current_Merged_Data.txt",header=T)

# data.table with the merged baskets
data.file <- "/net/junglebook/home/larsf/Projects/UKB/data/merged_baskets_20190924.txt"

ukb_data_columns <- scan(data.file,what=character(0),sep="\t",nlines=1,quiet=T)

# merge consecutive numbers to intervals and collapse
getRanges <- function(colnumbers,collapse=T){
	cx <- clusters(colnumbers,1)
	ux <- colnumbers[which(!colnumbers %in% unlist(cx))]
	ix <- sapply(cx,function(x) paste(range(x),collapse="-"))
	out <- c(ux,ix)
	out <- out[order(as.numeric(gsub("\\-.+","",out)))]
	paste(out,collapse=",")
}

reformatUKB <- function(fields,dataCoding=F){
	# print info about fields
	fieldinfo <- all.tables[which(all.tables$'Field ID' %in% fields & !is.na(all.tables$basket)),]

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

	colnumbers <- which(ukb_data_columns %in% selected_columns)
	reformatted <- fread(cmd=paste("cut -f",getRanges(colnumbers),data.file),header=T)
	
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

option_list <- list(
  make_option("--fields", type="character", default="",
    help="comma-separated fields"),
  make_option("--output", type="character", default="",
    help="Full path to output file")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

if(opt$fields != "" | opt$output != ""){
	output <- reformatUKB(as.integer(strsplit(opt$fields,",")[[1]]))
	fwrite(output,opt$output,sep="\t",col.names=T,row.names=F,quote=T)
} else {
	print("Not all parameters present")
}
