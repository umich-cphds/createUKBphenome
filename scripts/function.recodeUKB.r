recodeUKB <- function(ukb_data,datainfo){
	require("RCurl")
	require("XML")
	for(i in 1:nrow(datainfo)){
		fieldID <- as.character(datainfo[i,'Field ID'])
		newColumnName <-  as.character(datainfo[i,'Description'])
		URL <- datainfo[i,'Data Coding']	
		doc.raw <- getURL(URL)
		coding <- readHTMLTable(doc.raw)[[2]]	
		ukb_data[,newColmn:=ukb_data[,fieldID,with=F]]	
		ukb_data[,newColmn:=factor(newColmn,levels=coding$Coding,labels=coding$Meaning)]
		setnames(ukb_data,"newColmn",newColumnName)
	}
	return(ukb_data)
}

