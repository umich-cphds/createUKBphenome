harmonizeICD9 <- function(icd){
    if(is.na(icd) | icd == "" | icd == "NA" ) return(NA)
 
    nums <- 3        
    if(grepl("^[A-Z]+",icd)) {
        i0 <- gsub("^([A-Z]+).+","\\1",icd)
        if(i0 == "V") nums <- 2
    } else {
        i0 <- character(0)
    }
    icd <- gsub("^[A-Z]+","",icd)
    
    hasSuffix <- grepl("[A-Za-z]+$",icd)
    if(hasSuffix){
        i3 <- gsub(".+([A-Za-z]+)$","\\1",icd)
        icd <- gsub("[A-Za-z]+$","",icd)
    } else {
        i3 <- character(0)
    }

    if(nchar(icd)>nums & !grepl("\\.",icd)){
        icd <- paste0(substr(icd, 1, nums),".",substr(icd, nums+1, nchar(icd)))
    }
    itemp <- strsplit(icd,"\\.")[[1]]
    i1 <- formatC(as.integer(gsub("^[A-Z]+","",itemp[1])),width=nums,flag=0)

    if(length(itemp)==2){
        i2 <- paste0(".",itemp[2])
    } else {
        i2 <- character(0) 
    }
    paste0(i0,i1,i2)
}
