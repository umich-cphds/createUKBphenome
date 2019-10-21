expandPhecodes <- function(x,addIntegers=T){
    if(is.na(x) | x == "") return(NA)
    if(grepl("\\-",x)){
        # split range
        # character prefix
        i1 <- strsplit(x,"-")[[1]]
        
        # numeric length of digits before "."
        nprefix <- max(nchar(gsub("\\..+","",i1)))
        # numbers of digits
        ndigits <- max(c(nchar(gsub("^[0-9]+[\\.]{0,1}","",i1)),0))
        # add "." to length of formatted number if present
        addDot <- max(as.numeric(grepl("\\.",i1)))
        # create sequence of ICD codes
        seq1 <- seq(as.numeric(i1[1]),as.numeric(i1[2]),(1/10^ndigits))
        # format sequence to match intput
        seq1 <- formatC(seq1, format='f', digits=ndigits,width=nprefix+ndigits+addDot,flag=0)
        # add integers if within range 
        if(addIntegers) seq1 <- unique(sort(c(seq1,gsub("\\..+","",seq1[which(round(as.numeric(seq1)) == as.numeric(seq1))]))))
        
        if(ndigits == 2){
            seq2 <- seq(as.numeric(i1[1]),as.numeric(i1[2]),(1/10^(ndigits-1)))
            seq2 <- formatC(seq2, format='f', digits=ndigits-1,width=nprefix+ndigits+addDot-1,flag=0)
            seq1 <- unique(sort(c(seq1,seq2)))
        }      
        return(seq1)
    } else {
        return(x)
    }
}
