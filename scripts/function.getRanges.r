# merge consecutive integers to intervals and collapse

getRanges <- function(colnumbers,collapse=T){
	require("intervals")
	cx <- clusters(colnumbers,1)
	ux <- colnumbers[which(!colnumbers %in% unlist(cx))]
	ix <- sapply(cx,function(x) paste(range(x),collapse="-"))
	out <- c(ux,ix)
	out <- out[order(as.numeric(gsub("\\-.+","",out)))]
	paste(out,collapse=",")
}
