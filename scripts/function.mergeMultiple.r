require("data.table")
mymerge <- function(x,y) merge(x,y,by="id",all=TRUE)
mergeMultiple <- function(dtlist) {
	if(!is.list(dtlist)) stop("Not a list of data.tables / data.frames")
	Reduce(mymerge,dtlist)
}

