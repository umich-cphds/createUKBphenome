aggregateDT <- function(dataIn,field,fieldname,FUN,id="id"){
	out <- data.table(aggregate(x = dataIn[[field]],
                               by = list('id'=dataIn[[id]]),
                               FUN = FUN,na.rm=TRUE, na.action=NULL))
	setnames(out,c("id","x"),c(id,fieldname))
	return(out)
}
