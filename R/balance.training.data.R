# this is a short function for balancing the training data

balance.training.data<-function(x,repncolumnname="repn"){
  if(! "model.id" %in% colnames(x)){
    stop("no model.id column in the input data frame")
  }
  if(! repncolumnname %in% colnames(x)){
    stop("need a column indicating which rows are from the same species tree")
  }
  l<-c()
  model<-as.character(x$model.id)
  for( r in unique(x[,repncolumnname])){
    ils<-which(xx[,repncolumnname]==r & model=="ILS")
    hgt<-which(xx[,repncolumnname]==r & model=="HGT")
    if(length(ils)>length(hgt)){
      l<-c(l,ils[(length(ils)-length(hgt)):length(ils)])
    }else{
      l<-c(l,hgt[(length(hgt)-length(ils)):length(hgt)])
    }
  }

  x<-x[-l,]
  x$model.id<-as.factor(x$model.id)
  return(x)
}
