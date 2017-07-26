#this is a function calculating summary statistics for gene trees

sumstat<-function(repfolder=NULL,sfile=NULL,gfile=NULL,phylonet.path,max.freq=0.1,rf=T,tscore=T,mdcscore=T,subtree=T){
  require(ape)
  if(is.null(gfile)){
    gfile<-paste(repfolder,"/g_trees.trees",sep='')
  }
  g<-read.tree(gfile)
  x<-data.frame(gfile=rep(gfile,length(g)),stringsAsFactors = F)
  if(rf){
    print("calculating rf distances")
    x<-cbind(x,rfdistance(repfolder = repfolder,sfile=sfile,gfile=gfile))
  }
  if(tscore){
    print("calculating triplet score")
    x<-cbind(x,triplet.score(repfolder = repfolder,sfile=sfile,gfile=gfile))
  }
  if(mdcscore){
    print("calculating mdc score")
    x<-cbind(x,MDCcount(repfolder=repfolder,sfile=sfile,gfile=gfile,phylonet.path = phylonet.path))
  }
  if(subtree){
    print("calculating subtree frequencies")
    x<-cbind(x,subtree.freq(repfolder = repfolder,sfile = sfile,gfile=gfile,max.freq=max.freq))
  }
  colnames(x)[2:dim(x)[2]]<-paste("predictor.",colnames(x)[2:dim(x)[2]],sep='')
  return(x)
}
