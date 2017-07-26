# this is a script for getting the training data for DAPC

training.data<-function(repfolder=NULL,sfile=NULL,gfile=NULL,lfile=NULL,phylonet.path,max.freq=0.1,rf=T,tscore=T,mdcscore=T,subtree=T){
  y<-true.distance(repfolder=repfolder,sfile=sfile, gfile=gfile, lfile=lfile)
  x<-sumstat(repfolder=repfolder,sfile=sfile, gfile=gfile,phylonet.path=phylonet.path,max.freq=max.freq,rf=rf,tscore=tscore,mdcscore=mdcscore,subtree=subtree)
  x$model.id<-'ILS'
  x$model.id[y$dSL>0]<-'HGT'
  return(x)
}
