# this is a function to calculate the rf distance (z score) between gt and st

rfdistance<-function(repfolder=NULL,sfile=NULL, gfile=NULL){
  require(ape)
  require(phangorn)
  require(phytools)
  if(is.na(repfolder)){
    s<-read.tree(sfile)
    g<-read.tree(gfile)

  }else{
    s<-read.tree(paste(repfolder,"/s_tree.trees",sep=''))
    g<-read.tree(paste(repfolder,"/g_trees.trees",sep=''))
  }
  dsg<-rep(0,length(g))
  for (i in 1:length(g)){
    g[[i]]$tip.label<-sub("_0_0$","",g[[i]]$tip.label,perl =T )
    g[[i]]<-drop.tip(g[[i]],tip="outgroup")
    dsg[i]<-treedist(g[[i]],s)[1]
  }
  return(righttailp(dsg))
}
