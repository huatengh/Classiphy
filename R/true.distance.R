# this is a function calculate various true distances (i.e., not for empirical testing data just for training data):
# number of HGT events
# rf distance between st and lt-- whether HGT -induced topo change occurred
# rf distance between lt and gt-- the relative amount of coalescent variances
# by default saves a middle file truedistance.txt in the folder

true.distance<-function(repfolder=NULL,sfile=NULL,gfile=NULL,lfile=NULL,save.middle.file=T){
  require(ape)
  require(phangorn)
  require(phytools)
  if(is.null(repfolder)){
    l<-read.tree(lfile)
    s<-read.tree(sfile)
    g<-read.tree(gfile)

  }else{
    l<-read.tree(paste(repfolder,"/l_trees.trees",sep=''))
    s<-read.tree(paste(repfolder,"/s_tree.trees",sep=''))
    g<-read.tree(paste(repfolder,"/g_trees.trees",sep=''))
  }
  for (i in 1:length(l)){
    l[[i]]$tip.label<-sub("_0$","",l[[i]]$tip.label,perl =T )
    g[[i]]$tip.label<-gsub("_0_0$","",g[[i]]$tip.label,perl =T )
  }
  dis<-sapply(1:length(l),function(x){
    ll<-l[[x]]
    gg<-g[[x]]
    t<-ll$tip.label[grep("Rtransf",ll$tip.label)]
    if(length(t)>0){
      ll<-drop.tip(ll,t)
    }
    dsl<-treedist(ll,s)
    dlg<-treedist(ll,gg)
    c(length(t),dsl[1],dlg[1])
  })
  result=as.data.frame(as.matrix(t(dis)))
  colnames(result)=c("nHGT","dSL","dLG")
  outfile=paste(repfolder,'/truedistance.txt',sep='')
  if(save.middle.file){
    write.table(result,file=outfile,quote = F,row.names = F,sep = "\t")
  }
  return(result)
}
