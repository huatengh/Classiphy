#this is a function writing out input files for
#PhyloNet to computer MDC scores
#internal function used by the MDCcount function

write.phylonet.subtrees<-function(phy, subs,pre){
  require(ape)
  require(phangorn)
  require(phytools)
  require(stringr)
  l<-write.tree(phy,file='',tree.names = '')
  i=0;
  cat("Tree ",pre, i,' = ',l,"\n",sep='')
  treelist=c(paste(pre, i,sep=''))
  for (x in subs){
    i=i+1;
    if(length(x)>=length(phy$tip.label)-2){next}
    subphy<-drop.tip(phy,tip = x)
    #cat(length(subphy$tip.label))
    l<-write.tree(subphy,file='',tree.names = '')
    cat("Tree ",pre, i,' = ',l,"\n",sep='')
    treelist<-c(treelist,paste(pre, i,sep=''))
  }
  return(treelist)
}
