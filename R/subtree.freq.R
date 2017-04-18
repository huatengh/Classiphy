#this is a script to calculate the subtree frequency for gts
#take the folder name as input
#could specify the maximum frequency,default is 10% of the gene tree
#i.e., we are not interested in the number of subtrees that appears in more than 10% of the gts (max.freq=0.1)
#return a matrix
#by default also save a output file called subtreefreq.txt in the folder (save.middle.file=T)
subtree.freq<-function(repfolder=NULL,sfile=NULL, gfile=NULL,max.freq=0.1,save.middle.file=T){
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
  outfile=paste(repfolder,'/subtreefreq.txt',sep = '')
  for (i in 1:length(g)){
    g[[i]]$tip.label<-sub("_0_0$","",g[[i]]$tip.label,perl =T )
    g[[i]]<-drop.tip(g[[i]],tip="outgroup")
  }
  sub.tree<-c()
  counts<-c()
  for (gg in g){
    x<-subtrees(gg)
    for(xx in x){
      y<-as.integer(xx$tip.label)
      y<-y[order(y)]
      t<-paste(y,collapse = '_',sep = '')
      l<-which(sub.tree==t)
      if(length(l)==0){
        sub.tree<-c(sub.tree,t)
        counts<-c(counts,1)
      }else{
        counts[l]=counts[l]+1
      }
    }
  }
  gfreq<-matrix(0,nrow = 0,ncol = length(subtrees(s)))
  for(gg in g){
    x<-subtrees(gg)
    gfreq<-rbind(gfreq,sapply(x,function(xx){
      y<-as.integer(xx$tip.label)
      y<-y[order(y)]
      t<-paste(y,collapse = '_',sep = '')
      counts[which(sub.tree==t)]/length(g)}))
  }

  for(i in 1:dim(gfreq)[1]){
    t<-gfreq[i,]
    gfreq[i,]<-t[order(t)]
  }
  if(save.middle.file){
    write.table(gfreq,file=outfile,sep='\t',row.names = F,col.names = F)
  }

  j<-ceiling(max.freq*length(g))
  ks<-apply(gfreq,MARGIN = 1,function(l){sapply(1:j,function(i) sum(l<=i/dim(xx)[1]))})
  ks<-t(ks)

  Q75<-apply(ks,MARGIN = 2,function(i)quantile(i,0.75))
  Q25<-apply(ks,MARGIN = 2,function(i)quantile(i,0.25))
  IQR<-Q75-Q25
  IQR[IQR==0]<-0.5
  k<-matrix(0,nrow = dim(ks)[1],ncol = dim(ks)[2])
  for(i in 1:dim(xx)[1]){
    l<-ks[i,]-Q75
    k[i,]<-as.vector(as.numeric(l/IQR))
  }

  colnames(k)<-paste('subtreescore',1:j,sep='')

  return(k)
}
