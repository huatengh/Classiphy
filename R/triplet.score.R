# this is a function calculating the triplet score for gts
# have the option to save middle file (large RData file, default is false)
# depends on the number of tips, could need large memory

triplet.score<-function(repfolder=NULL,sfile=NULL, gfile=NULL,save.middle.file=F){
  require(ape)
  require(phangorn)
  require(phytools)
  if(is.null(repfolder)){
    s<-read.tree(sfile)
    g<-read.tree(gfile)
    repfolder=sub(basename(gfile),'',gfile,perl = T)

  }else{
    s<-read.tree(paste(repfolder,"/s_tree.trees",sep=''))
    g<-read.tree(paste(repfolder,"/g_trees.trees",sep=''))
  }

  phy_parents<-function(phy){
    nodelist<-c()
    nodeparent<-list()
    for (node in (length(phy$tip.label)+2):max(phy$edge)){
      p<-c(node)
      while (p[length(p)]>length(phy$tip.label)+1){
        p<-c(p,phy$edge[phy$edge[,2]==p[length(p)],1])
        if (p[length(p)] %in% nodelist){
          p<-c(p,nodeparent[[p[length(p)]]])
        }
        nodeparent[[node]]<-p[-1]
        nodelist<-c(nodelist,node)
      }
    }
    return(nodeparent)

  }
  sm<-findMRCA(s,type="node")
  #  write.table(sm, file=paste(repfolder,"/speciestree_mrca.txt",sep=''),sep="\t",quote = F)
  triplet<-matrix(0,nrow=dim(combn(dim(sm)[1],3))[2],ncol = 3)
  gtopo<-matrix(0,nrow = length(g),ncol=dim(combn(dim(sm)[1],3))[2])
  #  tripletlist<-matrix(0,nrow=dim(combn(dim(sm)[1],3))[2],ncol = 3)
  nodes<-matrix(0,nrow = dim(triplet)[1],ncol=2)
  correcttopo<-rep(0,dim(triplet)[1])


  for (i in 1:length(g)){
    g[[i]]$tip.label<-gsub("_0_0$","",g[[i]]$tip.label,perl =T )
    gm<-findMRCA(g[[i]],type="node")
    gm<-gm[row.names(sm),colnames(sm)]
    m=1;
    for(j in 1:(dim(sm)[1]-2)){
      for (k in (j+1):(dim(sm)[1]-1)){
        for (l in (k+1):(dim(sm)[1])){
          x<-c(gm[j,k],gm[j,l],gm[l,k])
          n<-which(x==max(x))
          triplet[m,n]<-triplet[m,n]+1
          gtopo[i,m]<-n
          #	  tripletlist[m,]<-row.names(sm)[c(j,k,l)]
          x<-c(sm[j,k],sm[j,l],sm[l,k])
          nodes[m,]=c(max(x),min(x))
          correcttopo[m]<-which(x==max(x))
          m<-m+1
        }
      }
    }
    #    if (i %% 200==0){
    #      cat(i,"\n")
    #    }
  }
  #  write.table(triplet, file=paste(repfolder,"/triplettopo_freq.txt",sep=''),sep="\t",quote = F,row.names = F,col.names = F)
  #  write.table(gtopo, file=paste(repfolder,"/genetree_triplettopo.txt",sep=''),sep="\t",quote = F,row.names = F,col.names = F)
  #  write.table(tripletlist, file=paste(repfolder,"/tripletlist.txt",sep=''),sep="\t",quote = F,row.names = F,col.names = F)
  #  write.table(correcttopo, file=paste(repfolder,"/triplet_correcttopo.txt",sep=''),sep="\t",quote = F,row.names = F,col.names = F)


  ps<-phy_parents(s)
  toposcore<-matrix(0,nrow = dim(triplet)[1],ncol=3)
  #  topoconsidered<-rep('',dim(triplet)[1])
  for (i in 1:dim(triplet)[1]){
    n_high<-nodes[i,2]
    if(n_high==length(s$tip.label)+1){next}
    p<-ps[[n_high]]
    l<-which(nodes[,1]==nodes[i,1] & nodes[,2] %in% p)
    f<-triplet[i,][-correcttopo[i]]
    f<-mean(f)
    score<-triplet[l,,drop=F]-f
    score[score<0]<-0
    toposcore[l,]=toposcore[l,]+score
    #    topoconsidered[l]<-paste(topoconsidered[l],i,sep=' ')
    #    if(i %% 5000==0){cat(i,"\n")}
  }
  for (i in 1:dim(toposcore)[1]){
    toposcore[i,correcttopo[i]]=0
  }
  #  write.table(topoconsidered, file=paste(repfolder,"/topoconsidered.txt",sep=''),sep="\t",quote = F,row.names = F,col.names = F)
  #  write.table(toposcore, file=paste(repfolder,"/triplet_toposcore.txt",sep=''),sep="\t",quote = F,row.names = F,col.names = F)

  gscore<-apply(gtopo,MARGIN = 1,function(x){
    gs<-0;
    for(i in 1:length(x)){
      gs=gs+toposcore[i,x[i]]
    }
    gs
  })
  #write.table(cbind(dis), file=paste(repfolder,"/genetree_H3scoresum.txt",sep=''),sep="\t",quote = F,row.names = F,col.names = F)
  if(save.middle.file){
    save(gscore,triplet,toposcore, file=paste(repfolder,"/toposcore.RData",sep=''))
  }
  return(data.frame(gscore=gscore,stringsAsFactors = F))
}
