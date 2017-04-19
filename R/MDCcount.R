#this is a function for calculating MDC between gene tree and the species tree
#Not only the MDC counts for the whole tree
#but also MDC counts for each branch
#need to specify the location of the PhyloNet program
#This function return a matrix with two columns, first is the MDC-Zscore for the whole tree, the second is the MDC-bscore based on branches' MDC count, also by default save a MDCcounts.txt file in the folder for each branch.
MDCcount<-function(repfolder=NULL,sfile=NULL, gfile=NULL,phylonet.path,save.middle.file=T){
  require(ape)
  require(phangorn)
  require(phytools)
  require(stringr)
  outfile=paste(repfolder,'/phylonetinput.nexus',sep = '')
  resultfile=paste(repfolder,'/phylonetoutput.txt',sep = '')
  matrixfile=paste(repfolder,'/MDCcounts.txt',sep = '')

  if(is.null(repfolder)){
    s<-read.tree(sfile)
    g<-read.tree(gfile)
    repfolder=sub(basename(gfile),'',gfile,perl = T)
  }else{
    s<-read.tree(paste(repfolder,"/s_tree.trees",sep=''))
    s$edge.length<-NULL
    g<-read.tree(paste(repfolder,"/g_trees.trees",sep=''))
    for (i in 1:length(g)){
      g[[i]]$tip.label<-sub("_0_0$","",g[[i]]$tip.label,perl =T )
      g[[i]]$edge.length<-NULL
      g[[i]]<-drop.tip(g[[i]],tip="outgroup")
    }
  }

  #get all the subtrees of s
  subs<-cut.phy.branch(s)
  resultm<-matrix(0,ncol = length(subs[[1]])+1,nrow = length(g))
  for (gg in 1:length(g)){
    sink(outfile)
    cat("#NEXUS\n\nBEGIN TREES;\n\n",sep = '')
    splist<-write.phylonet.subtrees(s,subs[[1]],'speciesTree')
    glist<-write.phylonet.subtrees(g[[gg]],subs[[1]],'geneTree')
    cat("END;\n\nBEGIN PHYLONET;\n\n",sep = '')
    for (i in 1:length(splist)){
      cat("DeepCoalCount_tree {",splist[i],"} (",glist[i],");\n",sep='')
    }
    cat("END;\n",sep = '')
    sink()
    system(command = paste('java -jar ',phylonet.path,' ',outfile,' >',resultfile,sep=''),wait = T)
    con = file(resultfile, "r")
    while ( TRUE ) {
      line = readLines(con, n = 1)
      if ( length(line) == 0 ) {
        break
      }else if(length(grep("DeepCoalCount",line))>0){
        line=sub(".*geneTree(\\d+).*","\\1",line)
        x=as.integer(line)
      }else if(length(grep("Total number of extra lineages",line))>0){
        resultm[gg,x+1]=as.numeric(sub("\\D+(\\d+.*\\d*)","\\1",line))
      }

    }
    close(con)
    #if(gg %% 100 ==0){
    #  cat(gg,"\n")
    #}
  }
  resultm[,2:dim(resultm)[2]]<-resultm[,1]-resultm[,2:dim(resultm)[2]]

  bresult<-matrix(0,ncol = length(subs[[1]])+1,nrow = length(g))
  bresult[,1]<-resultm[,1]
  for( i in 1:length(subs[[2]])){
    if(is.null(subs[[2]][[i]])){
      bresult[,i+1]<-resultm[,i+1]
    }else{
      bresult[,i+1]=resultm[,i+1]-resultm[,subs[[2]][[i]][1]+1]-resultm[,subs[[2]][[i]][2]+1]
    }
  }
  bresult[bresult<0]<-0
  unlink(outfile)
  unlink(resultfile)
  if(save.middle.file){
    write.table(bresult,file = matrixfile,quote = F,sep = "\t",row.names = F,col.names = F)
  }

  xx<-bresult[,-1]
  Q95<-apply(xx,MARGIN = 2,function(i)quantile(i,0.95))
  Q05<-apply(xx,MARGIN = 2,function(i)quantile(i,0.05))
  IQR<-Q95-Q05
  IQR[IQR==0]<-0.5
  k<-matrix(0,nrow = dim(xx)[1],ncol = dim(xx)[2])
  for(i in 1:dim(xx)[1]){
    l<-xx[i,]-Q95
    l[l<0]=0
    k[i,]<-as.vector(as.numeric(l/IQR))
  }
  MDCbscore<-apply(k+1,MARGIN = 1,prod)
  MDCzscore<-righttailp(bresult[,1])

  MDCresult<-data.frame(MDCzscore=MDCzscore,MDCbscore=MDCbscore,stringsAsFactors = F)
  return(MDCresult)
}
