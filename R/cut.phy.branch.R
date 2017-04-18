# this is a function for generating all the possible trimming
# of a phylogeny (trim one clade at a time)
# Internal function used by the MDCcount function

cut.phy.branch<-function(phy){
  require(ape)
  require(phangorn)
  require(phytools)
  require(stringr)
  phy<-reorder.phylo(phy,order = "post")
  tipnum<-length(phy$tip.label)
  trimtips<-vector("list",2*tipnum-1)
  downnode<-vector("list",2*tipnum-1)
  for( i in 1:dim(phy$edge)[1]){
    if(phy$edge[i,2]<=tipnum){
      trimtips[[phy$edge[i,2]]]<-phy$tip.label[phy$edge[i,2]]
      #downnode[[phy$edge[i,2]]]<-c()
    }
    #if(is.null(trimtips[[phy$edge[i,1]]])){
    #  trimtips[[phy$edge[i,1]]]<-c()
    #  downnode[[phy$edge[i,1]]]<-c()
    #}
    trimtips[[phy$edge[i,1]]]<-c(trimtips[[phy$edge[i,1]]],trimtips[[phy$edge[i,2]]])
    downnode[[phy$edge[i,1]]]<-c(downnode[[phy$edge[i,1]]],phy$edge[i,2])
  }
  return(list(trimtips, downnode))
}
