# function to calculate Z score using median and variance from the left half of the distribution
righttailp<-function(x){
  l<-density(x)

  m<-l$x[which(l$y==max(l$y))]
  y<-x[x<=m]
  sidesd<-sqrt(sum(sapply(y,function(yy){(yy-m)^2}))/length(y))
  z<-(x-m)/sidesd
  return(z)
}
