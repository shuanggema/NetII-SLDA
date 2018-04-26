afunN1<-function(alpha,x.org){
  a<-diag(1,ncol(x.org))
  for(j in 1:ncol(x.org))
  {
    for(k in j:ncol(x.org))
    {
      a[j,k]<-abs(cor(x.org[,j],x.org[,k]))^alpha
      a[k,j]<-a[j,k]
    }
  }
 
  return(a)
}