gen.data<-function(n,beta,sigma,p)
{
  y<-rbinom(n,size=1,prob=0.5)
  y<-sort(y)
  n1<-sum(y==0)
  n2<-n-n1
  y[y==0]<-(-n/n1)
  y[y==1]<-n/n2
  y1<-y[1:n1]
  y2<-y[(n1+1):n]
  mu1<-rep(0,p) 
  mu2<-sigma%*%beta
  x1<-mvrnorm(n1,mu1,sigma1)
  x2<-mvrnorm(n2,mu2,sigma)
  x.sub<-rbind(x1,x2)
  return(list(y=y,y1=y1,y2=y2,x1=x1,x2=x2,x.sub=x.sub))
}