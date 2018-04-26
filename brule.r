
brule<-function(x,y,x.new)
{
  n1<-sum(y<0)
  n2<-sum(y>0)
  n<-length(y)
  q<-c(n1/n,n2/n)
  x1<-x[which(y<0),]
  x2<-x[which(y>0),]
  mu1<-apply(x1,2,mean)
  mu2<-apply(x2,2,mean)
  s1<-(n1-1)*var(x1)
  s2<-(n2-1)*var(x2)
  s<-(s1+s2)/(n-2)
  f1<-as.vector(log(q[1])-t(mu1)%*%solve(s)%*%mu1/2)+(t(x.new)%*%solve(s))%*%mu1
  f2<-as.vector(log(q[2])-t(mu2)%*%solve(s)%*%mu2/2)+(t(x.new)%*%solve(s))%*%mu2
  value<-f1-f2
  value<-ifelse(value>0,0,1)
  return(list(result=value,f_value=c(f1,f2)))
}
