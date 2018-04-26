
timestart<-Sys.time()

library(MASS)
library(glmnet)
library(penalizedLDA)
library(pamr)

B=100
gamma=6
rou<-0.5
alpha<-3

M0=3
n1<-200
n2<-200
n3<-200
n<-c(n1,n2,n3)
n0<-sum(n)
n1.test<-300
n2.test<-300
n3.test<-300


p1<-600
p2<-600
p3<-600
p=c(p1,p2,p3)
p0<-sum(p)
group<-rep(1:max(p),3)

beta.best<-matrix(0,nrow=p0,ncol=B)
beta.DSDA1<-matrix(0,nrow=p0,ncol=B)
cv.N1=rep(0,B)
cv.N2=rep(0,B)

NSC.error=rep(0,B)
L1PLD.error=rep(0,B)
DSDA.error=rep(0,B)
bayes.error=rep(0,B)

TP.N1=rep(0,B)
TP.N2=rep(0,B)

TP.NSC=rep(0,B)
TP.L1PLD=rep(0,B)
TP.DSDA=rep(0,B)

FP.N1=rep(0,B)
FP.N2=rep(0,B)

FP.NSC=rep(0,B)
FP.L1PLD=rep(0,B)
FP.DSDA=rep(0,B)
#########################################################

lam1<-2^(seq(-5,-1,by=1/5))
lam2<-2^(seq(-5,-1.5,by=1/5))
sigma1<-diag(1,p1)
for(i in 1:p1)
{for(j in i:p1)
{sigma1[i,j]=rou^abs(i-j)
 sigma1[j,i]=sigma1[i,j]
}
}
beta1<-c(ifelse(rbinom(8,1,0.5)>0,1,-1)*runif(8,0.4,0.8),rep(0,(p1-8)))
beta2<-c(0,0,ifelse(rbinom(8,1,0.5)>0,1,-1)*runif(8,0.4,0.8),rep(0,(p2-10)))
beta3<-c(rep(0,4),ifelse(rbinom(8,1,0.5)>0,1,-1)*runif(8,0.4,0.8),rep(0,(p3-12)))

sigma2<-sigma1
sigma3<-sigma1
beta.cv<-array(0,dim<-c(p0,length(lam2),length(lam1)))

cv<-matrix(0,length(lam2),length(lam1))
beta.best0<-rep(0,p0)
error<-0
l1norm.beta<-rep(0,max(p))
l1norm.beta<-rep(0,max(p))
M<-rep(M0,max(p))
d<-max(M)


for(k in 1:B)
{
  data1<-gen.data(n1,beta1,sigma1,p1)
  y1<-data1$y
  x.sub1<-data1$x.sub
  x11<-data1$x1
  x12<-data1$x2
  
  ##############################end#################################
    data2<-gen.data(n2,beta2,sigma2,p2)
  y2<-data2$y
  x.sub2<-data2$x.sub
  x21<-data2$x1
  x22<-data2$x2
  ##############################end#######################################
  
  data3<-gen.data(n3,beta3,sigma3,p3)
  y3<-data3$y
  x.sub3<-data3$x.sub
  x31<-data3$x1
  x32<-data3$x2
  ##############################end#######################################
  x1<-cbind(scale(x.sub1,scale=FALSE),matrix(0,n1,p2+p3))
  x2<-cbind(matrix(0,n2,p1),scale(x.sub2,scale=FALSE),matrix(0,n2,p3))
  x3<-cbind(matrix(0,n3,p1+p2),scale(x.sub3,scale=FALSE))
  x=rbind(x1,x2,x3)
  x.org<-rbind(x.sub1,x.sub2,x.sub3)
  beta=c(beta1,beta2,beta3)
  y<-c(y1,y2,y3)
  
  test.data1<-gen.data(n1.test,beta1,sigma1,p1)
  y1.test<-ifelse(test.data1$y>0,1,0)
  x1.test<-test.data1$x.sub
  
  test.data2<-gen.data(n2.test,beta2,sigma2,p2)
  y2.test<-ifelse(test.data2$y>0,1,0)
  x2.test<-test.data2$x.sub
  
  test.data3<-gen.data(n3.test,beta3,sigma3,p3)
  y3.test<-ifelse(test.data3$y>0,1,0)
  x3.test<-test.data3$x.sub
  
  n.test=c(n1.test,n2.test,n3.test)
  x1.tran<-cbind(x1.test,matrix(0,n1.test,p2+p3))
  x2.tran<-cbind(matrix(0,n2.test,p1),x2.test,matrix(0,n2.test,p3))
  x3.tran<-cbind(matrix(0,n3.test,p1+p2),x3.test)
  x.test=rbind(x1.tran,x2.tran,x3.tran)
  y.test=c(y1.test,y2.test,y3.test)

  
  #######################################
  ####method 1--bayes rule####
  index1<-which(beta1!=0)
  index2<-which(beta2!=0)
  index3<-which(beta3!=0)
  
  fit1<-brule(x.sub1[,index1],y1,t(x1.test[,index1]))
  fit2<-brule(x.sub2[,index2],y2,t(x2.test[,index2]))
  fit3<-brule(x.sub3[,index3],y3,t(x3.test[,index3]))
  res1<-fit1$result
  res2<-fit2$result
  res3<-fit3$result
  pred1<-c(res1,res2,res3)
  TP<-sum(y.test==1&pred1==1)
  TN<-sum(as.numeric(y.test==0&pred1==0))
  FP<-sum(y.test==0&pred1==1)
  FN<-sum(y.test==1&pred1==0)
  accuracy<-(TP+TN)/(length(y.test))
  error.bayes<-1-accuracy
  bayes.error[k]<-error.bayes
  
  a=afunN1(alpha=alpha,x.org)
  
  for(j in 1:length(lam1))
  {
    
    
    lambda1<-lam1[j]
    for(i in 1:length(lam2))
    { 
      s<-0
      lambda2<-lam2[i]
      fit<-cv.glmnet(x,y,family="gaussian",alpha=0.1,standardize=FALSE)
      beta0<-coef(fit)[-1]
      
      repeat
      {
        beta.old<-beta0
        beta0<-round(as.vector(update.beta(beta.old,lambda1,lambda2,gamma)),3)
        s<-s+1
        if(s>1000||max((beta0-beta.old))<0.0001) break
      }
      s
      which(beta0!=0)
      beta0[which(beta0!=0)]
      beta10<-beta0[1:p1]
      beta20<-beta0[(p1+1):(p1+p2)]
      beta30<-beta0[(p2+p1+1):length(beta0)]
      if((colMeans(x12)-colMeans(x11))%*%beta10<0)
      {
        beta10=-beta10
      }
      if((colMeans(x22)-colMeans(x21))%*%beta20<0)
      {
        beta20=-beta20
      }
      if((colMeans(x32)-colMeans(x31))%*%beta30<0)
      {        
        beta30=-beta30
      }
      beta0<-c(beta10,beta20,beta30)
      s1<-((sum(y1<0)-1)*var(x11)+(n1 -sum(y1<0)-1)*var(x12))/(n1 -2)
      s2<-((sum(y2<0)-1)*var(x21)+(n2 -sum(y2<0)-1)*var(x22))/(n2 -2)
      s3<-((sum(y3<0)-1)*var(x31)+(n3 -sum(y3<0)-1)*var(x32))/(n3 -2)
      intercept1<--(colMeans(x11)+colMeans(x12))%*%beta10/2+
        t(beta10)%*%(s1)%*%beta10*{(t(colMeans(x11)-colMeans(x12))%*%beta10)^(-1)}*log(sum(y1>0)/(n1 -sum(y1>0)))
      intercept2<--(colMeans(x21)+colMeans(x22))%*%beta20/2+
        t(beta20)%*%(s2)%*%beta20*{(t(colMeans(x21)-colMeans(x22))%*%beta20)^(-1)}*log(sum(y2>0)/(n2 -sum(y2>0)))
      intercept3<--(colMeans(x31)+colMeans(x32))%*%beta30/2+
        t(beta30)%*%(s3)%*%beta30*{(t(colMeans(x31)-colMeans(x32))%*%beta30)^(-1)}*log(sum(y3>0)/(n3 -sum(y3>0)))
      
      
      error<-pred(x1.test,y1.test,beta0,x2.test,y2.test,x3.test,y3.test)
      cv[i,j]<-(error$error1*n1.test+error$error2*n2.test+error$error3*n3.test)/(n1.test+n2.test+n3.test)
      beta.cv[,i,j]<-beta0
      message(i,'th lambda1 and ',j,'th lambda2')
    }
  }
  ind<-which(cv==min(cv),arr.ind=TRUE)
  I<-ind[1,1]
  J<-ind[1,2]
  beta.best0<-beta.cv[,I,J]
  beta.best[,k]<-beta.best0
  cv.N1[k]<-min(cv)
  TP.N1[k]<-sum(beta.best0!=0&beta!=0)
  FP.N1[k]<-sum(beta.best0!=0&beta==0)
}
###mean
m.error<-c(mean(bayes.error),mean(cv.N1))
m.TP<-c(0,mean(TP.N1))
m.FP<-c(0,mean(FP.N1))


#median
median.error<-c(median(bayes.error),median(cv.N1))
median.TP<-c(0,median(TP.N1))
median.FP<-c(0,median(FP.N1))

##sd
sd.error<-c(sd(bayes.error),sd(cv.N1))
sd.TP<-c(0,sd(TP.N1))
sd.FP<-c(0,sd(FP.N1))

result<-rbind(median.TP,sd.TP,median.FP,sd.FP,median.error,sd.error)
colnames(result)<-c('Bayes','N1')
rownames(result)<-c('error','error.sd','TP','TP.sd','FP','FP.sd')
write.csv(result,file<-'result for AR(0.5).csv')

timeend<-Sys.time()
run.t<-timeend-timestart
run.t
