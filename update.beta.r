update.beta<-function(beta0,lambda1,lambda2,gamma)
{
 beta.new<-beta0
 beta1<-beta0[1:p1]
 beta2<-beta0[(p1+1):(p1+p2)]
 beta3<-beta0[(p2+p1+1):length(beta0)]
 beta.sup.new<-rbind(beta1,beta2,beta3)
 k2<-matrix(0,M0,max(p))
 k1<-rep(0,max(p))
  residual<-y-x%*%beta0
 l1norm.beta<-colSums(abs(beta.sup.new))
 
  for(j in 1:max(p))
  {
    k1[j]=lambda2*d/M[j]*(sum(a[j,])-a[j,j])
    for(m in 1:M0)
    {
      
      beta.mj<-beta0[group==j][m]
      lambda.j<-mcp(l1norm.beta[j],sqrt(M[j])*lambda1,gamma)
      value1<-sum(abs(beta0[group==j][-m]))/M[j]
      value2<-l1norm.beta[-j]/(sqrt(M[-j])*sqrt(M[j]))
      part2<-lambda2*d*(a[j,-j]%*%(value1-value2))
      k2[m,j]<-lambda.j+part2?
      k2.mj<-k2[m,j]
      xsup.j<-x[,seq(j,M0*max(p),by=max(p))]
      xmj<-xsup.j[,m]
    
      beta.mj.new<-Sfun(t(xmj)%*%residual/n0+beta.mj*(t(xmj)%*%xmj)/n0,k2.mj)/((t(xmj)%*%xmj)/n0+k1[j])
            residual<-residual-(xmj)%*%(beta.mj.new-beta.mj)
      beta.sup.new[m,j]<-beta.mj.new
    }
    beta.new[group==j]<-beta.sup.new[,j]
      
  }
  return(beta.new)
}
