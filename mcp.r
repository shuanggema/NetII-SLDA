mcp<-function(t,lambda,gamma)
{
  value<-lambda*max((1-abs(t)/(lambda*gamma)),0)
  return(value)
}

mcpfun<-function(t,lambda,gamma)
{
  value<-ifelse(abs(t)>gamma*lambda,gamma*lambda^2/2,
                lambda*abs(t)-t^2/(2*gamma))
  return(value)
}