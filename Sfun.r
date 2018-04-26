Sfun<-function(z,c){
  value<-ifelse(abs(z)>c,sign(z)*(abs(z)-c),0)
  return(value)
}