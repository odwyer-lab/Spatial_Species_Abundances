#Numerical method to calculate an SAD in 2d, given a choice of speciation rate and dispersal length scale
sadfreq<-function(n,nu,L,A){
  r<-1
  fac<-exp(-n*log(r))
  cauchy<-function(t){
    return(fac*(exp(-n*t)*(-sad.func2d(r*exp(t),nu,L,A))))
  }
  m<-max(n+1,8)
  tol<-10^(-10)
  im<-complex(real=0,imaginary=1)
  s<-cauchy(2*im*pi*(1:m)/m)
  val1<-mean(s)
  err1<-NaN
  while(m<=2^20){
    #cat(m,val1,"\n")
    m<-2*m
    s<- c(s,cauchy(2*im*pi*seq(1,m,2)/m))
    val<-mean(s)
    kappa<-mean(abs(s))/abs(val)
    err0<-abs(val-val1)/abs(val)
    err<-((err0/err1)^2)*err0
    if(!is.na(err1)){
      if(err <=kappa*tol|!is.finite(kappa)){
        m<-2^20+1
      }
    }
    val1<-val
    err1<-err0
  }
  #if(!is.na(val1)){
  #print(val1)
  if(-log10(abs(Re(val1)))+log10(kappa)<10&Re(val1)>0){
    return(Re(val1))
  }else{
    return(NA)
  }
}