#numerical implementation of the SAD generating function in O'Dwyer & Cornell
library(Bessel)
sad.func2d <- function(z,nu,L,A){
  m<-(1/L)
  #nu<-nu/(1-nu)
  nueff<-nu/((1 - nu)*(1 - z))*log((1 - (1 - nu)*z)/(nu))
  beta<-1 - z*(1 - nu)
  R<-sqrt(A/pi)
  #print(R)
  phi <-nueff*(1-z)*A+(L)*(1-z)*(1/sqrt(beta))*(1-nueff)*(2*pi*R)*BesselI((1/L)*R*sqrt(beta),1,expon.scaled=TRUE)/(BesselI((1/L)*R*sqrt(beta),0,expon.scaled=TRUE)+BesselI((1/L)*R*sqrt(beta),1,expon.scaled=TRUE)*BesselK((1/L)*R*sqrt(nueff*beta),0,expon.scaled=TRUE)/(BesselK((1/L)*R*sqrt(nueff*beta),1,expon.scaled=TRUE)*sqrt(nueff)))
  return(phi)
}