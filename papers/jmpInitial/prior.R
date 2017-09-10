library('MCMCpack')
library('truncnorm')
#sample approach
M=10000
s2=rinvgamma(M,3,.05)
nu=rnorm(M,0,.2)
gamma1=rnorm(M,mean=nu,sd=sqrt(s2))
gamma2=rnorm(M,mean=nu,sd=sqrt(s2))
dens=kde2d(gamma1, gamma2,n=100,lims=c(c(-1,1),c(-1,1)))
contour(dens,levels=seq(.1,5,.1))
abline(v=0)
abline(h=0)

gamma=seq(0,1,.01)
exGamma=c(-.05,0,gamma)
f=dtruncnorm(gamma,0,Inf,0,.2)
plot(exGamma,c(0,0,f),typ='l')
