#### (5/27/2015; JET) Combined prior density plots; derived from "prior.R"
#### with updated prior information from the analyses

## MCMC approach
library('MCMCpack')
library('truncnorm')
# Generate Samples 
M=10000
s2=rinvgamma(M,3,.02) #invgamma(3,.02) - mean is .01, sqrt of this is .1
nu=rnorm(M,0,.01) #normal(0,1e-4), sqrt of 1e-4 is .01
gamma1=rnorm(M,mean=nu,sd=sqrt(s2))
gamma2=rnorm(M,mean=nu,sd=sqrt(s2))
dens=kde2d(gamma1, gamma2,n=100,lims=c(c(-1,1),c(-1,1)))
contour(dens,levels=seq(.1,5,.1))
abline(v=0)
abline(h=0)

gamma=seq(0,1,.01)
exGamma=c(-.05,0,gamma)
f=dtruncnorm(gamma,0,Inf,0,.1)# mean is mean of nu, sd is sqrt of mean of s2
plot(exGamma,c(0,0,f),typ='l')

## Preview of Combined Plot
x11(width=10,height=5)
par(par.2pan)
contour(dens,levels=2^seq(-5,5,by=.5))
abline(v=0)
abline(h=0)
mtext("A.",adj=0,cex=1.3)

plot(exGamma,c(0,0,f),typ='l')
mtext("B.",adj=1,cex=1.3)


## Final Setup
# Load Packages
library('MCMCpack')
library('truncnorm')

# Samples
M=1000000
s2=rinvgamma(M,3,.02) #invgamma(3,.02) - mean is .01, sqrt of this is .1
nu=rnorm(M,0,.01) #normal(0,1e-4), sqrt of 1e-4 is .01
gamma1=rnorm(M,mean=nu,sd=sqrt(s2))
gamma2=rnorm(M,mean=nu,sd=sqrt(s2))
dens=kde2d(gamma1, gamma2,n=100,lims=c(c(-1,1),c(-1,1)))

gamma=seq(0,1,.001)
exGamma=c(-.01,0,gamma)
f=dtruncnorm(gamma,0,Inf,mean(gamma1),sd(gamma1))

# Plots (Updated 5/29/2015; JET to reflect actual variation)
#x11(width=10,height=5)
pdf("gammaPriors.pdf",width=10,height=5)
par(par.2pan)

plot(exGamma,c(0,0,f),typ='l',xlim=c(0,.5),
xlab=expression(paste("Parameter ", gamma[0]," (seconds)")),ylab="Density")
mtext("A.",adj=0,cex=1.3)


contour(dens,#levels=seq(.2,8,by=.2),
levels=2^seq(-3,4,by=1),xlim=c(-.5,.5),ylim=c(-.5,.5),
xlab=expression(paste("Parameter ",gamma[i]," (seconds)")),ylab=expression(paste("Parameter ",gamma[j]," (seconds)")))
abline(v=0)
abline(h=0)
mtext("B.",adj=1,cex=1.3)

dev.off()


