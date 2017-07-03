source('../../share/sysFacLib.R')

dat1a=clean.s2()
dat1b=clean.s5()
dat2a=clean.s8()
dat2b=clean.s7()


dat=dat1a[dat1a$bothDiff==1,]
y=dat$rt
sub=as.integer(as.factor(dat$sub))
a=dat$sizeChange-1
b=dat$angleChange-1
N=length(a)
I=max(sub)
variances=tapply(y,list(a,b,sub),var)

meanvar=apply(variances,c(1,2),mean)
meanvar
sqrt(meanvar)
