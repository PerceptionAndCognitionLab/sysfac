source('../../share/sysFacLib.R')
library(retimes)
library(shape)

dat=clean.s8()

dat=dat[dat$bothDiff==1,]
y=dat$rt
sub=as.integer(as.factor(dat$sub))
a=dat$sizeChange-1
b=dat$angleChange-1
N=length(a)
I=max(sub)	


### true vals

set.seed(123)
R=50
M=5000
serial=matrix(ncol=6,nrow=R)
for (r in 1:R)
{
	t.mu=rnorm(I,.900,.2)
	t.sigma=.120
	t.tau=array(dim=c(2,2,I))
	t.tau[1,1,]=rlnorm(I,log(.05),.2)
	t.tau[1,2,]=t.tau[1,1,]+rlnorm(I,log(.12),.3)
	t.tau[2,1,]=t.tau[1,1,]+rlnorm(I,log(.12),.3)
	t.tau[2,2,]=t.tau[1,2,]+t.tau[2,1,]-t.tau[1,1,]
	scale=t.tau[cbind(a,b,sub)]
	dat$rt=rnorm(N,t.mu[sub],t.sigma)+rexp(N,1/scale)
	serial[r,]=c(1,1/doBayes(dat,r.mean=.15,r.indv=.10,my.iter=M)$BF)
	print(r)
}

coactive=matrix(ncol=6,nrow=R)
for (r in 1:R)
{
	t.mu=rnorm(I,.900,.2)
	t.sigma=.120
	t.tau=array(dim=c(2,2,I))
	t.tau[1,1,]=rlnorm(I,log(.05),.2)
	t.tau[1,2,]=t.tau[1,1,]+rlnorm(I,log(.1),.3)
	t.tau[2,1,]=t.tau[1,1,]+rlnorm(I,log(.1),.3)
	t.tau[2,2,]=t.tau[1,2,]+t.tau[2,1,]-t.tau[1,1,]+rlnorm(I,log(.08),.1)
	scale=t.tau[cbind(a,b,sub)]
	dat$rt=rnorm(N,t.mu[sub],t.sigma)+rexp(N,1/scale)
	coactive[r,]=c(1,1/doBayes(dat,r.mean=.15,r.indv=.10,my.iter=M)$BF)
	print(r)
}
win=coactive[,5]
coac=coactive/win


parallel=matrix(ncol=6,nrow=R)
for (r in 1:R)
{
	t.mu=rnorm(I,.900,.2)
	t.sigma=.120
	t.tau=array(dim=c(2,2,I))
	t.tau[1,1,]=rlnorm(I,log(.05),.2)
	t.tau[1,2,]=t.tau[1,1,]+rlnorm(I,log(.1),.3)
	t.tau[2,1,]=t.tau[1,1,]+rlnorm(I,log(.1),.3)
	t.tau[2,2,]=t.tau[1,2,]+t.tau[2,1,]-t.tau[1,1,]-rlnorm(I,log(.08),.1)
	scale=t.tau[cbind(a,b,sub)]
	dat$rt=rnorm(N,t.mu[sub],t.sigma)+rexp(N,1/scale)
	parallel[r,]=c(1,1/doBayes(dat,r.mean=.15,r.indv=.10,my.iter=M)$BF)
	print(r)
}
win=parallel[,3]
pl=parallel/win


unconstrained=matrix(ncol=6,nrow=R)
for (r in 1:R)
{
	t.mu=rnorm(I,.900,.2)
	t.sigma=.120
	t.tau=array(dim=c(2,2,I))
	t.tau[1,1,]=rlnorm(I,log(.05),.2)
	t.tau[1,2,]=t.tau[1,1,]+rlnorm(I,log(.1),.3)
	t.tau[2,1,]=t.tau[1,1,]+rlnorm(I,log(.1),.3)
	t.tau[2,2,]=pmax(.005,t.tau[1,2,]+t.tau[2,1,]-t.tau[1,1,]+rnorm(I,0,.08))
	scale=t.tau[cbind(a,b,sub)]
	dat$rt=rnorm(N,t.mu[sub],t.sigma)+rexp(N,1/scale)
	unconstrained[r,]=c(1,1/doBayes(dat,r.mean=.15,r.indv=.10,my.iter=M)$BF)
	print(r)
}
win=unconstrained[,6]
u=unconstrained/win


pdf('simFig.pdf',width=6,height=11)

layout(matrix(nrow=3,ncol=3,c(0,1,2,0,3,4,0,5,6),byrow=T),widths=c(.1,1,1))
par(mar=c(4,1,1,1),mgp=c(2,1,0),cex=1.0)
t=seq(0,2.5,.01)
f1=dexgauss(t,.9,.12,.1)
f2=dexgauss(t,.9,.12,.2)
f3=dexgauss(t,.9,.12,.3)
matplot(t,cbind(f1,f2,f3),ty='l',col='black',lty=3:1,lwd=2,axes=F,
	xlab="Time (sec)",ylab="")
axis(1)
par(xpd=NA)
legend(1.6,2.5,legend=c("(1,1)","(1,2), (2,1)","(2,2)"),lty=1:3, lwd=2,cex=.8)
par(xpd=F)
mtext(side=3,adj=.5,cex=1.2,"A.")

t=seq(-.15,.15,.0002)
unc=dnorm(t,0,.08)
coa=dlnorm(t,log(.08),.1)
plot(t,coa,typ='l',ylim=c(0,80),axes=F,ylab="",xlab="Interaction (sec)")
lines(t,unc)
Arrows(0,0,0,80)
axis(1)
mtext(side=3,adj=.5,cex=1.2,"B.")


myPanel=function(x,title="")
{
x=ifelse(x<.001,.0001,x)
p=round(mean(apply(x,1,max)==1,na.rm=T),2)
par(mar=c(4,4,1,1))
myLim=c(-4.3,2)
matplot(log10(t(x)),typ='p',pch=19,col=rgb(1,0,0,.2),axes=F,ylim=myLim,xlab="",ylab="")
abline(h=0,col=rgb(.5,.5,.5,.5))
axis(1, at=1:6, labels = FALSE)
axis(2,at=seq(-4,2,2),label=c(".0001",".01","1","100"),las=0,cex=.8)
lab=c("Serial ","Par-1","Par-2","Coact-1","Coact-2","General")
text(1:6, par("usr")[3], labels = lab, srt = 45, adj = c(1.15,1.15), xpd = TRUE, cex=1)
mtext(side=3,adj=.5,cex=1.2,title,line=-1)
#mtext(side=3,adj=0,cex=1.2,line=-1,round(p,2))
return(p)
}

p=1:4
p[1]=myPanel(serial,"Serial")
mtext(side=2,adj=.5,line=3,"Bayes Factors")
p[2]=myPanel(pl,"Parallel")
p[3]=myPanel(coac, "Coactive")
mtext(side=2,adj=.5,line=3,"Bayes Factors")
p[4]=myPanel(u, "Unconstrained")
dev.off()

