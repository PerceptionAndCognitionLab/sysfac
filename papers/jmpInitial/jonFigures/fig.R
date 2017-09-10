#### (6/3/2015; JET) Plots for the paper

### ANOVA-like and Surface CI plots (6/3/2015; JET)
## Read in Data
# SysFac02
dat=read.table("~/pcl/active/jonDissertation/analysis/all.dat")
sfdat=read.table("~/pcl/active/jonDissertation/analysis/analyt.dat")
# SysFac05
dat2=read.table("~/pcl/active/jonDissertation/analysis/all2.dat")
sf2dat=read.table("~/pcl/active/jonDissertation/analysis/analyt2append.dat")
# SysFac07
dat3=read.table("~/pcl/active/jonDissertation/analysis/all3.dat")
sf3dat=read.table("~/pcl/active/jonDissertation/analysis/aanalyt3.dat")
# SysFac08
dat4=read.table("~/pcl/active/jonDissertation/analysis/all4.dat")
sf4dat=read.table("~/pcl/active/jonDissertation/analysis/aanalyt4.dat")

names(dat)=names(dat2)=names(dat3)=names(dat4)=names(sfdat)[1:25]

## Surface MIC Plot components
MICcont=c(1/4,-1/4,-1/4,1/4)
gam02=MICcont%*%tapply(sfdat$rt,list(sfdat$chchng,sfdat$nsub),mean)
gam05=MICcont%*%tapply(sf2dat$rt,list(sf2dat$chchng,sf2dat$nsub),mean)
gam07=MICcont%*%tapply(sf3dat$rt,list(sf3dat$chchng,sf3dat$nsub),mean)
gam08=MICcont%*%tapply(sf4dat$rt,list(sf4dat$chchng,sf4dat$nsub),mean)

## Gamma MSEs
sfMSE=rep(NA,max(sfdat$nsub))
sf2MSE=rep(NA,max(sf2dat$nsub))
sf3MSE=rep(NA,max(sf3dat$nsub))
sf4MSE=rep(NA,max(sf4dat$nsub))
for (i in 1:max(sfdat$nsub)){
sfMSE[i]=anova(lm(rt~f1chng*f2chng,sfdat[sfdat$nsub==i,]))[3,3]
}
for (i in 1:max(sf2dat$nsub)){
sf2MSE[i]=anova(lm(rt~f1chng*f2chng,sf2dat[sf2dat$nsub==i,]))[3,3]
}
for (i in 1:max(sf3dat$nsub)){
sf3MSE[i]=anova(lm(rt~f1chng*f2chng,sf3dat[sf3dat$nsub==i,]))[3,3]
}
for (i in 1:max(sf4dat$nsub)){
sf4MSE[i]=anova(lm(rt~f1chng*f2chng,sf4dat[sf4dat$nsub==i,]))[3,3]
}

## Gamma Errors (Misspecified, refer to Hays pp.429-433)
#serr=sqrt(sfMSE/table(sfdat$nsub))
#serr2=sqrt(sf2MSE/table(sf2dat$nsub))
#serr3=sqrt(sf3MSE/table(sf3dat$nsub))
#serr4=sqrt(sf4MSE/table(sf4dat$nsub))
cellN=table(sfdat$chchng,sfdat$nsub)
cell2N=table(sf2dat$chchng,sf2dat$nsub)
cell3N=table(sf3dat$chchng,sf3dat$nsub)
cell4N=table(sf4dat$chchng,sf4dat$nsub)
#MSScale=(MICcont^2)%*%(1/cellN)
#MS2Scale=(MICcont^2)%*%(1/cell2N)
#MS3Scale=(MICcont^2)%*%(1/cell3N)
#MS4Scale=(MICcont^2)%*%(1/cell4N)
#serr=sqrt(sfMSE/MSScale)
#serr2=sqrt(sf2MSE/MS2Scale)
#serr3=sqrt(sf3MSE/MS3Scale)
#serr4=sqrt(sf4MSE/MS4Scale)
serr=sqrt((MICcont^2)%*%(tapply(sfdat$rt,list(sfdat$chchng,sfdat$nsub),var)/cellN))
serr2=sqrt((MICcont^2)%*%(tapply(sf2dat$rt,list(sf2dat$chchng,sf2dat$nsub),var)/cell2N))
serr3=sqrt((MICcont^2)%*%(tapply(sf3dat$rt,list(sf3dat$chchng,sf3dat$nsub),var)/cell3N))
serr4=sqrt((MICcont^2)%*%(tapply(sf4dat$rt,list(sf4dat$chchng,sf4dat$nsub),var)/cell4N))



## Confidence Interval Bounds (Technically Misspecified, but the d.f. is 
## generally large enough to be irrelevant)
cLim=qnorm(.9)
lconf2=as.vector(gam02)-cLim*serr
uconf2=as.vector(gam02)+cLim*serr
lconf5=as.vector(gam05)-cLim*serr2
uconf5=as.vector(gam05)+cLim*serr2
lconf7=as.vector(gam07)-cLim*serr3
uconf7=as.vector(gam07)+cLim*serr3
lconf8=as.vector(gam08)-cLim*serr4
uconf8=as.vector(gam08)+cLim*serr4

## Sort indices for Surface MICs
sort02=rep(NA,length(gam02))
for (i in 1:length(gam02)){sort02[i]=which(gam02==sort(gam02)[i])}

sort05=rep(NA,length(gam05))
for (i in 1:length(gam05)){sort05[i]=which(gam05==sort(gam05)[i])}

sort07=rep(NA,length(gam07))
for (i in 1:length(gam07)){sort07[i]=which(gam07==sort(gam07)[i])}

sort08=rep(NA,length(gam08))
for (i in 1:length(gam08)){sort08[i]=which(gam08==sort(gam08)[i])}

## Ranges for all MIC CIs
range(c(lconf2,uconf2))
range(c(lconf5,uconf5))
range(c(lconf7,uconf7))
range(c(lconf8,uconf8))

## By-Condition Overall Accuracies and (post-cleaning) RTs 
e1Exclusion=as.numeric(setdiff(names(table(dat$sub)),names(table(sfdat$sub))))
e2Exclusion=as.numeric(setdiff(names(table(dat2$sub)),names(table(sf2dat$sub))))
e3Exclusion=as.numeric(setdiff(names(table(dat3$sub)),names(table(sf3dat$sub))))
e4Exclusion=as.numeric(setdiff(names(table(dat4$sub)),names(table(sf4dat$sub))))

## Mean RTs
e1RTCond=(dat$block!=0 & !is.element(dat$sub,e1Exclusion) & dat$ans==1 & dat$rt<5 & dat$rt>.3)
e1RT=as.data.frame(t(tapply(dat$rt[e1RTCond],dat$chchng[e1RTCond],mean)))
e2RTCond=(dat2$block!=0 & !is.element(dat2$sub,e2Exclusion) & dat2$ans==1 & dat2$rt<5 & dat2$rt>.3)
e2RT=as.data.frame(t(tapply(dat2$rt[e2RTCond],dat2$chchng[e2RTCond],mean)))
e3RTCond=(dat3$block!=0 & !is.element(dat3$sub,e3Exclusion) & dat3$ans==1 & dat3$rt<5 & dat3$rt>.3)
e3RT=as.data.frame(t(tapply(dat3$rt[e3RTCond],dat3$chchng[e3RTCond],mean)))
e4RTCond=(dat4$block!=0 & !is.element(dat4$sub,e4Exclusion) & dat4$ans==1 & dat4$rt<5 & dat4$rt>.3)
e4RT=as.data.frame(t(tapply(dat4$rt[e4RTCond],dat4$chchng[e4RTCond],mean)))

## Participant Counts
e1Sub=length(sort02)
e2Sub=length(sort05)
e3Sub=length(sort07)
e4Sub=length(sort08)

## Accuracy Plots
# Plot Parameters
xAlign=c(3,2,1,3,2,1,3,2,1)
pColors=rgb(rep(1:0,c(3,6)),rep(c(0,1,0),each=3),rep(0:1,c(6,3)))
nRing=3

### General Range parameters
range(c(lconf2,lconf5,lconf7,lconf8,uconf2,uconf5,uconf7,uconf8))
# -.3 to .3 for all MIC plots (JET; 6/3/2015)
range(c(e1RT,e2RT,e3RT,e4RT))
# .8 to 2
## Limits and Axis Notation
sicYL=c(-.3,.3)
xCondAx=c("Ø","L","H")
yRTAx=seq(.8,2,by=.2)
rtYL=c(.8,2)


### Experiment 1 Plot
red02=as.numeric(lconf2>0)/2
green02=as.numeric(uconf2<0)/2
blue02=.5-red02-green02

red05=as.numeric(lconf5>0)/2
green05=as.numeric(uconf5<0)/2
blue05=.5-red05-green05

#x11(width=10,height=10)
#pdf("../JE1surf.pdf",width=10,height=10)
par(par.4pan,cex=1.2)

#matplot(matrix(e1RT,nrow=3)[3:1,],col=2:4,typ='o',lty=1,pch=19,
#xlab="Change in Angle",ylab="Reaction Time",axes=F,ylim=rtYL)
matplot(matrix(e1RT,nrow=3)[2:1,1:2],col=c("darkred","darkgreen"),typ='o',
lty=1:2,pch=19,xlab="Change in Angle",ylab="Mean RT (seconds)",axes=F,ylim=rtYL,lwd=2)
#axis(1,at=1:3,labels=xCondAx)
axis(1,at=1:2,labels=xCondAx[2:3])
axis(2,at=yRTAx)
box()
mtext("A.",adj=0,cex=1.5)
#legend(2.1,1.25,c("Ø","L","H"),pch=19,col=4:2,lty=1,
legend(1,1.15,c("Low","High"),pch=19,col=c("darkgreen","darkred"),lty=2:1,lwd=2,
title="Change in Size")
#for (i in 1:nRing){
#points(c(3,2,3,2),e1RT[c(1,2,4,5)],pch=21,cex=1+(i/10))
#}
#matpoints(matrix(e1RT,nrow=3)[3:1,],col=2:4,typ='o',lty=1,pch=19)

plot(sort(gam02),pch=19,ylab="Observed Interaction Contrast (seconds)",col=rgb(red02[sort02],green02[sort02],blue02[sort02],.7),ylim=sicYL)
arrows(1:e1Sub,lconf2[sort02],y1=uconf2[sort02],length=.05,angle=90,code=3,
col=rgb(red02[sort02],green02[sort02],blue02[sort02],.7))
abline(h=0,lwd=2)
mtext("B.",adj=1,cex=1.5)

#matplot(matrix(e2RT,nrow=3)[3:1,],col=2:4,typ='o',lty=1,pch=19,
#xlab="Change in Angle",ylab="Reaction Time",axes=F,ylim=rtYL)
matplot(matrix(e2RT,nrow=3)[2:1,1:2],col=c("darkred","darkgreen"),typ='o',
lty=1,pch=19,xlab="Change in Angle",ylab="Reaction Time",axes=F,ylim=rtYL)
#axis(1,at=1:3,labels=xCondAx)
axis(1,at=1:2,labels=xCondAx[2:3])
axis(2,at=yRTAx)
box()
mtext("C.",adj=0,cex=1.5)
#for (i in 1:nRing){
#points(c(3,2,3,2),e2RT[c(1,2,4,5)],pch=21,cex=1+(i/10))
#}
#matpoints(matrix(e2RT,nrow=3)[3:1,],col=2:4,typ='o',lty=1,pch=19)

plot(sort(gam05),pch=19,ylab="MIC",col=rgb(red05[sort05],green05[sort05],blue05[sort05],.7),ylim=sicYL)
arrows(1:e2Sub,lconf5[sort05],y1=uconf5[sort05],length=.05,angle=90,code=3,
col=rgb(red05[sort05],green05[sort05],blue05[sort05],.7))
abline(h=0,lwd=2)
mtext("D.",adj=1,cex=1.5)

dev.off()



### Experiment 2 Plot
red08=as.numeric(lconf8>0)/2
green08=as.numeric(uconf8<0)/2
blue08=.5-red08-green08

red07=as.numeric(lconf7>0)/2
green07=as.numeric(uconf7<0)/2
blue07=.5-red07-green07

#x11(width=10,height=10)
#pdf("../Exp2SurfPlotV2.pdf",width=10,height=10)
par(par.4pan)

matplot(matrix(e4RT,nrow=3)[2:1,1:2],col=c("darkred","darkgreen"),typ='o',
lty=1,pch=19,xlab="Change in Second Digit",ylab="Reaction Time",axes=F,ylim=rtYL)
axis(1,at=1:2,labels=xCondAx[2:3])
axis(2,at=yRTAx)
box()
mtext("A.",adj=0,cex=1.5)
legend(1,2,c("L","H"),pch=19,col=c("darkgreen","darkred"),lty=1,
title="Change in First Digit")

plot(sort(gam08),pch=19,ylab="MIC",col=rgb(red08[sort08],green08[sort08],blue08[sort08],.7),ylim=sicYL)
arrows(1:e4Sub,lconf8[sort08],y1=uconf8[sort08],length=.05,angle=90,code=3,
col=rgb(red08[sort08],green08[sort08],blue08[sort08],.7))
abline(h=0,lwd=2)
mtext("B.",adj=1,cex=1.5)

matplot(matrix(e3RT,nrow=3)[2:1,1:2],col=c("darkred","darkgreen"),typ='o',
lty=1,pch=19,xlab="Change in Second Digit",ylab="Reaction Time",axes=F,ylim=rtYL)
axis(1,at=1:2,labels=xCondAx[2:3])
axis(2,at=yRTAx)
box()
mtext("C.",adj=0,cex=1.5)

plot(sort(gam07),pch=19,ylab="MIC",col=rgb(red07[sort07],green07[sort07],blue07[sort07],.7),ylim=sicYL)
arrows(1:e3Sub,lconf7[sort07],y1=uconf7[sort07],length=.07,angle=90,code=3,
col=rgb(red07[sort07],green07[sort07],blue07[sort07],.7))
abline(h=0,lwd=2)
mtext("D.",adj=1,cex=1.5)

dev.off()



### Prior plots (6/3/2015; JET)
library(truncnorm)

M=10000
s2=rinvgamma(M,3,.02) #invgamma(3,.02) - mean is .01, sqrt of this is .1
nu=rnorm(M,0,.01) #normal(0,1e-4), sqrt of 1e-4 is .01
gamma1=rnorm(M,mean=nu,sd=sqrt(s2))
gamma2=rnorm(M,mean=nu,sd=sqrt(s2))
dens1=kde2d(gamma1, gamma2,n=100,lims=c(c(-1,1),c(-1,1)))
gamma3=2*gamma1
gamma4=2*gamma2
dens3=kde2d(gamma3, gamma4,n=100,lims=c(c(-1,1),c(-1,1)))

gamma=seq(0,1,.001)
exGamma=c(-.01,0,gamma)
f1=dtruncnorm(gamma,0,Inf,mean(gamma1),sd(gamma1))
f3=dtruncnorm(gamma,0,Inf,mean(gamma3),sd(gamma3))

# Plots (Updated 5/29/2015; JET to reflect actual variation)
#x11(width=10,height=10)
#pdf("../gammaPriorsV2.pdf",width=10,height=10)
par(par.4pan)
contour(dens1,levels=2^seq(-5,4,by=1),
xlab=expression(gamma[i]),ylab=expression(gamma[j]))
abline(v=0)
abline(h=0)
mtext("A.",adj=0,cex=1.3)

plot(exGamma,c(0,0,f1),typ='l',
xlab=expression(gamma[i]),ylab="Density")
mtext("B.",adj=1,cex=1.3)

contour(dens3,levels=2^seq(-5,4,by=1),
xlab=expression(gamma[i]),ylab=expression(gamma[j]))
abline(v=0)
abline(h=0)
mtext("C.",adj=0,cex=1.3)

plot(exGamma,c(0,0,f3),typ='l',
xlab=expression(gamma[i]),ylab="Density")
mtext("D.",adj=1,cex=1.3)

dev.off()




### Relative BF barchart (6/3/2015; JET)
sf02NBF=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/logBF06242015.dat")
sf02GBF=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/logBF02272015b.dat")
sf02SNBF=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/logBF05292015.dat")

sf05NBF=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/logBF03022015a.dat")
sf05GBF=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/logBF03032015b.dat")
sf05SNBF=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/SimplifiedRepository/logBF05292015.dat")

sf07NBF=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/logBF02262015a.dat")
sf07GBF=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/logBF03052015a.dat")
sf07SNBF=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/SimplifiedRepository/logBF05292015.dat")

sf08NBF=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/logBF03132015a.dat")
sf08GBF=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/logBF03172015a.dat")
sf08SNBF=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/SimplifiedRepository/logBF05292015.dat")


sf02NormBFs=c(0,-c(sf02NBF[2,1],sf02NBF[2,2],sf02SNBF[2,2],sf02NBF[2,3],sf02SNBF[2,3]))
sf02GammBFs=c(0,-c(sf02GBF[1,1],sf02GBF[1,2],sf02GBF[1,3]))

sf05NormBFs=c(0,-c(sf05NBF[2,1],sf05NBF[2,2],sf05SNBF[2,2],sf05NBF[2,3],sf05SNBF[2,3]))
sf05GammBFs=c(0,-c(sf05GBF[1,1],sf05GBF[1,2],sf05GBF[1,3]))

sf07NormBFs=c(0,-c(sf07NBF[2,1],sf07NBF[2,2],sf07SNBF[2,2],sf07NBF[2,3],sf07SNBF[2,3]))
sf07GammBFs=c(0,-c(sf07GBF[1,1],sf07GBF[1,2],sf07GBF[1,3]))

sf08NormBFs=c(0,-c(sf08NBF[2,1],sf08NBF[2,2],sf08SNBF[2,2],sf08NBF[2,3],sf08SNBF[2,3]))
sf08GammBFs=c(0,-c(sf08GBF[1,1],sf08GBF[1,2],sf08GBF[1,3]))

range(sf02NormBFs,sf05NormBFs,sf07NormBFs,sf08NormBFs,sf02GammBFs,sf05GammBFs,sf07GammBFs,sf08GammBFs)

Nbarnames=c("Ser","Gen","Par-1","Par-2","Coa-1","Coa-2")
Gbarnames=c("Ser","Gen","Par-1","Coa-1")

is.max=function(x){return(x==max(x))}
is.min=function(x){return(x==min(x))}

x11(width=7,height=5)
par(par.1pan)

barplot(20+sf02NormBFs,names.arg=Nbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf02NormBFs),col='black',xlab="Model",
main="Exp. 1 Perc Norm")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf05NormBFs,names.arg=Nbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf05NormBFs),col='black',xlab="Model",
main="Exp. 1 Mem Norm")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf07NormBFs,names.arg=Nbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf07NormBFs),col='black',xlab="Model",
main="Exp. 2 Mem Norm")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf08NormBFs,names.arg=Nbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf08NormBFs),col='black',xlab="Model",
main="Exp. 2 Perc Norm")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))

barplot(20+sf02GammBFs,names.arg=Gbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf02GammBFs),col='black',xlab="Model",
main="Exp. 1 Perc Gam")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf05GammBFs,names.arg=Gbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf05GammBFs),col='black',xlab="Model",
main="Exp. 1 Mem Gam")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf07GammBFs,names.arg=Gbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf07GammBFs),col='black',xlab="Model",
main="Exp. 2 Mem Gam")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf08GammBFs,names.arg=Gbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf08GammBFs),col='black',xlab="Model",
main="Exp. 2 Perc Gam")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))



x11(width=14,height=10)
#pdf("../E1BFBars.pdf",width=14,height=10)
par(par.4pan)

barplot(20+sf02NormBFs,names.arg=Nbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf02NormBFs),col='black',xlab="Model",
main="Exp. 1 Perc Norm")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf05NormBFs,names.arg=Nbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf05NormBFs),col='black',xlab="Model",
main="Exp. 1 Mem Norm")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf02GammBFs,names.arg=Gbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf02GammBFs),col='black',xlab="Model",
main="Exp. 1 Perc Gam")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf05GammBFs,names.arg=Gbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf05GammBFs),col='black',xlab="Model",
main="Exp. 1 Mem Gam")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))

dev.off()


#pdf("../E2BFBars.pdf",width=14,height=10)
par(par.4pan)

barplot(20+sf08NormBFs,names.arg=Nbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf08NormBFs),col='black',xlab="Model",
main="Exp. 2 Perc Norm")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf07NormBFs,names.arg=Nbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf07NormBFs),col='black',xlab="Model",
main="Exp. 2 Mem Norm")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf08GammBFs,names.arg=Gbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf08GammBFs),col='black',xlab="Model",
main="Exp. 2 Perc Gam")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf07GammBFs,names.arg=Gbarnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf07GammBFs),col='black',xlab="Model",
main="Exp. 2 Mem Gam")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))

dev.off()



### Dot plot(s) (6/4/2015; JET)
norm02Gam=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/gammaN03022015.dat")#sfSample.gamma
gamm02Gam=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/gamma02272015.dat")#gamEst
normS02Gam=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/gammaN05292015.dat")
#gPriorD02=read.table("sf02Analysis/ChainRepository/priorD02272015.dat")
#gPostD02=read.table("sf02Analysis/ChainRepository/postD02272015.dat")
nBrn=55000
gBrn=5000
nThin=50
nM=dim(norm02Gam)[1]
#gM=dim(gammGam)[1]
nKeep=c(-(1:nBrn),-setdiff((nBrn+1):nM,seq(nThin,nM,by=nThin)))
gKeep=-(1:gBrn)

# Surface-Normal Plot
sMICs=gam02#tapply(sfdat$rt,list(sfdat$nsub,sfdat$chchng),mean)%*%c(1,-1,-1,1)
sMIC02=sMICs
nMICBs=matrix(apply(norm02Gam[nKeep,,],2,mean),ncol=4)
nMICSs=matrix(apply(normS02Gam[nKeep,,],2,mean),ncol=4)
normMICs=cbind(nMICBs[,2:4],nMICSs[,3:4])

#x11(width=7,height=7)
#pdf("../Exp1NormEsts.pdf",width=7,height=7)
par(par.1pan)
matplot(t(sMICs),normMICs,pch=19,xlab="Surface MICs (s)",ylab="Normal MICs (s)",
main="Exp. 1 Surface-Normal Plot",col=c(2:4,"darkgreen","darkblue"))
abline(0,1)
legend(.049,-.049,legend=c("General","Parallel-1","Coactive-1","Parallel-2","Coactive-2"),pch=19,col=c(2:4,"darkgreen","darkblue"))
dev.off()

# Surface-Gamma Plot
gMICs=2*matrix(apply(gamm02Gam,2,mean),ncol=4)

#x11(width=7,height=7)
#pdf("../Exp1GammEsts.pdf",width=7,height=7)
#par(par.1pan)
#matplot(t(sMICs),gMICs[,-1],pch=19,xlab="Surface MICs (s)",
#ylab="Gamma MICs (s)",main="Exp. 1 Surface-Gamma Plot",col=2:4)
#abline(0,1)
#legend(-.001,-.05,legend=c("General","Parallel","Coactive"),pch=19,
#col=2:4)
#dev.off()

# Combined Dot Plot(s)
range(sMICs)
range(nMICs)
range(normMICs,gMICs)

#pdf("../Exp1DotPlots.pdf",width=14,height=7)
#x11(width=14,height=7)
axlims=c(-.11,.11)
par(par.2pan)

matplot(t(sMICs),normMICs,pch=19,xlab="Surface MICs (s)",ylab="Normal MICs (s)",
main="Exp. 1 Surface-Normal Plot",col=c(2:4,"darkgreen","darkblue"),
ylim=axlims)
abline(0,1)
legend(.049,-.049,legend=c("General","Parallel-1","Coactive-1","Parallel-2","Coactive-2"),pch=19,col=c(2:4,"darkgreen","darkblue"))
mtext("A.",adj=0,cex=1.3)

matplot(t(sMICs),gMICs[,-1],pch=19,xlab="Surface MICs (s)",ylim=axlims,
ylab="Gamma MICs (s)",main="Exp. 1 Surface-Gamma Plot",col=2:4)
abline(0,1)
mtext("B.",adj=1,cex=1.3)

dev.off()





### Trace plot(s)
subs=33:36
rS=1:4==2
gS=1:4==3
bS=1:4==4

#x11(width=10,height=10)
#pdf("../GeneralTraceV1.pdf",width=10,height=10)
layout(matrix(c(1,2,1,3),nrow=2))
par(par.4pan[-1])
matplot(norm02Gam[-(1:nBrn),subs],typ='l',lty=1,col=rgb(rS,gS,bS,.15),
xlab="Iteration (Post-Burn-In)",ylab="Chain Estimate (s/4)")
mtext("A.",adj=0,cex=1.3)
matplot(norm02Gam[nKeep,subs],typ='l',lty=1,col=rgb(rS,gS,bS,.3),
xlab="Iteration (Post-Burn-In-and-Thin)",ylab="Chain Estimate (s/4)")
mtext("B.",adj=0,cex=1.3)
matplot(gamm02Gam[,subs],typ='l',lty=1,col=rgb(rS,gS,bS,.3),
xlab="Iteration (Post-Burn-In)",ylab="Chain Estimate (s/8)")
mtext("C.",adj=1,cex=1.3)

dev.off()



### Experiment Schematic
screwheads=function(xc,yc,r,theta,inches=F,...){
symbols(xc,yc,circles=r,inches=inches,...)
segments(xc-r*cos(theta),yc-r*sin(theta),xc+r*cos(theta),yc+r*sin(theta),...)
}
x11()
screwheads(c(.25,.75),c(.25,.75),c(.1,.2),c(0,pi/2),inches=F,xlim=c(0,1),ylim=c(0,1),lwd=2)
x11()
screwheads(c(.25,.25,.75,.75),c(.25,.75,.25,.75),c(.1,.125,.15,.2),c(0,pi/4,3*pi/4,pi/2),inches=F,xlim=c(0,1),ylim=c(0,1),lwd=2)
x11()
screwheads(c(.25,.25,.75,.75),c(.25,.75,.25,.75),c(.1,.125,.15,.2),c(0,pi/4,3*pi/4,pi/2),inches=F,xlim=c(0,1),ylim=c(0,1),lwd=2,axes=F,xlab="",ylab="")

symbols(c(.1,.2,.3),c(.1,.2,.3),rectangles=matrix(c(.1,.1,.1,.1,.1,.1),nrow=3),inches=F,lwd=2,xlim=c(0,1),ylim=c(0,1))

#symbols(c(.1,.3,.5),c(.9,.7,.5),rectangles=matrix(c(.15,.15,.15,.2,.2,.2),nrow=3),inches=F,lwd=2,add=T)

symbols(c(.2,.4,.6),c(.9,.65,.4),rectangles=matrix(c(.4,.4,.4,.3,.3,.3),nrow=3),inches=F,lwd=2,add=T,bg="white")

exper.frames=function(xc1,yc1,xw,yh,nframe,xoff=1/2,yoff=5/6,inches=F,bg="white",...){
ord=0:(nframe-1)
xcs=xc1+ord*xw*xoff
ycs=yc1-ord*yh*yoff
rectMat=matrix(rep(c(xw,yh),each=nframe),nrow=nframe)
symbols(xcs,ycs,rectangles=rectMat,inches=inches,bg=bg,...)
}

exper.frames(.2,.8,.4,.3,4,yoff=5/6,ylim=c(-.2,1),xlim=c(0,2),lwd=2,axes=F,xlab="",ylab="")
exper.frames(1.2,.8,.4,.3,3,yoff=5/6,lwd=2,add=T)

exper.frames.span=function(xw,yh,nframe,xoff,yoff){
xextent=nframe*xoff*xw+(1-xoff)*xw
yextent=nframe*yoff*yh+(1-yoff)*yh
return(list(xextent=xextent,yextent=yextent))
}

exper.frames.centers=function(xc1,yc1,xw,yh,nframe,xoff,yoff){
ord=0:(nframe-1)
xcs=xc1+ord*xw*xoff
ycs=yc1-ord*yh*yoff
return(list(xcens=xcs,ycens=ycs))
}

fix.cross=function(xcs,ycs,alen,...){
segments(xcs,ycs-alen,xcs,ycs+alen,...)
segments(xcs-alen,ycs,xcs+alen,ycs,...)
}


#x11(width=10,height=6)
#pdf("../experimentFrames.pdf",width=10,height=6)
exper.frames(.4,.65,.4,.3,2,xlim=c(0,2),ylim=c(0,1.1),lwd=2,axes=F,xlab="",ylab="")
exper.frames(1.2,.9,.4,.3,4,lwd=2,add=T)
fix.cross(c(.4,1.2),c(.65,.9),c(.02,.02),lwd=2)
screwheads(c(.5,.7,1.4,1.8),c(.4,.4,.65,.15),c(.02,.03,.02,.03),c(pi/3,2*pi/3,pi/3,2*pi/3),lwd=2,add=T)
text(c(0,2),c(1.1,1.1),c("A.","B."),cex=1.3)
text(c(.6,.8,.8,1.4,1.6,1.8,2,2)-.05,c(.75,.55,.5,1,.75,.5,.3,.25)+.1,c("2000 ms","Until","Response","500 ms","200 ms","2000 ms","Until","Response"),cex=1.3)

dev.off()



#############################################################################
#### (6/18/2015; JET) New Figures (mainly new dot plots)
### Dot Plots
## Chains
# Exp. 1A
norm02Gam=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/gammaN03022015.dat")
gamm02Gam=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/gamma02272015.dat")
normS02Gam=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/gammaN05292015.dat")
gammS02Gam=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/gamma06152015.dat")

# Exp. 1B
norm05Gam=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/gammaN03022015.dat")
gamm05Gam=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/gamma03032015.dat")
normS05Gam=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/SimplifiedRepository/gammaN05292015.dat")
gammS05Gam=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/SimplifiedRepository/gamma06152015.dat")

# Exp. 2A
norm07Gam=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/gammaN02262015.dat")
gamm07Gam=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/gamma03052015.dat")
normS07Gam=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/SimplifiedRepository/gammaN05292015.dat")
gammS07Gam=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/SimplifiedRepository/gamma06162015.dat")

# Exp. 2B
norm08Gam=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/gammaN03132015.dat")
gamm08Gam=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/gamma03172015.dat")
normS08Gam=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/SimplifiedRepository/gammaN05292015.dat")
gammS08Gam=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/SimplifiedRepository/gamma06162015.dat")



## MIC Estimates
# Exp. 1A
sMIC02=gam02
nMICB02s=matrix(apply(norm02Gam[nKeep,,],2,mean),ncol=4)
nMICS02s=matrix(apply(normS02Gam[nKeep,,],2,mean),ncol=4)
normMIC02s=cbind(nMICB02s[,2:4],nMICS02s[,3:4])
gMICB02s=2*matrix(apply(gamm02Gam,2,mean),ncol=4)
gMICS02s=2*matrix(apply(gammS02Gam,2,mean),ncol=4)
gammMIC02s=cbind(gMICB02s[,2:4],nMICS02s[,3:4])

# Exp. 1B
sMIC05=gam05
nMICB05s=matrix(apply(norm05Gam[nKeep,,],2,mean),ncol=4)
nMICS05s=matrix(apply(normS05Gam[nKeep,,],2,mean),ncol=4)
normMIC05s=cbind(nMICB05s[,2:4],nMICS05s[,3:4])
gMICB05s=2*matrix(apply(gamm05Gam,2,mean),ncol=4)
gMICS05s=2*matrix(apply(gammS05Gam,2,mean),ncol=4)
gammMIC05s=cbind(gMICB05s[,2:4],nMICS05s[,3:4])

# Exp. 2A
sMIC07=gam07
nMICB07s=matrix(apply(norm07Gam[nKeep,,],2,mean),ncol=4)
nMICS07s=matrix(apply(normS07Gam[nKeep,,],2,mean),ncol=4)
normMIC07s=cbind(nMICB07s[,2:4],nMICS07s[,3:4])
gMICB07s=2*matrix(apply(gamm07Gam,2,mean),ncol=4)
gMICS07s=2*matrix(apply(gammS07Gam,2,mean),ncol=4)
gammMIC07s=cbind(gMICB07s[,2:4],nMICS07s[,3:4])

# Exp. 2B
sMIC08=gam08
nMICB08s=matrix(apply(norm08Gam[nKeep,,],2,mean),ncol=4)
nMICS08s=matrix(apply(normS08Gam[nKeep,,],2,mean),ncol=4)
normMIC08s=cbind(nMICB08s[,2:4],nMICS08s[,3:4])
gMICB08s=2*matrix(apply(gamm08Gam,2,mean),ncol=4)
gMICS08s=2*matrix(apply(gammS08Gam,2,mean),ncol=4)
gammMIC08s=cbind(gMICB08s[,2:4],nMICS08s[,3:4])

# Range Checks
range(sMIC02)
range(sMIC05)
range(normMIC02s,normMIC05s,gammMIC02s,gammMIC05s)

range(sMIC07)
range(sMIC08)
range(normMIC07s,normMIC08s,gammMIC07s,gammMIC08s)


#pdf("../Exp1FullDotPlots.pdf",width=14,height=14)
#x11(width=14,height=14)
axlims=c(-.15,.11)
aylims=c(-.15,.11)
par(par.4pan)

matplot(t(sMIC02),normMIC02s,pch=19,xlab="Surface MICs (s)",
ylab="Normal MICs (s)",main="Exp. 1A Surface-Normal Plot",
col=c(2:4,"darkgreen","darkblue"),xlim=axlims,ylim=aylims)
abline(0,1)
legend(.039,-.08,legend=c("General","Parallel-1","Coactive-1","Parallel-2","Coactive-2"),pch=19,col=c(2:4,"darkgreen","darkblue"))
mtext("A.",adj=0,cex=1.3)

matplot(t(sMIC02),gammMIC02s,pch=19,xlab="Surface MICs (s)",xlim=axlims,
ylim=aylims,ylab="Gamma MICs (s)",main="Exp. 1A Surface-Gamma Plot",
col=c(2:4,"darkgreen","darkblue"))
abline(0,1)
mtext("B.",adj=1,cex=1.3)

matplot(t(sMIC05),normMIC05s,pch=19,xlab="Surface MICs (s)",
ylab="Normal MICs (s)",main="Exp. 1B Surface-Normal Plot",
col=c(2:4,"darkgreen","darkblue"),xlim=axlims,ylim=aylims)
abline(0,1)
mtext("C.",adj=0,cex=1.3)

matplot(t(sMIC05),gammMIC05s,pch=19,xlab="Surface MICs (s)",xlim=axlims,
ylim=aylims,ylab="Gamma MICs (s)",main="Exp. 1B Surface-Gamma Plot",
col=c(2:4,"darkgreen","darkblue"))
abline(0,1)
mtext("D.",adj=1,cex=1.3)

dev.off()


#pdf("../Exp2FullDotPlots.pdf",width=14,height=14)
#x11(width=14,height=14)
axlims=c(-.13,.16)
aylims=c(-.13,.16)
par(par.4pan)

matplot(t(sMIC08),normMIC08s,pch=19,xlab="Surface MICs (s)",
ylab="Normal MICs (s)",main="Exp. 2A Surface-Normal Plot",
col=c(2:4,"darkgreen","darkblue"),xlim=axlims,ylim=aylims)
abline(0,1)
legend(.08,-.052,legend=c("General","Parallel-1","Coactive-1","Parallel-2","Coactive-2"),pch=19,col=c(2:4,"darkgreen","darkblue"))
mtext("A.",adj=0,cex=1.3)

matplot(t(sMIC08),gammMIC08s,pch=19,xlab="Surface MICs (s)",xlim=axlims,
ylim=aylims,ylab="Gamma MICs (s)",main="Exp. 2A Surface-Gamma Plot",
col=c(2:4,"darkgreen","darkblue"))
abline(0,1)
mtext("B.",adj=1,cex=1.3)

matplot(t(sMIC07),normMIC07s,pch=19,xlab="Surface MICs (s)",
ylab="Normal MICs (s)",main="Exp. 2B Surface-Normal Plot",
col=c(2:4,"darkgreen","darkblue"),xlim=axlims,ylim=aylims)
abline(0,1)
mtext("C.",adj=0,cex=1.3)

matplot(t(sMIC07),gammMIC07s,pch=19,xlab="Surface MICs (s)",xlim=axlims,
ylim=aylims,ylab="Gamma MICs (s)",main="Exp. 2B Surface-Gamma Plot",
col=c(2:4,"darkgreen","darkblue"))
abline(0,1)
mtext("D.",adj=1,cex=1.3)

dev.off()






### Relative BF barcharts
I02=32
I05=30
I07=I08=28
sf02NNum=matrix(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/lSDNumerN03022015.dat")),nrow=nM,ncol=3)
sf02NDen=matrix(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/lSDDenomN03022015.dat")),nrow=nM,ncol=3)
sf02NCor=matrix(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/lSDCorrN03022015.dat")),nrow=nM,ncol=3)
sf02NNumVec=array(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/lSDNumerVecN03022015.dat")),dim=c(nM,I02,3))
sf02NDenVec=array(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/lSDDenomVecN03022015.dat")),dim=c(nM,I02,3))
sf02NCorVec=array(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/lSDCorrVecN03022015.dat")),dim=c(nM,I02,3))

sf02SNNum=matrix(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/lSDNumerN05292015.dat")),nrow=nM,ncol=3)
sf02SNDen=matrix(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/lSDDenomN05292015.dat")),nrow=nM,ncol=3)
sf02SNCor=matrix(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/lSDCorrN05292015.dat")),nrow=nM,ncol=3)
sf02SNNumVec=array(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/lSDNumerVecN05292015.dat")),dim=c(nM,I02,3))
sf02SNDenVec=array(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/lSDDenomVecN05292015.dat")),dim=c(nM,I02,3))
sf02SNCorVec=array(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/lSDCorrVecN05292015.dat")),dim=c(nM,I02,3))

apply(sf02NNum[nKeep,],2,logMeanExpLogs)-apply(sf02NDen[nKeep,]+sf02NCor[nKeep,],2,logMeanExpLogs)
apply(sf02SNNum[nKeep,],2,logMeanExpLogs)-apply(sf02SNDen[nKeep,]+sf02SNCor[nKeep,],2,logMeanExpLogs)


x11(width=9,height=6)
par(par.6pan)
matplot(sf02NNum,typ='l',ylim=c(-50,200))
abline(h=apply(sf02NNum[nKeep,],2,logMeanExpLogs),col=4)
matplot(sf02NDen,typ='l',ylim=c(-50,200))
abline(h=apply(sf02NDen[nKeep,],2,logMeanExpLogs),col=4)
matplot(sf02NCor,typ='l',ylim=c(-50,200))
abline(h=apply(sf02NCor[nKeep,],2,logMeanExpLogs),col=4)
matplot(sf02SNNum,typ='l',ylim=c(-50,200))
abline(h=apply(sf02SNNum[nKeep,],2,logMeanExpLogs),col=4)
matplot(sf02SNDen,typ='l',ylim=c(-50,200))
abline(h=apply(sf02SNDen[nKeep,],2,logMeanExpLogs),col=4)
matplot(sf02SNCor,typ='l',ylim=c(-50,200))
abline(h=apply(sf02SNCor[nKeep,],2,logMeanExpLogs),col=4)


par.8pan=par.6pan
par.8pan$mfrow=c(2,4)
x11(width=12,height=6)
par(par.8pan)
matplot(sf02NNum,typ='l',ylim=c(-50,200))
abline(h=apply(sf02NNum[nKeep,],2,logMeanExpLogs),col=4)
matplot(sf02NDen,typ='l',ylim=c(-50,200))
abline(h=apply(sf02NDen[nKeep,],2,logMeanExpLogs),col=4)
matplot(sf02NCor,typ='l',ylim=c(-50,200))
abline(h=apply(sf02NCor[nKeep,],2,logMeanExpLogs),col=4)
matplot(sf02NDen+sf02NCor,typ='l',ylim=c(-50,200))
abline(h=apply((sf02NDen+sf02NCor)[nKeep,],2,logMeanExpLogs),col=4)
matplot(sf02SNNum,typ='l',ylim=c(-50,200))
abline(h=apply(sf02SNNum[nKeep,],2,logMeanExpLogs),col=4)
matplot(sf02SNDen,typ='l',ylim=c(-50,200))
abline(h=apply(sf02SNDen[nKeep,],2,logMeanExpLogs),col=4)
matplot(sf02SNCor,typ='l',ylim=c(-50,200))
abline(h=apply(sf02SNCor[nKeep,],2,logMeanExpLogs),col=4)
matplot(sf02SNDen+sf02SNCor,typ='l',ylim=c(-50,200))
abline(h=apply((sf02SNDen+sf02SNCor)[nKeep,],2,logMeanExpLogs),col=4)

apply(sf02NNum[nKeep,],2,logMeanExpLogs)
apply(sf02NDen[nKeep,],2,logMeanExpLogs)
apply(sf02NCor[nKeep,],2,logMeanExpLogs)
apply((sf02NDen+sf02NCor)[nKeep,],2,logMeanExpLogs)

apply(sf02SNNum[nKeep,],2,logMeanExpLogs)
apply(sf02SNDen[nKeep,],2,logMeanExpLogs)
apply(sf02SNCor[nKeep,],2,logMeanExpLogs)
apply((sf02SNDen+sf02SNCor)[nKeep,],2,logMeanExpLogs)


## Check of correction factor to see if it is part of the error
## (6/18/2015; JET)
sf08NCor=matrix(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/lSDCorrN03132015.dat")),nrow=nM,ncol=3)
sf08SNCor=matrix(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/SimplifiedRepository/lSDCorrN05292015.dat")),nrow=nM,ncol=3)
plot(sf08NCor,sf08SNCor)

sf08GCor=matrix(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/priorC03172015.dat")),nrow=10000,ncol=4)
sf08SGCor=matrix(unlist(read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/SimplifiedRepository/priorC06162015.dat")),nrow=100,ncol=4)
#plot(sf08GCor,sf08SGCor)

## Result: Correction factor was not needed. The error was in the M-H sampler
## (6/24/2015; JET). Conclusion: recalculate here.

I02=32
I05=30
I07=I08=28

sf02NNum=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/lSDNumerN06242015.dat")[nKeep,]
sf02NDen=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/lSDDenomN06242015.dat")[nKeep,]

sf02SNNum=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/lSDNumerN05292015.dat")[nKeep,]
sf02SNDen=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/lSDDenomN05292015.dat")[nKeep,]

sf02GNum=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/postD02272015.dat")[-1,]
sf02GDen=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/priorD02272015.dat")[-1,]

sf02SGNum=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/postD06152015.dat")
sf02SGDen=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/priorD06152015.dat")

sf02GNum=apply(array(unlist(sf02GNum),dim=c(5000,I02,4)),c(1,3),prod)
sf02SGNum=apply(array(unlist(sf02SGNum),dim=c(100,I02,4)),c(1,3),prod)

sf02NBF=c(log10(exp(apply(sf02NNum,2,logMeanExpLogs)-apply(sf02NDen,2,logMeanExpLogs))),log10(exp(apply(sf02SNNum,2,logMeanExpLogs)-apply(sf02SNDen,2,logMeanExpLogs)))[2:3])
sf02GBF=c(log10(apply(sf02GNum,2,mean)/apply(sf02GDen,2,mean))[2:4],log10(apply(sf02SGNum,2,mean)/apply(sf02SGDen,2,mean))[3:4])


sf05NNum=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/lSDNumerN06242015.dat")[nKeep,]
sf05NDen=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/lSDDenomN06242015.dat")[nKeep,]

sf05SNNum=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/SimplifiedRepository/lSDNumerN05292015.dat")[nKeep,]
sf05SNDen=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/SimplifiedRepository/lSDDenomN05292015.dat")[nKeep,]

sf05GNum=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/postD03032015.dat")[-(1:5000),]
sf05GDen=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/priorD03032015.dat")[-(1:5000),]

sf05SGNum=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/SimplifiedRepository/postD06152015.dat")
sf05SGDen=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/SimplifiedRepository/priorD06152015.dat")

sf05GNum=apply(array(unlist(sf05GNum),dim=c(5000,I05,4)),c(1,3),prod)
sf05SGNum=apply(array(unlist(sf05SGNum),dim=c(100,I05,4)),c(1,3),prod)

sf05NBF=c(log10(exp(apply(sf05NNum,2,logMeanExpLogs)-apply(sf05NDen,2,logMeanExpLogs))),log10(exp(apply(sf05SNNum,2,logMeanExpLogs)-apply(sf05SNDen,2,logMeanExpLogs)))[2:3])
sf05GBF=c(log10(apply(sf05GNum,2,mean)/apply(sf05GDen,2,mean))[2:4],log10(apply(sf05SGNum,2,mean)/apply(sf05SGDen,2,mean))[3:4])


sf08NNum=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/lSDNumerN06242015.dat")[nKeep,]
sf08NDen=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/lSDDenomN06242015.dat")[nKeep,]

sf08SNNum=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/SimplifiedRepository/lSDNumerN05292015.dat")[nKeep,]
sf08SNDen=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/SimplifiedRepository/lSDDenomN05292015.dat")[nKeep,]

sf08GNum=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/postD03172015.dat")[-(1:5000),]
sf08GDen=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/priorD03172015.dat")[-(1:5000),]

sf08SGNum=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/SimplifiedRepository/postD06162015.dat")
sf08SGDen=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/SimplifiedRepository/priorD06162015.dat")

sf08GNum=apply(array(unlist(sf08GNum),dim=c(5000,I08,4)),c(1,3),prod)
sf08SGNum=apply(array(unlist(sf08SGNum),dim=c(100,I08,4)),c(1,3),prod)

sf08NBF=c(log10(exp(apply(sf08NNum,2,logMeanExpLogs)-apply(sf08NDen,2,logMeanExpLogs))),log10(exp(apply(sf08SNNum,2,logMeanExpLogs)-apply(sf08SNDen,2,logMeanExpLogs)))[2:3])
sf08GBF=c(log10(apply(sf08GNum,2,mean)/apply(sf08GDen,2,mean))[2:4],log10(apply(sf08SGNum,2,mean)/apply(sf08SGDen,2,mean))[3:4])


sf07NNum=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/lSDNumerN06242015.dat")[nKeep,]
sf07NDen=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/lSDDenomN06242015.dat")[nKeep,]

sf07SNNum=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/SimplifiedRepository/lSDNumerN05292015.dat")[nKeep,]
sf07SNDen=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/SimplifiedRepository/lSDDenomN05292015.dat")[nKeep,]

sf07GNum=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/postD03052015.dat")[-(1:5000),]
sf07GDen=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/priorD03052015.dat")[-(1:5000),]

sf07SGNum=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/SimplifiedRepository/postD06162015.dat")
sf07SGDen=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/SimplifiedRepository/priorD06162015.dat")

sf07GNum=apply(array(unlist(sf07GNum),dim=c(5000,I07,4)),c(1,3),prod)
sf07SGNum=apply(array(unlist(sf07SGNum),dim=c(100,I07,4)),c(1,3),prod)

sf07NBF=c(log10(exp(apply(sf07NNum,2,logMeanExpLogs)-apply(sf07NDen,2,logMeanExpLogs))),log10(exp(apply(sf07SNNum,2,logMeanExpLogs)-apply(sf07SNDen,2,logMeanExpLogs)))[2:3])
sf07GBF=c(log10(apply(sf07GNum,2,mean)/apply(sf07GDen,2,mean))[2:4],log10(apply(sf07SGNum,2,mean)/apply(sf07SGDen,2,mean))[3:4])

range(sf02NBF,sf05NBF,sf07NBF,sf08NBF,sf02GBF,sf05GBF,sf07GBF,sf08GBF)

sf02NormBFs=c(0,-sf02NBF)
sf02GammBFs=c(0,-sf02GBF)
sf05NormBFs=c(0,-sf05NBF)
sf05GammBFs=c(0,-sf05GBF)
sf08NormBFs=c(0,-sf08NBF)
sf08GammBFs=c(0,-sf08GBF)
sf07NormBFs=c(0,-sf07NBF)
sf07GammBFs=c(0,-sf07GBF)

## Save these outputs! (6/26/2015; JET)
tabNames=c("0G","0P1","0C1","0P2","0C2")
names(sf02NBF)=tabNames
names(sf05NBF)=tabNames
names(sf08NBF)=tabNames
names(sf07NBF)=tabNames
names(sf02GBF)=tabNames
names(sf05GBF)=tabNames
names(sf08GBF)=tabNames
names(sf07GBF)=tabNames

write.table(sf02NBF,"Exp1ANormBF.dat")
write.table(sf05NBF,"Exp1BNormBF.dat")
write.table(sf08NBF,"Exp2ANormBF.dat")
write.table(sf07NBF,"Exp2BNormBF.dat")
write.table(sf02GBF,"Exp1AGammBF.dat")
write.table(sf05GBF,"Exp1BGammBF.dat")
write.table(sf08GBF,"Exp2AGammBF.dat")
write.table(sf07GBF,"Exp2BGammBF.dat")





range(sf02NormBFs,sf05NormBFs,sf07NormBFs,sf08NormBFs,sf02GammBFs,sf05GammBFs,sf07GammBFs,sf08GammBFs)

Barnames=c("Ser","Gen","Par-1","Coa-1","Par-2","Coa-2")

is.max=function(x){return(x==max(x))}
is.min=function(x){return(x==min(x))}

x11(width=14,height=10)
#pdf("../E1BFBars.pdf",width=14,height=10)
par(par.4pan)

barplot(20+sf02NormBFs,names.arg=Barnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf02NormBFs),col='black',xlab="Model",
main="Exp. 1 Perc Norm")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf05NormBFs,names.arg=Barnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf05NormBFs),col='black',xlab="Model",
main="Exp. 1 Mem Norm")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf02GammBFs,names.arg=Barnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf02GammBFs),col='black',xlab="Model",
main="Exp. 1 Perc Gam")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf05GammBFs,names.arg=Barnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf05GammBFs),col='black',xlab="Model",
main="Exp. 1 Mem Gam")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))

dev.off()


x11(width=14,height=10)
##pdf("../E2BFBars.pdf",width=14,height=10)
par(par.4pan)

barplot(20+sf08NormBFs,names.arg=Barnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf08NormBFs),col='black',xlab="Model",
main="Exp. 2 Perc Norm")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf07NormBFs,names.arg=Barnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf07NormBFs),col='black',xlab="Model",
main="Exp. 2 Mem Norm")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf08GammBFs,names.arg=Barnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf08GammBFs),col='black',xlab="Model",
main="Exp. 2 Perc Gam")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))
barplot(20+sf07GammBFs,names.arg=Barnames,ylab="log BF relative to Serial",
axes=F,density=10*is.max(20+sf07GammBFs),col='black',xlab="Model",
main="Exp. 2 Mem Gam")
axis(2,at=seq(0,20,by=5),labels=seq(-20,0,by=5))

dev.off()



###########################################################################
#### (6/25/2015; JET) New Figures (new dot plots and new trace plot) run in
#### light of the updated Normal analyses
### Dot Plots
## Chains
# Exp. 1A
norm02Gam=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/gammaN06242015.dat")
gamm02Gam=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/gamma02272015.dat")
normS02Gam=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/gammaN05292015.dat")
gammS02Gam=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/SimplifiedRepository/gamma06152015.dat")

# Exp. 1B
norm05Gam=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/gammaN06242015.dat")
gamm05Gam=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/gamma03032015.dat")
normS05Gam=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/SimplifiedRepository/gammaN05292015.dat")
gammS05Gam=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/SimplifiedRepository/gamma06152015.dat")

# Exp. 2B
norm07Gam=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/gammaN06242015.dat")
gamm07Gam=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/gamma03052015.dat")
normS07Gam=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/SimplifiedRepository/gammaN05292015.dat")
gammS07Gam=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/SimplifiedRepository/gamma06162015.dat")

# Exp. 2A
norm08Gam=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/gammaN06242015.dat")
gamm08Gam=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/gamma03172015.dat")
normS08Gam=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/SimplifiedRepository/gammaN05292015.dat")
gammS08Gam=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/SimplifiedRepository/gamma06162015.dat")



## MIC Estimates
# Exp. 1A
sMIC02=gam02
nMICB02s=matrix(apply(norm02Gam[nKeep,,],2,mean),ncol=4)
nMICS02s=matrix(apply(normS02Gam[nKeep,,],2,mean),ncol=4)
normMIC02s=cbind(nMICB02s[,2:4],nMICS02s[,3:4])
gMICB02s=2*matrix(apply(gamm02Gam,2,mean),ncol=4)
gMICS02s=2*matrix(apply(gammS02Gam,2,mean),ncol=4)
gammMIC02s=cbind(gMICB02s[,2:4],nMICS02s[,3:4])

# Exp. 1B
sMIC05=gam05
nMICB05s=matrix(apply(norm05Gam[nKeep,,],2,mean),ncol=4)
nMICS05s=matrix(apply(normS05Gam[nKeep,,],2,mean),ncol=4)
normMIC05s=cbind(nMICB05s[,2:4],nMICS05s[,3:4])
gMICB05s=2*matrix(apply(gamm05Gam,2,mean),ncol=4)
gMICS05s=2*matrix(apply(gammS05Gam,2,mean),ncol=4)
gammMIC05s=cbind(gMICB05s[,2:4],nMICS05s[,3:4])

# Exp. 2B
sMIC07=gam07
nMICB07s=matrix(apply(norm07Gam[nKeep,,],2,mean),ncol=4)
nMICS07s=matrix(apply(normS07Gam[nKeep,,],2,mean),ncol=4)
normMIC07s=cbind(nMICB07s[,2:4],nMICS07s[,3:4])
gMICB07s=2*matrix(apply(gamm07Gam,2,mean),ncol=4)
gMICS07s=2*matrix(apply(gammS07Gam,2,mean),ncol=4)
gammMIC07s=cbind(gMICB07s[,2:4],nMICS07s[,3:4])

# Exp. 2A
sMIC08=gam08
nMICB08s=matrix(apply(norm08Gam[nKeep,,],2,mean),ncol=4)
nMICS08s=matrix(apply(normS08Gam[nKeep,,],2,mean),ncol=4)
normMIC08s=cbind(nMICB08s[,2:4],nMICS08s[,3:4])
gMICB08s=2*matrix(apply(gamm08Gam,2,mean),ncol=4)
gMICS08s=2*matrix(apply(gammS08Gam,2,mean),ncol=4)
gammMIC08s=cbind(gMICB08s[,2:4],nMICS08s[,3:4])

# Range Checks
range(sMIC02)
range(sMIC05)
range(normMIC02s,normMIC05s,gammMIC02s,gammMIC05s)

range(sMIC07)
range(sMIC08)
range(normMIC07s,normMIC08s,gammMIC07s,gammMIC08s)


#pdf("../Exp1FullDotPlots.pdf",width=14,height=14)
#x11(width=14,height=14)
axlims=c(-.15,.11)
aylims=c(-.15,.11)
par(par.4pan)

matplot(t(sMIC02),normMIC02s,pch=19,xlab="Surface MICs (s)",
ylab="Normal MICs (s)",main="Exp. 1A Surface-Normal Plot",
col=c(2:4,"darkgreen","darkblue"),xlim=axlims,ylim=aylims)
abline(0,1)
legend(.039,-.08,legend=c("General","Parallel-1","Coactive-1","Parallel-2","Coactive-2"),pch=19,col=c(2:4,"darkgreen","darkblue"))
mtext("A.",adj=0,cex=1.3)

matplot(t(sMIC02),gammMIC02s,pch=19,xlab="Surface MICs (s)",xlim=axlims,
ylim=aylims,ylab="Gamma MICs (s)",main="Exp. 1A Surface-Gamma Plot",
col=c(2:4,"darkgreen","darkblue"))
abline(0,1)
mtext("B.",adj=1,cex=1.3)

matplot(t(sMIC05),normMIC05s,pch=19,xlab="Surface MICs (s)",
ylab="Normal MICs (s)",main="Exp. 1B Surface-Normal Plot",
col=c(2:4,"darkgreen","darkblue"),xlim=axlims,ylim=aylims)
abline(0,1)
mtext("C.",adj=0,cex=1.3)

matplot(t(sMIC05),gammMIC05s,pch=19,xlab="Surface MICs (s)",xlim=axlims,
ylim=aylims,ylab="Gamma MICs (s)",main="Exp. 1B Surface-Gamma Plot",
col=c(2:4,"darkgreen","darkblue"))
abline(0,1)
mtext("D.",adj=1,cex=1.3)

dev.off()


#pdf("../Exp2FullDotPlots.pdf",width=14,height=14)
#x11(width=14,height=14)
axlims=c(-.13,.16)
aylims=c(-.13,.16)
par(par.4pan)

matplot(t(sMIC08),normMIC08s,pch=19,xlab="Surface MICs (s)",
ylab="Normal MICs (s)",main="Exp. 2A Surface-Normal Plot",
col=c(2:4,"darkgreen","darkblue"),xlim=axlims,ylim=aylims)
abline(0,1)
legend(.08,-.052,legend=c("General","Parallel-1","Coactive-1","Parallel-2","Coactive-2"),pch=19,col=c(2:4,"darkgreen","darkblue"))
mtext("A.",adj=0,cex=1.3)

matplot(t(sMIC08),gammMIC08s,pch=19,xlab="Surface MICs (s)",xlim=axlims,
ylim=aylims,ylab="Gamma MICs (s)",main="Exp. 2A Surface-Gamma Plot",
col=c(2:4,"darkgreen","darkblue"))
abline(0,1)
mtext("B.",adj=1,cex=1.3)

matplot(t(sMIC07),normMIC07s,pch=19,xlab="Surface MICs (s)",
ylab="Normal MICs (s)",main="Exp. 2B Surface-Normal Plot",
col=c(2:4,"darkgreen","darkblue"),xlim=axlims,ylim=aylims)
abline(0,1)
mtext("C.",adj=0,cex=1.3)

matplot(t(sMIC07),gammMIC07s,pch=19,xlab="Surface MICs (s)",xlim=axlims,
ylim=aylims,ylab="Gamma MICs (s)",main="Exp. 2B Surface-Gamma Plot",
col=c(2:4,"darkgreen","darkblue"))
abline(0,1)
mtext("D.",adj=1,cex=1.3)

dev.off()



### Trace plot(s)
subs=33:36
rS=1:4==2
gS=1:4==3
bS=1:4==4

#x11(width=10,height=10)
#pdf("../GeneralTraceV1.pdf",width=10,height=10)
layout(matrix(c(1,2,1,3),nrow=2))
par(par.4pan[-1])
matplot(norm02Gam[-(1:nBrn),subs],typ='l',lty=1,col=rgb(rS,gS,bS,.15),
xlab="Iteration (Post-Burn-In)",ylab="Chain Estimate (s/4)")
mtext("A.",adj=0,cex=1.3)
matplot(norm02Gam[nKeep,subs],typ='l',lty=1,col=rgb(rS,gS,bS,.3),
xlab="Iteration (Post-Burn-In-and-Thin)",ylab="Chain Estimate (s/4)")
mtext("B.",adj=0,cex=1.3)
matplot(gamm02Gam[,subs],typ='l',lty=1,col=rgb(rS,gS,bS,.3),
xlab="Iteration (Post-Burn-In)",ylab="Chain Estimate (s/8)")
mtext("C.",adj=1,cex=1.3)

dev.off()


###########################################################################
#### Residual Plots (7/2/2015; JET)

## Read in Chains
# Exp. 1A
norm02Eta=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/etaN06242015.dat")
norm02Alp=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/alphaN06242015.dat")
norm02Bet=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/betaN06242015.dat")
norm02Sig=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/sigma2N06242015.dat")
gamm02Eta=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/mu02272015.dat")
gamm02Alp=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/alpha02272015.dat")
gamm02Bet=read.table("~/pcl/active/jonDissertation/analysis/sf02Analysis/ChainRepository/beta02272015.dat")

# Exp. 1B
norm05Eta=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/etaN06242015.dat")
norm05Alp=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/alphaN06242015.dat")
norm05Bet=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/betaN06242015.dat")
norm05Sig=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/sigma2N06242015.dat")
gamm05Eta=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/mu03032015.dat")
gamm05Alp=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/alpha03032015.dat")
gamm05Bet=read.table("~/pcl/active/jonDissertation/analysis/sf05Analysis/ChainRepository/beta03032015.dat")

# Exp. 2B
norm07Eta=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/etaN06242015.dat")
norm07Alp=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/alphaN06242015.dat")
norm07Bet=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/betaN06242015.dat")
norm07Sig=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/sigma2N06242015.dat")
gamm07Eta=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/mu03052015.dat")
gamm07Alp=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/alpha03052015.dat")
gamm07Bet=read.table("~/pcl/active/jonDissertation/analysis/sf07Analysis/ChainRepository/beta03052015.dat")

# Exp. 2A
norm08Eta=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/etaN06242015.dat")
norm08Alp=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/alphaN06242015.dat")
norm08Bet=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/betaN06242015.dat")
norm08Sig=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/sigma2N06242015.dat")
gamm08Eta=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/mu03172015.dat")
gamm08Alp=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/alpha03172015.dat")
gamm08Bet=read.table("~/pcl/active/jonDissertation/analysis/sf08Analysis/ChainRepository/beta03172015.dat")


## Posterior Means/Parameter Estimates
norm02Ests=cbind(matrix(apply(norm02Eta[nKeep,],2,mean),ncol=4)[,1],matrix(apply(norm02Alp[nKeep,],2,mean),ncol=4)[,1],matrix(apply(norm02Bet[nKeep,],2,mean),ncol=4)[,1])
gamm02Ests=cbind(matrix(apply(gamm02Eta[gKeep,],2,mean),ncol=4)[,1],matrix(apply(gamm02Alp[gKeep,],2,mean),ncol=4)[,1],matrix(apply(gamm02Bet[gKeep,],2,mean),ncol=4)[,1])

norm05Ests=cbind(matrix(apply(norm05Eta[nKeep,],2,mean),ncol=4)[,1],matrix(apply(norm05Alp[nKeep,],2,mean),ncol=4)[,1],matrix(apply(norm05Bet[nKeep,],2,mean),ncol=4)[,1])
gamm05Ests=cbind(matrix(apply(gamm05Eta[gKeep,],2,mean),ncol=4)[,1],matrix(apply(gamm05Alp[gKeep,],2,mean),ncol=4)[,1],matrix(apply(gamm05Bet[gKeep,],2,mean),ncol=4)[,1])

norm07Ests=cbind(matrix(apply(norm07Eta[nKeep,],2,mean),ncol=4)[,1],matrix(apply(norm07Alp[nKeep,],2,mean),ncol=4)[,1],matrix(apply(norm07Bet[nKeep,],2,mean),ncol=4)[,1])
gamm07Ests=cbind(matrix(apply(gamm07Eta[gKeep,],2,mean),ncol=4)[,1],matrix(apply(gamm07Alp[gKeep,],2,mean),ncol=4)[,1],matrix(apply(gamm07Bet[gKeep,],2,mean),ncol=4)[,1])

norm08Ests=cbind(matrix(apply(norm08Eta[nKeep,],2,mean),ncol=4)[,1],matrix(apply(norm08Alp[nKeep,],2,mean),ncol=4)[,1],matrix(apply(norm08Bet[nKeep,],2,mean),ncol=4)[,1])
gamm08Ests=cbind(matrix(apply(gamm08Eta[gKeep,],2,mean),ncol=4)[,1],matrix(apply(gamm08Alp[gKeep,],2,mean),ncol=4)[,1],matrix(apply(gamm08Bet[gKeep,],2,mean),ncol=4)[,1])

norm02SEst=apply(norm02Sig,2,mean)[1]
norm05SEst=apply(norm05Sig,2,mean)[1]
norm07SEst=apply(norm07Sig,2,mean)[1]
norm08SEst=apply(norm08Sig,2,mean)[1]

norm02Means=norm02Ests[sfdat$nsub,1]+((-1)^sfdat$f1chng)*norm02Ests[sfdat$nsub,2]+((-1)^sfdat$f2chng)*norm02Ests[sfdat$nsub,3]
norm02Resid=(sfdat$rt-norm02Means)/sqrt(norm02SEst)

norm05Means=norm05Ests[sf2dat$nsub,1]+((-1)^sf2dat$f1chng)*norm05Ests[sf2dat$nsub,2]+((-1)^sf2dat$f2chng)*norm05Ests[sf2dat$nsub,3]
norm05Resid=(sf2dat$rt-norm05Means)/sqrt(norm05SEst)

norm07Means=norm07Ests[sf3dat$nsub,1]+((-1)^sf3dat$f1chng)*norm07Ests[sf3dat$nsub,2]+((-1)^sf3dat$f2chng)*norm07Ests[sf3dat$nsub,3]
norm07Resid=(sf3dat$rt-norm07Means)/sqrt(norm07SEst)

norm08Means=norm08Ests[sf4dat$nsub,1]+((-1)^sf4dat$f1chng)*norm08Ests[sf4dat$nsub,2]+((-1)^sf4dat$f2chng)*norm08Ests[sf4dat$nsub,3]
norm08Resid=(sf4dat$rt-norm08Means)/sqrt(norm08SEst)

gamm02Scal=gamm02Ests[sfdat$nsub,1]+((-1)^sfdat$f1chng)*gamm02Ests[sfdat$nsub,2]+((-1)^sfdat$f2chng)*gamm02Ests[sfdat$nsub,3]
gamm02Resid=qnorm(pgamma(sfdat$rt,shape=2,scale=gamm02Scal))

gamm05Scal=gamm05Ests[sf2dat$nsub,1]+((-1)^sf2dat$f1chng)*gamm05Ests[sf2dat$nsub,2]+((-1)^sf2dat$f2chng)*gamm05Ests[sf2dat$nsub,3]
gamm05Resid=qnorm(pgamma(sf2dat$rt,shape=2,scale=gamm05Scal))

gamm07Scal=gamm07Ests[sf3dat$nsub,1]+((-1)^sf3dat$f1chng)*gamm07Ests[sf3dat$nsub,2]+((-1)^sf3dat$f2chng)*gamm07Ests[sf3dat$nsub,3]
gamm07Resid=qnorm(pgamma(sf3dat$rt,shape=2,scale=gamm07Scal))

gamm08Scal=gamm08Ests[sf4dat$nsub,1]+((-1)^sf4dat$f1chng)*gamm08Ests[sf4dat$nsub,2]+((-1)^sf4dat$f2chng)*gamm08Ests[sf4dat$nsub,3]
gamm08Resid=qnorm(pgamma(sf4dat$rt,shape=2,scale=gamm08Scal))

#x11(width=10,height=10)
#pdf("../Exp1AResiduals.pdf",width=10,height=10)
par(par.4pan)
plot(norm02Resid,ylab="Residuals",main="")
plot(gamm02Resid,ylab="Residuals",main="")
qqnorm(norm02Resid,main="Normal Model Q-Q Plot")
abline(0,1)
qqnorm(gamm02Resid,main="Gamma Model Q-Q Plot")
abline(0,1)

dev.off()


#x11(width=10,height=10)
#pdf("../Exp1BResiduals.pdf",width=10,height=10)
par(par.4pan)
plot(norm05Resid,ylab="Residuals",main="")
plot(gamm05Resid,ylab="Residuals",main="")
qqnorm(norm05Resid,main="Normal Model Q-Q Plot")
abline(0,1)
qqnorm(gamm05Resid,main="Gamma Model Q-Q Plot")
abline(0,1)

dev.off()


#x11(width=10,height=10)
#pdf("../Exp2BResiduals.pdf",width=10,height=10)
par(par.4pan)
plot(norm07Resid,ylab="Residuals",main="")
plot(gamm07Resid,ylab="Residuals",main="")
qqnorm(norm07Resid,main="Normal Model Q-Q Plot")
abline(0,1)
qqnorm(gamm07Resid,main="Gamma Model Q-Q Plot")
abline(0,1)

dev.off()


#x11(width=10,height=10)
#pdf("../Exp2AResiduals.pdf",width=10,height=10)
par(par.4pan)
plot(norm08Resid,ylab="Residuals",main="")
plot(gamm08Resid,ylab="Residuals",main="")
qqnorm(norm08Resid,main="Normal Model Q-Q Plot")
abline(0,1)
qqnorm(gamm08Resid,main="Gamma Model Q-Q Plot")
abline(0,1)

dev.off()


######################
#### (7/7/2015; JET) Additional Expoloration
x11(width=10,height=10)
par(par.4pan)
plot(pnorm(norm02Resid),pnorm(gamm02Resid),xlim=c(0,1),ylim=c(0,1))
plot(pnorm(norm05Resid),pnorm(gamm05Resid),xlim=c(0,1),ylim=c(0,1))
plot(pnorm(norm08Resid),pnorm(gamm08Resid),xlim=c(0,1),ylim=c(0,1))
plot(pnorm(norm07Resid),pnorm(gamm07Resid),xlim=c(0,1),ylim=c(0,1))










