#Jeff, Copied from Jon

### ANOVA-like and Surface CI plots (6/3/2015; JET)
## Read in Data
# SysFac02dat
dat=read.table("dat/all.dat")
sfdat=read.table("dat/analyt.dat")
# SysFac05
dat2=read.table("dat/all2.dat")
sf2dat=read.table("dat/analyt2append.dat")
# SysFac07
dat3=read.table("dat/all3.dat")
sf3dat=read.table("dat/aanalyt3.dat")
# SysFac08
dat4=read.table("dat/all4.dat")
sf4dat=read.table("dat/aanalyt4.dat")

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

cellN=table(sfdat$chchng,sfdat$nsub)
cell2N=table(sf2dat$chchng,sf2dat$nsub)
cell3N=table(sf3dat$chchng,sf3dat$nsub)
cell4N=table(sf4dat$chchng,sf4dat$nsub)

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

save(e1RT,gam02,lconf2,uconf2,e2RT,gam05,lconf5,uconf5,file='jE1.Rdat')

save(e4RT,gam08,lconf8,uconf8,e3RT,gam07,lconf7,uconf7,file='jE2.Rdat')



###########################################



# Exp. 1A
norm02Gam=read.table("dat/norm02Gam")

nBrn=55000
gBrn=5000
nThin=50
nM=dim(norm02Gam)[1]
nKeep=c(-(1:nBrn),-setdiff((nBrn+1):nM,seq(nThin,nM,by=nThin)))
gKeep=-(1:gBrn)



normS02Gam=read.table("dat/normS02Gam")
sMIC02=gam02
nMICB02s=matrix(apply(norm02Gam[nKeep,,],2,mean),ncol=4)
nMICS02s=matrix(apply(normS02Gam[nKeep,,],2,mean),ncol=4)
normMIC02s=cbind(nMICB02s[,2:4],nMICS02s[,3:4])

# Exp. 1B
norm05Gam=read.table("dat/norm05Gam")
normS05Gam=read.table("dat/normS05Gam")
sMIC05=gam05
nMICB05s=matrix(apply(norm05Gam[nKeep,,],2,mean),ncol=4)
nMICS05s=matrix(apply(normS05Gam[nKeep,,],2,mean),ncol=4)
normMIC05s=cbind(nMICB05s[,2:4],nMICS05s[,3:4])

save(file="jEst1.Rdat",sMIC02,normMIC02s,sMIC05,normMIC05s)

# Exp. 2A
norm08Gam=read.table("dat/norm08Gam")
normS08Gam=read.table("dat/normS08Gam")
sMIC08=gam08
nMICB08s=matrix(apply(norm08Gam[nKeep,,],2,mean),ncol=4)
nMICS08s=matrix(apply(normS08Gam[nKeep,,],2,mean),ncol=4)
normMIC08s=cbind(nMICB08s[,2:4],nMICS08s[,3:4])

# Exp. 2B
norm07Gam=read.table("dat/norm07Gam")
normS07Gam=read.table("dat/normS07Gam")
sMIC07=gam07
nMICB07s=matrix(apply(norm07Gam[nKeep,,],2,mean),ncol=4)
nMICS07s=matrix(apply(normS07Gam[nKeep,,],2,mean),ncol=4)
normMIC07s=cbind(nMICB07s[,2:4],nMICS07s[,3:4])

save(file="jEst2.Rdat",sMIC08,normMIC08s,sMIC07,normMIC07s)


