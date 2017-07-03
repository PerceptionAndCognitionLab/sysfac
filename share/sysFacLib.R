library('BayesFactor')
library('MCMCpack')

clean.s2=function()
{
inDat=read.table(url('https://raw.githubusercontent.com/PerceptionCognitionLab/data1/master/sysfactorial/SysFac02/SysFac02.all'))

names=c("sub","trial","block","trialBlock","sizeChange","angleChange","changeType","bothDiff",
"size1","angle1","size2","angle2","resp","acc","rt","bad","badCount",
"tS","tF","tT","tB","tP","tR","percMem","feedback")
colnames(inDat)=names

bad1 = inDat$sub %in% c(4,21,33)
bad2 = inDat$rt<.3 | inDat$rt>5
bad3 = inDat$block==0
bad4 = inDat$acc==0
bad=(bad1 | bad2 | bad3 |bad4)
return(inDat[!bad,])	
}

clean.s5=function(hard=F)
{
inDat=read.table(url('https://raw.githubusercontent.com/PerceptionCognitionLab/data1/master/sysfactorial/SysFac05/SysFac05.all'))

names=c("sub","trial","block","trialBlock","sizeChange","angleChange","changeType","bothDiff",
"size1","angle1","size2","angle2","resp","acc","rt","bad","badCount",
"tS","tF","tT","tB","tP","tR","percMem","feedback")
colnames(inDat)=names

easyTask=as.numeric(inDat$sub>21) #as per Jon's notes in database
inDat=inDat[easyTask==1,]

bad1=inDat$sub %in% c(43,24,39,26,33,30,35,49)
bad2 = inDat$rt<.3 | inDat$rt>5
bad3 = inDat$block==0
bad4 = inDat$acc==0
bad=(bad2 | bad3 |bad4)
if (hard) bad=bad|bad1
return(inDat[!bad,])
}

#size is 10s, angle is 1s

clean.s7=function()
{
inDat=read.table(url('https://raw.githubusercontent.com/PerceptionCognitionLab/data1/master/sysfactorial/SysFac07/SysFac07.all'))
names=c("sub","trial","block","trialBlock","sizeChange","angleChange","changeType","bothDiff",
"size1","angle1","size2","angle2","resp","acc","rt","bad","badCount",
"tS","tF","tT","tB","tP","tR","percMem","feedback")
colnames(inDat)=names

bad2 = inDat$rt<.3 | inDat$rt>5
bad3 = inDat$block==0
bad4 = inDat$acc==0
bad=(bad2 | bad3 |bad4)
return(inDat[!bad,])
}

clean.s8=function()
{
inDat=read.table(url('https://raw.githubusercontent.com/PerceptionCognitionLab/data1/master/sysfactorial/SysFac08/SysFac08.all'))
names=c("sub","trial","block","trialBlock","sizeChange","angleChange","changeType","bothDiff",
"size1","angle1","size2","angle2","resp","acc","rt","bad","badCount",
"tS","tF","tT","tB","tP","tR","percMem","feedback")
colnames(inDat)=names

bad1= inDat$sub == 27
bad2 = inDat$rt<.3 | inDat$rt>5
bad3 = inDat$block==0
bad4 = inDat$acc==0
bad=(bad1 | bad2 | bad3 |bad4)
return(inDat[!bad,])
}



doBayes=function(dat,r.mean=.167,r.indv=.111,my.iter=10000)
{
	dat=dat[dat$bothDiff==1,]
	y=dat$rt
	sub=as.integer(as.factor(dat$sub))
	a=dat$sizeChange-1
	b=dat$angleChange-1
	N=length(a)
	I=max(sub)
	###Unstructured
	X=matrix(0,nrow=N,ncol=4*I+3)
	for (n in 1:N)
	{
		X[n,sub[n]]=1   #grand mean
		X[n,(I+sub[n])]=2*(a[n]-3/2)  #size deviate
		X[n,(2*I+sub[n])]=2*(b[n]-3/2) #angle deviate
		X[n,(3*I+sub[n])]=4*(a[n]-3/2)*(b[n]-3/2) #interaction deviate
		X[n,(4*I+1)]=2*(a[n]-3/2)  #size mean
		X[n,(4*I+2)]=2*(b[n]-3/2)  #angle mean
		X[n,(4*I+3)]=4*(a[n]-3/2)*(b[n]-3/2) #interaction mean	
	}
	g=c(rep(0:3,each=I),4,5,6)
	gen=nWayAOV(y,X,g,rscale=c(1,rep(r.indv,3),rep(r.mean,3)),iterations=my.iter)
	genSamp=nWayAOV(y,X,g,rscale=c(1,rep(r.indv,3),rep(r.mean,3)),posterior=T,iterations=my.iter)
	###Serial
	Xs=matrix(0,nrow=N,ncol=3*I+2)
	for (n in 1:N)
	{
		Xs[n,sub[n]]=1   #grand mean
		Xs[n,(I+sub[n])]=2*(a[n]-3/2)  #size deviate
		Xs[n,(2*I+sub[n])]=2*(b[n]-3/2) #angle deviate
		Xs[n,(3*I+1)]=2*(a[n]-3/2)  #size mean
		Xs[n,(3*I+2)]=2*(b[n]-3/2)  #angle mean
	}
	gs=c(rep(0:2,each=I),3,4)
	serial=nWayAOV(y,Xs,gs,rscale=c(1,rep(r.indv,2),rep(r.mean,2)),iterations=my.iter)
	BFs.u=exp(serial$bf-gen$bf)
	# Common Interaction Effect
	X1=matrix(0,nrow=N,ncol=3*I+3)
	for (n in 1:N)
	{
		X1[n,sub[n]]=1   #grand mean
		X1[n,(I+sub[n])]=2*(a[n]-3/2)  #size deviate
		X1[n,(2*I+sub[n])]=2*(b[n]-3/2) #angle deviate
		X1[n,(3*I+1)]=2*(a[n]-3/2)  #size mean
		X1[n,(3*I+2)]=2*(b[n]-3/2)  #angle mean
		X1[n,(3*I+3)]=4*(a[n]-3/2)*(b[n]-3/2) #interaction mean	
	}
	g1=c(rep(0:2,each=I),3:5)
	one=nWayAOV(y,X1,g1,rscale=c(1,rep(r.indv,2),rep(r.mean,3)),iterations=my.iter)
	oneSamp=nWayAOV(y,X1,g1,rscale=c(1,rep(r.indv,2),rep(r.mean,3)),posterior=T,iterations=my.iter)
	oneGamma=oneSamp[,3*I+4]
	posMult=(mean(oneGamma>0))/.5
	negMult=(mean(oneGamma<0))/.5
	BFcov2.o=posMult
	BFpar.o=negMult
	BFs.cov2=exp(serial$bf-one$bf-log(posMult))
	BFs.par2=exp(serial$bf-one$bf-log(negMult))
	#prior computation for par1, cov1
	M=100000
	all=1:M
	mu.gamma=rcauchy(M,0,r.mean)
	s2.gamma=rinvgamma(M,.5,.5*r.indv*r.indv)
	for (m in 1:M)
	{
		gamma=rnorm(I,mu.gamma[m],sqrt(s2.gamma[m]))
		all[m]=mean(gamma>0)==1
	}
	priorVal=mean(all)
	myIntDev=genSamp[,(2+3*I):(4*I+1)]
	myIntMean=genSamp[,(4*I)+4]
	myInt=myIntDev+myIntMean
	pos=myInt>0
	neg=myInt<0
	posPost=mean((apply(pos,1,mean)==1))
	negPost=mean((apply(neg,1,mean)==1))
	BFpar1.u=negPost/priorVal
	BFcov1.u=posPost/priorVal
	BFs.par1=BFs.u/BFpar1.u
	BFs.cov1=BFs.u/BFcov1.u	
	BF=c(BFs.par1,BFs.par2,BFs.cov1,BFs.cov2,BFs.u)
	names(BF)=c('s.par1','s.par2','s.cov1','s.cov2','s.u')
	return(list(BF=BF,mcmc=myInt))
}

doMeanPlotA=function(dat,...)
{
	dat=dat[dat$bothDiff==1,]
	y=dat$rt
	sub=as.integer(as.factor(dat$sub))
	a=dat$sizeChange-1
	b=dat$angleChange-1
	N=length(a)
	I=max(sub)
	mrt=tapply(y,list(a,b),mean)
    se=sqrt(mean(tapply(y,list(a,b,sub),var))/I) 
	matplot(1:2,t(mrt),axes=F,ylab="Response Time (sec)",lty=1:2,lwd=2,col=c('black','darkblue'),typ='l',xlim=c(.9,2.1),xlab="Change in Angle",...)
	matpoints(1:2,t(mrt),pch=21,bg=c('white','blue'),col=1,cex=1.2)
    axis(2)
	axis(1,at=1:2,labels=c("Low","High"))
	box()
	m1=((mrt[1,1]-mrt[1,2])+(mrt[2,1]-mrt[2,2]))/2
	m2=((mrt[1,1]-mrt[2,1])+(mrt[1,2]-mrt[2,2]))/2
	return(c(m1,m2))
}

mySE=function(ao,I)
{
	SS=ao[[2]][[1]][[2]][2]+ao[[3]][[1]][[2]][2]+ao[[4]][[1]][[2]][2]
	df=ao[[2]][[1]][[1]][2]+ao[[3]][[1]][[1]][2]+ao[[4]][[1]][[1]][2]
	se=sqrt((SS/df)/I)
	return(se)
}

doMeanPlotAwErr=function(dat,...)
{
	dat=dat[dat$bothDiff==1,]
	mrt=tapply(dat$rt,list(dat$sizeChange,dat$angleChange,dat$sub),mean)
	newDat=as.data.frame.table(mrt)	
	colnames(newDat)=c("size","angle","sub","rt")
	ao=summary(aov(rt~size*angle+Error(sub/(size*angle)),data=newDat))
	I=length(levels(newDat$sub))
	eb=mySE(ao,I)*qt(.975,I-1)
	mrt=tapply(newDat$rt,list(newDat$size,newDat$angle),mean)
	matplot(1:2,t(mrt),axes=F,ylab="Response Time (sec)",lty=1:2,lwd=2,col=c('black','darkblue'),typ='l',xlim=c(.9,2.1),xlab="Change in Angle",...)
	arrows(1:2,t(mrt)[,1]-eb,1:2,t(mrt)[,1]+eb,code=3,angle=90,length=.1)
	arrows(1:2,t(mrt)[,2]-eb,1:2,t(mrt)[,2]+eb,code=3,angle=90,length=.1)
	matpoints(1:2,t(mrt),pch=21,bg=c('white','blue'),col=1,cex=1.2)
    axis(2)
	axis(1,at=1:2,labels=c("Low","High"))
	box()	
	m1=((mrt[1,1]-mrt[1,2])+(mrt[2,1]-mrt[2,2]))/2
	m2=((mrt[1,1]-mrt[2,1])+(mrt[1,2]-mrt[2,2]))/2
	return(c(m1,m2))
}


doMeanPlotB=function(dat,...)
{
	dat=dat[dat$bothDiff==1,]
	y=dat$rt
	sub=as.integer(as.factor(dat$sub))
	a=dat$sizeChange-1
	b=dat$angleChange-1
	N=length(a)
	I=max(sub)
	mrt=tapply(y,list(a,b),mean)
	matplot(1:2,t(mrt),axes=F,ylab="Response Time (sec)",lty=1:2,lwd=2,col=c('black','darkblue'),typ='l',xlim=c(.9,2.1),xlab="Change in 1s Digit",...)
	matpoints(1:2,t(mrt),pch=21,bg=c('white','blue'),col=1,cex=1.2)
    axis(2)
	axis(1,at=1:2,labels=c("Low","High"))
	box()
	m1=((mrt[1,1]-mrt[1,2])+(mrt[2,1]-mrt[2,2]))/2
	m2=((mrt[1,1]-mrt[2,1])+(mrt[1,2]-mrt[2,2]))/2
	return(c(m1,m2))
}


doMeanPlotBwErr=function(dat,...)
{
	dat=dat[dat$bothDiff==1,]
	mrt=tapply(dat$rt,list(dat$sizeChange,dat$angleChange,dat$sub),mean)
	newDat=as.data.frame.table(mrt)	
	colnames(newDat)=c("size","angle","sub","rt")
	ao=summary(aov(rt~size*angle+Error(sub/(size*angle)),data=newDat))
	I=length(levels(newDat$sub))
	eb=mySE(ao,I)*qt(.975,I-1)
	mrt=tapply(newDat$rt,list(newDat$size,newDat$angle),mean)
	matplot(1:2,t(mrt),axes=F,ylab="Response Time (sec)",lty=1:2,lwd=2,col=c('black','darkblue'),typ='l',xlim=c(.9,2.1),xlab="Change in Digit",...)
	arrows(1:2,t(mrt)[,1]-eb,1:2,t(mrt)[,1]+eb,code=3,angle=90,length=.1)
	arrows(1:2,t(mrt)[,2]-eb,1:2,t(mrt)[,2]+eb,code=3,angle=90,length=.1)
	matpoints(1:2,t(mrt),pch=21,bg=c('white','blue'),col=1,cex=1.2)
    axis(2)
	axis(1,at=1:2,labels=c("Low","High"))
	box()	
	m1=((mrt[1,1]-mrt[1,2])+(mrt[2,1]-mrt[2,2]))/2
	m2=((mrt[1,1]-mrt[2,1])+(mrt[1,2]-mrt[2,2]))/2
	return(c(m1,m2))
}

	
doEmpiricalPlot=function(dat,plot=T,...)
{
	dat=dat[dat$bothDiff==1,]
	y=dat$rt
	sub=as.integer(as.factor(dat$sub))
	a=dat$sizeChange-1
	b=dat$angleChange-1
	N=length(a)
	I=max(sub)	
	mrts=tapply(y,list(sub,a,b),mean)	
	myMIC=function(mat) {return((mat[1,1]+mat[2,2]-mat[1,2]-mat[2,1])/4)}
	mic=apply(mrts,1,myMIC)
	sd=dof=1:I
	for (i in 1:I)
	{
        indvDat=dat[sub==i,]
        out=summary(aov(rt~sizeChange*angleChange,data=indvDat))
        sd[i]=sqrt(out[[1]][3][4,1])
        dof[i]=out[[1]][1][4,1]
	}
	c80=qt(.9,dof)*(sd/sqrt(dof-1))
	o=order(mic)
	if (plot){
		plot(1:I,mic[o],ylab="Observed MIC (sec)",typ='n',xlab="Participants",...)
		abline(h=0)
		arrows(1:I,mic[o]-c80[o],1:I,mic[o]+c80[o],code=3,angle=90,length=.1)
		hi=mic[o]-c80[o]>0
		lo=mic[o]+c80[o]<0
		myCol=rep("white",I)
		myCol[lo]="darkred"
		myCol[hi]="yellow"	
		points(1:I,mic[o],bg=myCol,pch=21,cex=1.2)}
	return(mic)
}

polyCI <- function(upper, lower, col){
  len <- length(upper)
  polygon(x = c(rev(1 : len), 1 : len)
          , y = c(rev(lower), upper)
          , col = col
          , border = NA)}


doEstimatePlot=function(mic,est,inRange=c(0,0))
{
	I=length(mic)
	pm=apply(est,2,mean)
	o=order(mic)
	qL=apply(est,2,quantile,.025)
	qU=apply(est,2,quantile,.975)
	range=ifelse(inRange==c(0,0),c(min(mic,qL),max(mic,qU)),inRange)
	plot(1:I,pm,typ='n',xlab="Participants",ylab="MIC Estimate (sec)",ylim=range)
	polyCI(qU[o],qL[o],'grey90')
	abline(h=0)
	lines(1:I,mic[o],lty=2)
	points(1:I,mic[o],pch=21,bg='white')
	lines(1:I,pm[o])
	points(1:I,pm[o],pch=21,bg='red')
}
	
