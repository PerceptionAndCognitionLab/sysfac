#source('jLoad.R')
load('jE1.Rdat')
load('jE2.Rdat')


meanRTplot=function(eRT,label,leg=T,leg.y=1.5,ylimits=c(1.2,2.0))
{
matplot(matrix(eRT,nrow=3)[2:1,1:2],col=c("darkslateblue","black"),typ='o',
lty=1:2,pch=c(19,1),xlab="Change in Angle",ylab="Mean RT (seconds)",axes=F,ylim=ylimits,lwd=2)
axis(1,at=1:2,labels=c("Low","High"))
axis(2)
box()
mtext(label,adj=0,cex=1.4)
if (leg){
legend(1,leg.y,c("Low     ","High"),pch=c(1,19),col=c("black","darkslateblue"),lty=2:1,lwd=2,
title="Change in Size")}
}

ciPlot=function(mic,lo,hi,label,leg=T,ylimits=c(-.25,.25))
{
mycol=c("white","darkred","blue")
mypch=21
o=order(mic)
mic=mic[o]
lo=lo[o]
hi=hi[o]
sub=1:length(mic)
ciType=rep(1,length(mic))
ciType[lo>0]=3
ciType[hi<0]=2
plot(mic,pch=19,ylab=expression(paste("Observed MIC ",hat(M)," (seconds)")),ylim=ylimits,typ='n',xlab="Participants")
abline(h=0,lty=2)
arrows(sub,lo,sub,hi,angle=90,code=3,length=.05,lwd=2)
points(sub,mic,pch=mypch,bg=mycol[ciType],cex=1.2)
text(sub[ciType==2],lo[ciType==2]-.04,"*",adj=c(.5,.5),cex=1.4)
text(sub[ciType==3],hi[ciType==3]+.04,"*",adj=c(.5,.5),cex=1.4)
if(leg) legend(5,.3,legend=c("Serial","Parallel","Coactive"),pch=mypch,pt.bg=mycol)
mtext(label,adj=1,cex=1.4)
}

pdf("jE1.pdf",width=10,height=10)
par(mfrow=c(2,2),mar=c(4,4,3,2),cex=1.2,mgp=c(2,1,0))
meanRTplot(e1RT,"A.")
mtext("Perception",side=3,adj=1.4,cex=1.4,line=1)
ciPlot(as.vector(gam02),lconf2,uconf2,"B.",leg=F)
meanRTplot(e2RT,"C.",leg=F,ylimits=c(1,1.3))
mtext("Memory",side=3,adj=1.4,cex=1.4,line=1)
ciPlot(as.vector(gam05),lconf5,uconf5,"D.",leg=F)
dev.off()


pdf("jE2.pdf",width=10,height=10)
par(mfrow=c(2,2),mar=c(4,4,3,2),cex=1.2,mgp=c(2,1,0))
meanRTplot(e4RT,"A.",leg.y=1.45,ylimits=c(1,1.5))
mtext("Perception",side=3,adj=1.4,cex=1.4,line=1)
ciPlot(as.vector(gam08),lconf8,uconf8,"B.",leg=F,ylimits=c(-.25,.28))
meanRTplot(e3RT,"C.",leg=F,ylimits=c(.9,1.2))
mtext("Memory",side=3,adj=1.4,cex=1.4,line=1)
ciPlot(as.vector(gam07),lconf7,uconf7,"D.",leg=F,ylimits=c(-.25,.2))
dev.off()

