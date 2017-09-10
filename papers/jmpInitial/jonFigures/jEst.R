#source('jLoad.R')
load('jEst1.Rdat')
load('jEst2.Rdat')


axlims=c(-.15,.15)
aylims=c(-.15,.15)
mycol=c('black','darkred','black','darkred','darkblue')
mybcol=c('white','darkred','blue')

estPlot=function(sMic,mic)
{
o=order(t(sMic))
par2=mean(mic[,4])
coa2=mean(mic[,5])
matplot(t(sMic)[o],mic[o,1:3],typ='l',lty=1,col=mycol, 
xlab=expression(paste("Observed MIC, ",hat(M)," (seconds)")),
ylab=expression(paste("Interaction Estimate,",hat(gamma)," (seconds)")),xlim=axlims,ylim=aylims)
matpoints(t(sMic),mic[,1:3],pch=21:24,bg=mybcol,col=mycol,cex=1.4)
abline(h=par2,col="darkred")
abline(h=0,lty=2)
abline(h=coa2,col="darkblue")
}

pdf('jEst1.pdf',width=12,height=6)
par(par.2pan)
estPlot(sMIC02,normMIC02s)
par(xpd=T)
legend(.02,-.05,legend=c("General","Parallel-1","Coactive-1","Parallel-2","Coactive-2"),pch=c(21:23,NA,NA),col=mycol,pt.bg=mybcol,lwd=c(NA,NA,NA,2,2),bg="white")
par(xpd=F)
mtext(side=3,"A.",adj=0,cex=1.3)
mtext(side=3,"Exp 1A\nPerception",adj=.95,cex=1.3,line=-1)
estPlot(sMIC05,normMIC05s)
mtext(side=3,"B.",adj=0,cex=1.3)
mtext(side=3,"Exp 1B\nMemory",adj=.95,cex=1.3,line=-1)
dev.off()

pdf('jEst2.pdf',width=12,height=6)
par(par.2pan)
estPlot(sMIC08,normMIC08s)
par(xpd=T)
legend(.02,-.05,legend=c("General","Parallel-1","Coactive-1","Parallel-2","Coactive-2"),pch=c(21:23,NA,NA),col=mycol,pt.bg=mybcol,lwd=c(NA,NA,NA,2,2),bg="white")
par(xpd=F)
mtext(side=3,"A.",adj=0,cex=1.3)
mtext(side=3,"Exp 2A\nPerception",adj=.95,cex=1.3,line=-1)
estPlot(sMIC07,normMIC07s)
mtext(side=3,"B.",adj=0,cex=1.3)
mtext(side=3,"Exp 2B\nMemory",adj=.95,cex=1.3,line=-1)
dev.off()
