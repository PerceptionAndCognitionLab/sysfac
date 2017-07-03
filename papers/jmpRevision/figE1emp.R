source('../../share/sysFacLib.R')

dat1a=clean.s2()
dat1b=clean.s5()
dat1bAlt=clean.s5(hard=T)


pdf('figE1emp.pdf',width=10,height=10)
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2,.7,0),cex=1.2)
me.1a=round(doMeanPlotA(dat1a,ylim=c(1.3,1.9)),3)
legend(1,1.5,title="Change in Size",bg="antiquewhite", legend=c("Low","High"),lty=c(1,2),pch=21,pt.bg=c("white","blue"))
mtext(side=3,adj=0,cex=1.3,"A.")
doEmpiricalPlot(dat1a,ylim=c(-.225,.225))
legend(1,.215,title="Classification",bg='antiquewhite',
legend=c('Parallel','Serial','Coactive'),pch=21,pt.bg=c('darkred','white','yellow'))
mtext(side=3,adj=1,cex=1.3,"B.")
me.1b=round(doMeanPlotA(dat1b,ylim=c(1.00,1.3)),3)
#legend(1,1.15,title="Change in Size",legend=c("Low","High"),lty=c(1,2),pch=21,pt.bg=c("white","blue"))
mtext(side=3,adj=0,cex=1.3,"C.")
doEmpiricalPlot(dat1b,ylim=c(-.225,.225))
mtext(side=3,adj=1,cex=1.3,"D.")
dev.off()
