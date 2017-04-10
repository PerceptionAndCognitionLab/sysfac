source('../../share/sysFacLib.R')

dat2a=clean.s8()
dat2b=clean.s7()


pdf('figE2emp.pdf',width=10,height=10)
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2,.7,0),cex=1.2)
me.2a=round(doMeanPlotB(dat2a,ylim=c(1.0,1.3)),3)
legend(1,1.08,title="Change in 10s Digit",bg="antiquewhite", legend=c("Low","High"),lty=c(1,2),pch=21,pt.bg=c("white","blue"))
mtext(side=3,adj=0,cex=1.3,"A.")
doEmpiricalPlot(dat2a,ylim=c(-.235,.235))
legend(1,.21,title="Classification",bg='antiquewhite',
legend=c('Parallel','Serial','Coactive'),pch=21,pt.bg=c('darkred','white','yellow'))
mtext(side=3,adj=1,cex=1.3,"B.")
me.2b=round(doMeanPlotB(dat2b,ylim=c(.85,1.15)),3)
#legend(1,.95,title="Change in 10s Digit",legend=c("Low","High"),lty=c(1,2),pch=21,pt.bg=c("white","blue"))
mtext(side=3,adj=0,cex=1.3,"C.")
doEmpiricalPlot(dat2b,ylim=c(-.2,.2))
mtext(side=3,adj=1,cex=1.3,"D.")
dev.off()


