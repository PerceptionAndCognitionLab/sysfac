z=seq(-.3,.3,.001)
prior=dnorm(z,0,.2)

plotter=function(muD)
{
post=dnorm(z,muD,.03)
plot(z,post,typ='l',lwd=2,col='darkblue',axes=F,xlab=expression(paste("Interaction contrast ",gamma)),ylab="")
lines(z,prior,lwd=2,lty=2,col='darkgreen')
abline(v=0,col='grey')
points(c(0,0),dnorm(c(0,0),c(0,muD),c(.2,.03)),pch=21,cex=1.5,bg=c('darkgreen','darkblue'))
axis(1)
print(dnorm(0,0,.2)/dnorm(0,muD,.03))
}

pdf('savageDickey.pdf',width=10,height=5)
par(mfrow=c(1,2),mar=c(4,1,1,1))
plotter(.03)
mtext(side=3,adj=.2,cex=1.4,"A.",line=-2)
plotter(.08)
mtext(side=3,adj=.2,cex=1.4,"B.",line=-2)
dev.off()