source('../../share/sysFacLib.R')
dat1b=clean.s5()
dat2b=clean.s7()

myRun=function(dat,r.mean=c(.15,.30,.075,.15,.15),r.indv=c(.1,.2,.05,.15,.075))
{
	R=length(r.mean)	
	out=matrix(ncol=5,nrow=R)
	for (r in 1:R) out[r,]=doBayes(dat,r.mean=r.mean[r],r.indv=r.indv[r])$BF
	return(out)
}

out1b=cbind(rep(1,5),myRun(dat1b))
out2b=cbind(rep(1,5),myRun(dat2b))

vals=rbind(out1b,out2b)

library(papaja)

rowMin=apply(vals,1,min)

o=c(1,3,2,5,4,6)

output=matrix(ncol=6,nrow=10,paste('1-to-',round(vals[,o]/rowMin,1),sep=''))
output=ifelse(output=='1-to-1','*',output)
output=ifelse(output=="1-to-Inf", "$ \\approx  0$",output)
colnames(output)=c("Serial","Parallel-1","Parallel-2","Coative-1","Coactive-2","General")
rownames(output)=rep(c("Our Choice (.15,.10)","Double (.30,.15)","Half(.075,.05)","Large Ratio (.15,.15)","Small Ratio (.15,.05)"),2)
apa_table(list(
	"Experiment 1B"=output[1:5,],
	"Experiment 2B"=output[6:10,]),align=c("l",rep("c",6)),caption="Effects of Prior Specification on Bayes Factors.\\label{robust}")

