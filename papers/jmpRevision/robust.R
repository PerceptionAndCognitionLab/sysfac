source('../../share/sysFacLib.R')
dat2b=clean.s7()



base=doBayes(dat2b,r.mean=.15,r.indv=.10)
double=doBayes(dat2b,r.mean=.30,r.indv=.20)
half=doBayes(dat2b,r.mean=.075,r.indv=.05)
ratDouble=doBayes(dat2b,r.mean=.15,r.indv=.15)
ratHalf=doBayes(dat2b,r.mean=.15,r.indv=.075)

save.image(file='robustTemp.RData')
#load('robustTemp.RData')
  
library(papaja)
valsTemp=rbind(base$BF,double$BF,half$BF,ratDouble$BF,ratHalf$BF)
vals=cbind(rep(1,5),valsTemp)
rowMin=apply(vals,1,min)

o=c(1,3,2,5,4,6)

output=matrix(ncol=6,nrow=5,paste('1-to-',round(vals[,o]/rowMin,1),sep=''))
output=ifelse(output=='1-to-1','*',output)
output=ifelse(output=="1-to-Inf", "$ \\approx  0$",output)
colnames(output)=c("Serial","Parallel-1","Parallel-2","Coative-1","Coactive-2","General")
rownames(output)=c("Our Choice (.15,.10)","Double (.30,.15)","Half(.075,.05)","Large Ratio (.15,.15)","Small Ratio (.15,.05)")
apa_table(output,align=c("l",rep("c",6)),caption="Effects of Prior Specification on Bayes Factors.\\label{robust}")