library(papaja)
load('bayesTemp.RData')
valsTemp=rbind(out1a$BF,out1b$BF,out2a$BF,out2b$BF)
vals=cbind(rep(1,4),valsTemp)
rowMin=apply(vals,1,min)

o=c(1,3,2,5,4,6)

output=matrix(ncol=6,nrow=4,paste('1-to-',round(vals[,o]/rowMin,1),sep=''))
output=ifelse(output=='1-to-1','*',output)
output=ifelse(output=="1-to-Inf", "$ \\approx  0$",output)
colnames(output)=c("Serial","Parallel-1","Parallel-2","Coative-1","Coactive-2","General")
rownames(output)=c("Expt. 1A","Expt. 1B", "Expt. 2A","Expt. 2B")
apa_table(output,align=c("l",rep("c",6)),caption="Bayes factor values.\\label{bf}")