pdf('compute.pdf')
par(cex=1.2,mar=c(.1,.1,.1,.1))

library('diagram')

names=c('Parallel-2','General','Coactive-2','Parallel-1','One Effect','Coactive-1','Serial')

myCol=c('peachpuff','thistle','lightgoldenrod')

pos=matrix(nrow=7,ncol=2)
pos[,1]=c(.15,.5,.85,.15,.5,.85,.5)
pos[,2]=c(.75,.85,.75,.4,.5,.4,.15)

M=matrix(nrow=length(names),byrow=F,ncol=length(names),data=0)
M[1,2]=M[2,3]='E'
M[2,5]=M[5,7]='BayesFactor'
M[4,5]=M[5,6]='E'

plotmat(M,pos,name=names,curve=0,box.type="round",box.size=.06,box.prop=1,box.col=myCol[c(1,2,3,1,2,3,2)],arr.length=0)
box()
dev.off()