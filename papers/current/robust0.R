source('../../share/sysFacLib.R')
dat1b=clean.s5()

base=doBayes(dat1b,r.mean=.15,r.indv=.10)
double=doBayes(dat1b,r.mean=.30,r.indv=.20)
half=doBayes(dat1b,r.mean=.075,r.indv=.05)
ratDouble=doBayes(dat1b,r.mean=.15,r.indv=.15)
ratHalf=doBayes(dat1b,r.mean=.15,r.indv=.075)


