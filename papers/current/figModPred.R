library("spatialfil")
library("tmvtnorm")
library("msm")

sd0 <- .07
eta = .7

gamma <- seq(-.25, .25, .0025)

kern <- convKernel(sigma = 7, k = "gaussian")

nrmlz <- function(mat)
{
  tot <- sum(mat)
  mat/tot
}


#Conditional model specification
norm0 <- function(theta1, theta2, Sigma) dnorm(theta1, 0,Sigma) * dnorm(theta2, 0, Sigma)
norm <- function(theta1, theta2, Sigma) dmvnorm(cbind(theta1, theta2), c(0,0), Sigma)
normT1 <- function(theta1, theta2, Sigma, l, u) dtmvnorm(cbind(theta1, theta2)
                                                   , c(0,0)
                                                   , Sigma
                                                   , lower = rep(l, 2)
                                                   , upper = rep(u, 2))
normT <- function(theta1, theta2, Sigma, l , u){
  dtnorm(theta1, 0, Sigma, lower = l, upper = u) * dtnorm(theta2, 0, Sigma, lower = l, upper = u)
}

Serial <- outer(gamma, gamma, norm0, Sigma = .002)
Serial <- nrmlz(Serial)
Parallel1 <- outer(gamma
                   , gamma
                   , normT1
                   , Sigma = matrix(c(sd0^2, sd0^2.001, sd0^2.001, sd0^2)
                                    , nrow = 2)
                   , l = -Inf
                   , u = 0) 
Parallel1 <- nrmlz(Parallel1)
Parallel2 <- outer(gamma
                   , gamma
                   , normT
                   , sd0
                   , l = -Inf
                   , u = 0)
Parallel2 <- nrmlz(Parallel2)
Coactive1 <- outer(gamma
                   , gamma
                   , normT1
                   , Sigma = matrix(c(sd0^2, sd0^2.001, sd0^2.001, sd0^2)
                                    , nrow = 2)
                   , l = 0
                   , u = Inf)
Coactive1 <- nrmlz(Coactive1)
Coactive2 <- outer(gamma
                   , gamma
                   , normT
                   , sd0
                   , l = 0
                   , u = Inf)
Coactive2 <- nrmlz(Coactive2)
General <- outer(gamma
                 , gamma
                 , norm
                 , Sigma = matrix(c(sd0^2, 0, 0, sd0^2)
                                  , nrow = 2))
General <- nrmlz(General)

#Marginal model specification
GeneralH <- outer(gamma
                  , gamma
                  , norm
                  , Sigma = matrix(c(sd0^2, eta*sd0^2, eta*sd0^2, sd0^2)
                                     , nrow = 2))
GeneralH <- nrmlz(GeneralH)

Parallel2H <- 4 * GeneralH
index <- gamma > 0
Parallel2H[index, ] <- 0
Parallel2H[, index] <- 0
Parallel2H <- nrmlz(Parallel2H)

Coactive2H <- 4 * GeneralH
index <- gamma < 0
Coactive2H[index, ] <- 0
Coactive2H[, index] <- 0
Coactive2H <- nrmlz(Coactive2H)

#Model Predictions
SerialP <- nrmlz(applyFilter(Serial, kern))
Parallel1P <- nrmlz(applyFilter(Parallel1, kern))
Parallel2P <- nrmlz(applyFilter(Parallel2H, kern))
Coactive1P <- nrmlz(applyFilter(Coactive1, kern))
Coactive2P <- nrmlz(applyFilter(Coactive2H, kern))
GeneralP <- nrmlz(applyFilter(GeneralH, kern))

#####Figure
top1 <- max(Parallel1, Coactive2H)
top2 <- max(Parallel1P)
top3 <- max(SerialP)

modFig <- function(mat, par, ylabel, xlabel, main, top, mod){
  image(par
        , par
        , mat
        , col = grey((256:0)/256)
        , zlim = c(0, top)
        , axes = FALSE
        , ylab = ylabel
        , xlab = xlabel
        , bty = "n"
        , main = ""
        , cex.lab = 1.6)
  axis(1, at = seq(-.2, .2, .2), cex.axis = 1.7)
  axis(2, at = seq(-.2, .2, .2), cex.axis = 1.7)
  abline(h = 0, col = "gray70")
  abline(v = 0, col = "gray70")
  mtext(mod, side = 2, line = 4, cex = 1.6)
  mtext(main, side = 3, line = 1, cex = 1.6)
}

pdf('figModPred.pdf',width=10,height=20)
layout(matrix(1:18, ncol = 3), widths = c(.35, .31, .34))

par(mar=c(3,6,3,1), mgp = c(2.4,.9,0))
#models
modFig(Serial, gamma
       , ylabel = expression(paste(gamma[2])), xlabel = ""
       , main = "Conditional Models", top = top1, mod = "Serial")
points(0, 0, pch = 19)

par(mar=c(3,6,1,1), mgp = c(2.4,.9,0))
modFig(Parallel1, gamma
       , ylabel = expression(paste(gamma[2])), xlabel = ""
       , mod = "Parallel 1", top = top1, main = "")

modFig(Parallel2, gamma
       , ylabel = expression(paste(gamma[2])), xlabel = ""
       , mod = "Parallel 2", top = top2, main = "")

modFig(Coactive1, gamma
       , ylabel = expression(paste(gamma[2])), xlabel = ""
       , mod = "Coactive 1", top = top1, main = "")

modFig(Coactive2, gamma
       , ylabel = expression(paste(gamma[2])), xlabel = ""
       , mod = "Coactive 2", top = top2, main = "")

par(mar=c(3.5,6,1,1), mgp = c(2.4,.9,0))

modFig(General, gamma
       , ylabel = expression(paste(gamma[2])), xlabel = expression(paste(gamma[1]))
       , mod = "General", top = top2, main = "")

#marginal
par(mar=c(3,2,3,1), mgp = c(2.4,.9,0))
modFig(Serial, gamma
       , ylabel = "", xlabel = ""
       , main = "Marginal Models", top = top1, mod = "")
points(0, 0, pch = 19)

par(mar=c(3,2,1,1), mgp = c(2.4,.9,0))
modFig(Parallel1, gamma
       , ylabel = "", xlabel = ""
       , mod = "", top = top1, main = "")

modFig(Parallel2H, gamma
       , ylabel = "", xlabel = ""
       , mod = "", top = top2, main = "")

modFig(Coactive1, gamma
       , ylabel = "", xlabel = ""
       , mod = "", top = top1, main = "")

modFig(Coactive2H, gamma
       , ylabel = "", xlabel = ""
       , mod = "", top = top2, main = "")

par(mar=c(3.5,2,1,1), mgp = c(2.4,.9,0))

modFig(GeneralH, gamma
       , ylabel = "", xlabel = expression(paste(gamma[1]))
       , mod = "", top = top2, main = "")

#predictions
par(mar=c(3,4.5,3,1), mgp = c(2.4,.9,0))
modFig(SerialP, gamma
       , ylabel = expression(paste("Observed MIC, ", hat(gamma)[2])), xlabel = ""
       , main = "Predictions", top = top3, mod = "")
par(mar=c(3,4.5,1,1), mgp = c(2.4,.9,0))
modFig(Parallel1P, gamma
       , ylabel = expression(paste("Observed MIC, ", hat(gamma)[2])), xlabel = ""
       , mod = "", top = top2 + .001, main = "")

modFig(Parallel2P, gamma
       , ylabel = expression(paste("Observed MIC, ", hat(gamma)[2])), xlabel = ""
       , mod = "", top = top2, main = "")

modFig(Coactive1P, gamma
       , ylabel = expression(paste("Observed MIC, ", hat(gamma)[2])), xlabel = ""
       , mod = "", top = top2 + .001, main = "")

modFig(Coactive2P, gamma
       , ylabel = expression(paste("Observed MIC, ", hat(gamma)[2])), xlabel = ""
       , mod = "", top = top2, main = "")

par(mar=c(3.5,4.5,1,1), mgp = c(2.4,.9,0))
modFig(GeneralP, gamma
       , ylabel = expression(paste("Observed MIC, ", hat(gamma)[2]))
       , xlabel = expression(paste("Observed MIC, ", hat(gamma)[1]))
       , mod = "", top = top2, main = "")

dev.off()

##Prediction for -.1 and -.09
# Parallel1P[61, 63]/Parallel2P[61, 63]


