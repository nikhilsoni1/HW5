save(list=ls(all=T), file="HW5.RData")
setwd("~/Google Drive/Purdue University/Academics/Sem-3/STAT545/HW5")
load("~/Google Drive/Purdue University/Academics/Sem-3/STAT545/HW5/HW5.RData")
library(audio)
library(VaRES)
library(expm)
# functions----
hyperbolic_secant<-function(x, log=FALSE)
{
  pdf = x
  pdf[log == FALSE] = 1/(pi*cosh(pi * x/2))
  pdf[log == TRUE] = -log(pi) - log(cosh(pi * x/2))
  return(pdf)
}

gradient<-function(X, Y, A)
{
  summation<-(1/dim(X)[2])*(-tanh(Y) %*% t(X))
  return(t(A)+summation)
}

gradient_ascent<-function(X, thres=1e-03, seed=9)
{
  set.seed(seed)
  A<-matrix(runif(9),nrow=3)
  W<-solve(A)
  Y<-t(A) %*% X
  k<-1
  while(TRUE)
  {
    ayta<-1/(1+k)
    k<-k+1
    W<-solve(A)
    Y<-t(A) %*% X
    W.new <- W + ayta*(gradient(X,Y,A))
    A<-solve(W.new)
    flag<-max(abs((W.new-W)))
    if(flag < thres)
    {
      break
    }
  }
  return(W.new)
}

cov_matrix<-function(X)
  {
  mat<-matrix(0 ,nrow = dim(X)[1], ncol = dim(X)[1])
  for (i in 1:dim(X)[1]){
    for (j in 1:dim(X)[1]){
      mat[i,j]<-cov(X[i,], X[j,])
    }
  }
  return(mat)
}

cov_plot<-function(X, r=2, c=2)
{
  par(mfrow = c(r,c))
  for(i in 1:dim(X)[1])
  {
    for (j in 1:dim(X)[1])
    {
      if(i < j)
      {
        plot(X[i,], X[j,], ylab=j, xlab=i)
      }
      else{
      }
    }
  }
}

# data----
X<-matrix(c(load.wave("mike1.wav"), load.wave("mike2.wav"), load.wave("mike3.wav")), nrow=3, byrow = T)

# operations----

W<-gradient_ascent(X)
cov<-cov_matrix(X)

png("plots/cov.png", width=1000, height=1000, units="px")
cov_plot(X)
dev.off()

X.white<-solve(sqrtm(cov)) %*% X
W_hat<-gradient_ascent(X.white)
Y.white<- W_hat %*% X.white
W<-solve(sqrtm(cov)) %*% W_hat
Y <- W %*% X
cov_matrix(Y.white)

png("plots/cov_Y_white.png", width=1000, height=1000, units="px")
cov_plot(Y.white)
dev.off()

png("plots/cov_Y.png", width=1000, height=1000, units="px")
cov_plot(Y)
dev.off()