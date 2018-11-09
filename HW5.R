# header----
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

gradient_ascent<-function(X, thres=1e-03, seed=NULL)
{
  if(!is.null(seed))
  {
    set.seed(seed)
    ITER.SEED<-seed
  }
  else
  {
    
    ITER.SEED <- sample(1:1000, 1)
    set.seed(ITER.SEED)
  }
  store.ll<-list()
  A<-matrix(rnorm(9),nrow=3)
  W<-solve(A)
  Y<-W %*% X
  k<-1
  while(TRUE)
  {
    ayta<-1/(1+k)
    k<-k+1
    W.new <- W + ayta*(gradient(X,Y,A))
    flag<-max(abs((W.new-W)))
    if(flag < thres)
    {
      break
    }
    W<-W.new
    A<-solve(W)
    Y<-W %*% X
    store.ll<-c(store.ll, log_likelihood(W, Y))
  }
  store.ll<-unlist(store.ll)
  return(list(W=W.new, iter_seed = ITER.SEED, iter = k, ll = store.ll))
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
    }
  }
}

log_likelihood<-function(W, Y)
{
  summation<-sum(-log(pi*cosh(pi*Y/2)))
  log_term<-dim(Y)[2]*log(abs(det(W)))
  return(summation+log_term)
}


# Simulated Data----

Y.SIM<-matrix(runif(3*490000), nrow=3)
Y.SIM<-hyperbolic_secant(Y.SIM)
A.SIM<-matrix(runif(9),nrow=3)
W.SIM<-solve(A.SIM)
X.SIM<-W.SIM %*% Y.SIM
W.SIM.GRADIENT.ASCENT<-gradient_ascent(X.SIM)
A.SIM.GRADIENT.ASCENT<-solve(W.SIM.GRADIENT.ASCENT)
Y.SIM.GRADIENT.ASCENT<-A.SIM.GRADIENT.ASCENT %*% X.SIM

png("plots/cov_Y_SIM.png", width=1000, height=1000, units="px")
cov_plot(Y.SIM)
dev.off()

png("plots/cov_Y_GA_SIM.png", width=1000, height=1000, units="px")
cov_plot(Y.SIM.GRADIENT.ASCENT)
dev.off()

# data----
X<-matrix(c(load.wave("mike1.wav"), load.wave("mike2.wav"), load.wave("mike3.wav")), nrow=3, byrow = T)

# operations----

W<-gradient_ascent(X)
cov<-cov_matrix(X)

cov_matrix(X)
png("plots/cov.png", width=1000, height=1000, units="px")
cov_plot(X)
dev.off()

sqcov.inv<-solve(sqrtm(cov))
X.white<-sqcov.inv %*% X
cov_matrix(X.white)
diag(dim(X.white)[1])

W_hat<-gradient_ascent(X.white, seed = 720)
Y.white<- W_hat$W %*% X.white
W<-sqcov.inv %*% W_hat$W
A<-solve(W)
Y <- W %*% X

plot(W_hat$ll, ylab = "Log-Likelihood", xlab = "Index", main = "Log-Likelihood of Iterations")
W_hat$W
W

png("plots/cov_Y_white.png", width=1000, height=1000, units="px")
cov_plot(Y.white)
dev.off()

par(mfrow = c(2,3))
hist(X[1,], main = "mike1.wav")
hist(X[2,], main = "mike2.wav")
hist(X[3,], main = "mike3.wav")
hist(X.white[1,], main = "mike1.wav - Whitened")
hist(X.white[1,], main = "mike2.wav - Whitened")
hist(X.white[1,], main = "mike3.wav - Whitened")
par(mfrow = c(1,1))

par(mfrow = c(2,2))
hist(Y[1,], main = "Source 1")
hist(Y[2,], main = "Source 2")
hist(Y[3,], main = "Source 3")
par(mfrow = c(1,1))

cov_matrix(Y)
png("plots/cov_Y.png", width=1000, height=1000, units="px")
cov_plot(Y)
dev.off()

Y.OUT<-t(apply(Y, 1, function (x) x/(10*max(x))))
save.wave(Y.OUT[1,], "plots/src1.wav")
save.wave(Y.OUT[2,], "plots/src2.wav")
save.wave(Y.OUT[3,], "plots/src3.wav")

# Initializations----
# runif
runif_init<-list()
for(i in 1:10)
{
  print(i)
  W_hat<-gradient_ascent(X.white)
  Y.white<- W_hat$W %*% X.white
  W<-sqcov.inv %*% W_hat$W
  A<-solve(W)
  Y <- W %*% X
  temp_list<-list(iter_seed = W_hat$iter_seed, cov = cov_matrix(Y))
  runif_init<-c(runif_init, temp_list)
  rm(temp_list)
}
rm(runif.max.abs.cov)
runif.max.abs.cov<-list()
i<-1
while(i<=10)
{
  runif.max.abs.cov<-c(runif.max.abs.cov, max(abs(runif_init[[i*2]])))
  i = i+1
}
runif.max.abs.cov<-unlist(runif.max.abs.cov)
runif.cov.index<-which.min(runif.max.abs.cov)*2
runif.iter.seed<-runif_init[[runif.cov.index-1]]

# rnorm
rnorm_init<-list()
for(i in 1:10)
{
  print(i)
  W_hat<-gradient_ascent(X.white)
  Y.white<- W_hat$W %*% X.white
  W<-sqcov.inv %*% W_hat$W
  A<-solve(W)
  Y <- W %*% X
  temp_list<-list(iter_seed = W_hat$iter_seed, cov = cov_matrix(Y))
  rnorm_init<-c(rnorm_init, temp_list)
  rm(temp_list)
}
rnorm.max.abs.cov<-list()
i<-1
while(i<=10)
{
  rnorm.max.abs.cov<-c(rnorm.max.abs.cov, max(abs(rnorm_init[[i*2]])))
  i = i+1
}
rnorm.max.abs.cov<-unlist(rnorm.max.abs.cov)
rnorm.cov.index<-which.min(rnorm.max.abs.cov)*2
rnorm.iter.seed<-rnorm_init[[rnorm.cov.index-1]]
