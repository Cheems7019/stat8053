library(MASS)
source("https://raw.githubusercontent.com/hyunseungkang/invalidIV/master/TSHT.R")
source("https://raw.githubusercontent.com/hyunseungkang/invalidIV/master/JMCode.R")
source("PoSI.R",chdir = TRUE)

data.generation <- function(n,beta,gamma0,tau,px){
  pi.star <- c(rep(0,6),rep(tau*gamma0,2),-0.5,-1)
  gamma.star <- rep(gamma0,10)
  psi.star <- seq(from=0.5,by=0.1,length.out=px)
  phi.star <- seq(from=1.1,by=0.1,length.out=px)
  
  pz <- 10; p <- pz+px
  ZXSigma <- outer(1:p,1:p,function(x,y){0.5^(abs(x-y))})
  ZX <- mvrnorm(n,mu=rep(0,p),Sigma = ZXSigma)
  Z <- ZX[,1:pz,drop=FALSE]; X <- ZX[,(pz+1):p]
  
  epsilonSigma <- matrix(c(1,0.8,0.8,1),2,2)
  epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
  
  D = Z %*% gamma.star + X %*% psi.star + epsilon[,1]
  Y =  D * beta + Z %*% pi.star + X %*% phi.star + epsilon[,2]
  
  coef <- c(pi.star,phi.star,beta)
  convert <- rbind(diag(p),c(gamma.star,psi.star))
  ZXDSigma <- convert%*%ZXSigma%*%t(convert)
  return(list(Y=Y,Z=Z,D=D,X=X,coef=coef,ZXDSigma=ZXDSigma))
}

simu.single <- function(n,beta,gamma0,tau,px,Nsim.posi){
  ## Data generation
  dat <- data.generation(n,beta,gamma0,tau,px)
  Y <- dat$Y; Z <- dat$Z; D <- dat$D; X <- dat$X; coef <- dat$coef; ZXDSigma <- dat$ZXDSigma;
  
  ## TSHT-CI construction
  tsht.result <- TSHT(Y,D,Z,X)
  tsht.ci <- tsht.result$ci
  validiv <- tsht.result$VHat
  
  
  ## TSHT-POSI-CI construction
  dmodzx <- lm(D~Z+X+0)
  hat_D <- dmodzx$fitted.values
  hat_D <- matrix(hat_D,ncol=1)
  hat_D <- cbind(hat_D,X)
  posi.result <- PoSIK1(hat_D,Z,verbose = 0,modelSZ = 1:floor((ncol(Z)-1)/2),Nsim=Nsim.posi)
  posi.result <- summary(posi.result)
  posi.ci <- c(tsht.result$betaHat-posi.result[1,1]*sqrt(tsht.result$betaVarHat/n),tsht.result$betaHat+posi.result[1,1]*sqrt(tsht.result$betaVarHat/n))
  
  ## beta* in the sub-model
  idex <- nrow(ZXDSigma) - length(validiv)
  incre <- solve(ZXDSigma[-validiv,-validiv])%*%ZXDSigma[-validiv,validiv]%*%coef[validiv]
  beta.sub <- beta + incre[idex]
  
  return(c(beta,tsht.ci,posi.ci,beta.sub))
}

simu <- function(n,beta,gamma0,tau,px,Nsim.posi,Nsim){
  result.mat <- matrix(0,Nsim,6)
  for (i in 1:Nsim){
    result.single <- simu.single(n,beta,gamma0,tau,px,Nsim.posi)
    result.mat[i,] <- result.single
  }
  return(result.mat)
}

summary.simu <- function(mat){
  n <- nrow(mat)
  cover.tsht <- sum((mat[,2]<mat[1,1])*(mat[,3]>mat[1,1]))/n
  cover.posi <- sum((mat[,4]<mat[1,1])*(mat[,5]>mat[1,1]))/n
  cover.tsht.sub <- sum((mat[,2]<mat[,6])*(mat[,3]>mat[,6]))/n
  cover.posi.sub <- sum((mat[,4]<mat[,6])*(mat[,5]>mat[,6]))/n
  len.tsht <- mean(mat[,3]-mat[,2])
  len.posi <- mean(mat[,5]-mat[,4])
  return(c(cover.tsht,cover.posi,cover.tsht.sub,cover.posi.sub,len.tsht,len.posi))
}

## Summary of Results
results <- rbind(summary.simu(mat1),summary.simu(mat2),summary.simu(mat3))

## Setting 1
n = 500; beta = 1; gamma0 = 0.5; tau = 0.4; px=10; Nsim.posi = 10000; Nsim=1000
set.seed(0)
mat1 <- simu(n,beta,gamma0,tau,px,Nsim.posi,Nsim)
write.csv(mat1,"./simu_results/mat1.csv")

## Setting 2
n = 1000; beta = 1; gamma0 = 0.5; tau = 0.4; px=10; Nsim.posi = 10000; Nsim=1000
set.seed(0)
mat2 <- simu(n,beta,gamma0,tau,px,Nsim.posi,Nsim)
write.csv(mat2,"./simu_results/mat2.csv")

## Setting 3
n = 1500; beta = 1; gamma0 = 0.5; tau = 0.4; px=10; Nsim.posi = 10000; Nsim=1000
set.seed(0)
mat3 <- simu(n,beta,gamma0,tau,px,Nsim.posi,Nsim)
write.csv(mat3,"./simu_results/mat3.csv")



