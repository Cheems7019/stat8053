source("https://raw.githubusercontent.com/hyunseungkang/invalidIV/master/TSHT.R")
source("https://raw.githubusercontent.com/hyunseungkang/invalidIV/master/JMCode.R")
source("/Users/lichen/Dropbox/post_selection_inference/Code/PoSI.R",chdir = TRUE)

### Comparison among PoSI1 and PoSI

data(Boston, package="MASS")
UL.Boston1 <- PoSI(Boston[,-14]) 
summary(UL.Boston1)

UL.Boston2 <- PoSIK(Boston[,13,drop=FALSE],Boston[,-c(13,14)])
summary(UL.Boston2)

### Toy example (no covariates)

library(MASS)
set.seed(0)
n = 500; L = 10; s = 3; 
alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,L))
epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
Z = matrix(rnorm(n*L),n,L)

epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
D = Z %*% gamma + epsilon[,1]
Y = Z %*% alpha + D * beta + epsilon[,2]

tsht.example <- TSHT(Y,D,Z)

dmodz <- lm(D~Z+0)
hat_D <- dmodz$fitted.values
hat_D <- matrix(hat_D,ncol=1)
UL.example <- PoSIK(hat_D,Z,modelSZ = 1:floor((ncol(Z)-1)/2))
summary(UL.example)


### Toy example (with covariates)

library(MASS)
set.seed(0)
n = 500; L = 10; s = 3; p <- 3;
alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,L)); beta2 = rep(1,3);
epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
Z = matrix(rnorm(n*L),n,L)
X = matrix(rnorm(n*p),n,p)

epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
D = Z %*% gamma + epsilon[,1]
Y = Z %*% alpha + D * beta + X %*% beta2 + epsilon[,2]

tsht.example <- TSHT(Y,D,Z,X)

dmodz <- lm(D~Z+0)
hat_D <- dmodz$fitted.values
hat_D <- matrix(hat_D,ncol=1)
hat_D <- cbind(hat_D,X)
UL.example <- PoSIK1(hat_D,Z,modelSZ = 1:floor((ncol(Z)-1)/2))
summary(UL.example)





