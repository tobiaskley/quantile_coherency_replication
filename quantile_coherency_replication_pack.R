
###################################################################
## Code to replicate the results in the Barunik and Kley (2018)
###################################################################
##
## Auxiliary Functions definition 
##
## Authors: Jozef Barunik, and Tobias Kley (2018)
###################################################################

library(alphahull)
library(astsa)
library(copula)
library(latex2exp)
library(mAr)        # for simulation of VAR 
library(pbivnorm)
library(quantreg)
library(quantspec)
library(rugarch)
library(vars)       # for estimation of VAR
library(zoo)

# Define various Data Generating Processes

indep <- function(n, rho = 0) {
  Y1 <- rnorm(n)
  Y2 <- rho * Y1 + sqrt(1 - rho^2) * rnorm(n)
  res <- matrix(c(Y1, Y2), ncol = 2)
  return(as.matrix(res))
}

indepq <- function(n) {  
  y1 <- rnorm(n)
  y2 <- (y1^2 - 1) / 3
  return(cbind(y1, y2))  
}

indepqshift <- function(n) {
  y1 <- rnorm(n+1)
  y2 <- (y1[2:(n+1)]^2 - 1) / 3
  return(cbind(y1[1:n], y2))
}

transl <- function (rho,tau1,tau2){
  p12 <- pbivnorm(qnorm(tau1, 0, 1), qnorm(tau2, 0, 1), rho)
  return( (p12-tau1*tau2)/(sqrt(tau1*(1-tau1))*sqrt(tau2*(1-tau2))) )
}

indepR <- function(n, rho=0.6) {
  Y1 <- rnorm(n)
  Y2 <- rho * Y1 + sqrt(1 - rho^2) * rnorm(n)
  res <- matrix(c(Y1, Y2), ncol = 2)
  return(as.matrix(res))
}



quar31<-function (n, 
                  th11 = function(u) { 0 * ((u - 0.5))}, 
                  th12 = function(u) { 1.2 * ((u - 0.5))},
                  th21 = function(u) { 1.2 * ((u - 0.5))}, 
                  th22 = function(u) { 0 * ((u - 0.5))},
                  
                  ga11 = function(u) { 0 * ((u - 0.5))},
                  ga12 = function(u) { 0 * ((u - 0.5))},
                  ga21 = function(u) { 0 * ((u - 0.5))},
                  ga22 = function(u) { 0 * ((u - 0.5))},
                  overhead = 1000, th01 = qnorm,th02 = qnorm,rho = 0) 
{
  Y1 <- rep(0, n + overhead)
  Y2 <- rep(0, n + overhead)
  Y1[1:2] <- th01(runif(2))
  Y2[1:2] <- th02(runif(2))
  for (t in 3:(n + overhead)) {
    U1 <- runif(1)
    U2 <- runif(1)
    Y1[t] <- th11(U1) * Y1[t - 1] + th12(U1) * Y2[t - 1] +
             ga11(U1) * Y1[t - 2] + ga12(U1) * Y2[t - 2] + th01(U1)
    Y2[t] <- th21(U2) * Y1[t - 1]  + th22(U2) * Y2[t - 1] +
             ga21(U2) * Y1[t - 2]  + ga22(U2) * Y2[t - 2] + th02(U2)
  }
  Y <- cbind(Y1[(overhead + 1):(overhead + n)],
             Y2[(overhead + 1):(overhead + n)])
  return(as.matrix(Y))
}

quar32<-function (n, 
                  th11 = function(u) { 0 * ((u - 0.5))}, 
                  th12 = function(u) { 0 * ((u - 0.5))},
                  th21 = function(u) { 0 * ((u - 0.5))}, 
                  th22 = function(u) { 0 * ((u - 0.5))},
                  
                  ga11 = function(u) { 0 * ((u - 0.5))},
                  ga12 = function(u) { 1.2 * ((u - 0.5))},
                  ga21 = function(u) { 1.2 * ((u - 0.5))},
                  ga22 = function(u) { 0 * ((u - 0.5))},
                  overhead = 1000, th01 = qnorm, th02 = qnorm, rho = 0) 
{
  Y1 <- rep(0, n + overhead)
  Y2 <- rep(0, n + overhead)
  Y1[1:2] <- th01(runif(2))
  Y2[1:2] <- th02(runif(2))
  for (t in 3:(n + overhead)) {
    U1 <- runif(1)
    U2 <- runif(1)
    Y1[t] <- th11(U1) * Y1[t - 1] + th12(U1) * Y2[t - 1] +
             ga11(U1) * Y1[t - 2] + ga12(U1) * Y2[t - 2] + th01(U1)
    Y2[t] <- th21(U2) * Y1[t - 1]  + th22(U2) * Y2[t - 1] +
             ga21(U2) * Y1[t - 2]  + ga22(U2) * Y2[t - 2] + th02(U2)
  }
  Y <- cbind(Y1[(overhead + 1):(overhead + n)],
             Y2[(overhead + 1):(overhead + n)])
  return(as.matrix(Y))
}

quar33<-function (n, 
                  th11 = function(u) { 0 * ((u - 0.5))}, 
                  th12 = function(u) { 0 * ((u - 0.5))},
                  th21 = function(u) { 0 * ((u - 0.5))}, 
                  th22 = function(u) { 0 * ((u - 0.5))},
                  
                  ga11 = function(u) { 0 * ((u - 0.5))},
                  ga12 = function(u) { 1.2 * ((u - 0.5))},
                  ga21 = function(u) { 1.2 * ((u - 0.5))},
                  ga22 = function(u) { 0 * ((u - 0.5))},
                  overhead = 1000, th01 = qnorm, th02 = qnorm, rho = 0) 
{
  Y1 <- rep(0, n + overhead)
  Y2 <- rep(0, n + overhead)
  Y1[1:3] <- th01(runif(2))
  Y2[1:3] <- th02(runif(2))
  for (t in 4:(n + overhead)) {
    U1 <- runif(1)
    U2 <- runif(1)
    Y1[t] <- th11(U1) * Y1[t - 1] + th12(U1) * Y2[t - 1] +
             ga11(U1) * Y1[t - 3] + ga12(U1) * Y2[t - 3] + th01(U1)
    Y2[t] <- th21(U2) * Y1[t - 1]  + th22(U2) * Y2[t - 1] +
             ga21(U2) * Y1[t - 3]  + ga22(U2) * Y2[t - 3] + th02(U2)
  }
  Y <- cbind(Y1[(overhead + 1):(overhead + n)],
          Y2[(overhead + 1):(overhead + n)])
  return(as.matrix(Y))
}


###############################################################################
# The function below computes
#
#		(1) the traditional spectrum
#		(2) the quantile spectrum
#	
#	for a VAR(1) process with
#
# 	X_t = A X_{t-1} + eps_t
#
# 	where eps_t is white noise with Var(eps_t) = Sigma
# 
# Author: Tobias Kley
###############################################################################

getTradQuantCoherencyInVAR1 <- function(A, Sigma, levels, maxLag = 127) {
  
  D <- ncol(A)
  
  Gamma0 <- matrix(solve(diag(rep(1, D^2)) - A %x% A) %*% as.vector(Sigma),
                   ncol = D)
  
  ## determine all (traditional) Gamma_k
  
  Gamma <- array(0, dim = c(maxLag, D, D))
  
  Gamma[1,,] <- A %*% Gamma0 
  for (i in 2:maxLag) {
    Gamma[i,,] <- A %*% Gamma[i-1,,]
  }
  
  ## next compute copula cross-covariances
  
  K <- length(levels)
  
  qGamma0 <- array(0, dim = c(D, K, D, K))
  qGamma <- array(0, dim = c(maxLag, D, K, D, K)) 
  
  Rho0 <- array(0, dim = c(D, D))
  Rho <- array(0, dim = c(maxLag, D, D))
  
  for (i1 in 1:D) {
    for (i2 in 1:D) {
      Rho0[i1,i2] <- Gamma0[i1,i2] / sqrt(Gamma0[i1,i1] * Gamma0[i2,i2])  
      Rho[,i1,i2] <- Gamma[,i1,i2] / sqrt(Gamma0[i1,i1] * Gamma0[i2,i2]) 
    }
  }
  
  for (i1 in 1:D) {
    for (k1 in 1:K) {
      for (i2 in 1:D) {
        for (k2 in 1:K) {
          tau1 <- levels[k1]
          tau2 <- levels[k2]
          
          rho <- Rho0[i1,i2]
          p12 <- pbivnorm(qnorm(tau1,0,1), qnorm(tau2,0,1), rho)
          qGamma0[i1,k1,i2,k2] <- (p12-tau1*tau2) /
              (sqrt(tau1*(1-tau1))*sqrt(tau2*(1-tau2)))
          
          rho <- Rho[,i1,i2]
          p12 <- pbivnorm(qnorm(tau1,0,1), qnorm(tau2,0,1), rho)
          qGamma[,i1,k1,i2,k2] <- (p12-tau1*tau2) /
              (sqrt(tau1*(1-tau1))*sqrt(tau2*(1-tau2)))       
        }
      }
    }
  }
  
  ## now compute the spectra
  # (1) traditional spectra
  
  tradSpec <- array(0, dim=c(maxLag+1, D, D))
  quantSpec <- array(0, dim=c(maxLag+1, D, K, D, K))
  
  tradCoh <- array(0, dim=c(maxLag+1, D, D))
  quantCoh <- array(0, dim=c(maxLag+1, D, K, D, K))
  
  
  for (i1 in 1:D) {
    for (i2 in 1:D) {
      AA <- fft(c(Gamma0[i1,i2], Gamma[,i1,i2]))
      BB <- Conj(fft(c(Gamma0[i2,i1], Gamma[,i2,i1]))) - Gamma0[i2,i1]
      tradSpec[,i1,i2] <- AA + BB
      
      for (k1 in 1:K) {
        for (k2 in 1:K) {
          AA <- fft(c(qGamma0[i1,k1,i2,k2], qGamma[,i1,k1,i2,k2]))
          BB <- Conj(fft(c(qGamma0[i1,k1,i2,k2], qGamma[,i1,k1,i2,k2]))) -
              qGamma0[i1,k1,i2,k2]
          quantSpec[,i1,k1,i2,k2] <- AA + BB
        }
      }
    }
  }
  tradSpec <- tradSpec / (2*pi)
  quantSpec <- quantSpec / (2*pi)
  
  ## now compute coherency
  
  tradCoh <- array(0, dim=c(maxLag+1, D, D))
  quantCoh <- array(0, dim=c(maxLag+1, D, K, D, K))
  
  for (i1 in 1:D) {
    for (i2 in 1:D) {
      tradCoh[,i1,i2] <- tradSpec[,i1,i2] /
          sqrt(Re(tradSpec[,i1,i1]) * Re(tradSpec[,i2,i2]))
      for (k1 in 1:K) {
        for (k2 in 1:K) {
          
          quantCoh[,i1,k1,i2,k2] <- quantSpec[,i1,k1,i2,k2] /
              sqrt(Re(quantSpec[,i1,k1,i1,k1]) * Re(quantSpec[,i2,k2,i2,k2]))
        }
      }
    }
  }
  return( list(tradSpec = tradSpec, quantSpec = quantSpec,
               tradCoh = tradCoh, quantCoh = quantCoh) )
  
}


# functiono to determine points on the alpha convex hull:
aConvHull <- function(X, a = 0.1) {
  Xu <- unique(round(X,5))
  
  alsh <- ashape(x=Xu[,1], y=Xu[,2], alpha = a)
  
  Pts <- alsh$edges[,3:6]
  
  Y <-rbind(Pts[1,1:2], Pts[1,3:4])
  Pts <- Pts[-1,]
  
  repeat {
    
    if (sum((Pts[,1] == Y[nrow(Y),1]) + (Pts[,2] == Y[nrow(Y),2])) > 1) {
      j <- which(Pts[,1] == Y[nrow(Y),1] & Pts[,2] == Y[nrow(Y),2])
      Y <- rbind(Y, Pts[j[1],3:4])
    } else if (sum((Pts[,3] == Y[nrow(Y),1]) + (Pts[,4] == Y[nrow(Y),2])) > 1) {
      j <- which(Pts[,3] == Y[nrow(Y),1] & Pts[,4] == Y[nrow(Y),2])
      Y <- rbind(Y, Pts[j[1],1:2])
    } else {
      stop("something is wrong here!")
    }
    if (length(j) > 1) {break}
    
    Pts <- Pts[-j[1],]
    if (!is.matrix(Pts)) {break}
  }
  
  return (rbind(Y,Y[1,])) 
}
# END: (1) determine points on the alpha convex hull:




###############################################################################
## DEFINE AUXILIARY FUNCTIONS for estimation of QVAR
###############################################################################

##simulate a QVAR(1), WITHOUT the spacial component
rqvar1 <- function (n, th10, th11, th12, th20, th21, th22,
                    myCop, overhead = 1000) 
{
  Y1 <- rep(0, n + overhead)
  Y2 <- rep(0, n + overhead)
  Y1[1:2] <- th10(runif(2))
  Y2[1:2] <- th20(runif(2))
  
  innov = rCopula(n + overhead, myCop)
  
  U01 <- innov[,1]
  U02 <- innov[,2]
  
  for (t in 3:(n + overhead)) {
    U1 <- U01[t]
    U2 <- U02[t]
    Y1[t] <- th10(U1) + th11(U1) * Y1[t - 1] + th12(U1) * Y2[t - 1]
    Y2[t] <- th20(U2) + th21(U2) * Y1[t - 1] + th22(U2) * Y2[t - 1]
  }
  Y <- cbind(Y1[(overhead + 1):(overhead + n)],
             Y2[(overhead + 1):(overhead + n)])
  return(as.matrix(Y))
}

## simulate a QVAR(1), WITH spatial components in Y_1 and Y_2
rqvar2 <- function (n, th10, th120, th111, th121, th20, th210, th211, th221,
                    myCop, overhead = 1000) 
{
  Y1 <- rep(0, n + overhead)
  Y2 <- rep(0, n + overhead)
  Y1[1:2] <- 0 
  Y2[1:2] <- 0 
  
  innov = rCopula(n + overhead, myCop)
  
  U01 <- innov[,1]
  U02 <- innov[,2]
  
  for (t in 3:(n + overhead)) {
    U1 <- U01[t]
    U2 <- U02[t]
    
    d1 <- 1 - th120(U1) * th210(U2)
    b1 <- th10(U1) + th120(U1) * th20(U2)
    a11 <- th111(U1) + th120(U1) * th211(U2)
    a12 <- th121(U1) + th120(U1) * th221(U2)
    Y1[t] <- (b1 + a11 * Y1[t - 1] +  a12 * Y2[t - 1]) / d1
    
    d2 <- 1 - th210(U2) * th120(U1)
    b2 <- th20(U2) + th210(U2) * th10(U1)
    a21 <- th211(U2) + th210(U2) * th111(U1)
    a22 <- th221(U2) + th210(U2) * th121(U1)
    Y2[t] <- (b2 + a21 * Y1[t - 1] + a22 * Y2[t - 1]) / d2
  }
  Y=cbind(Y1[(overhead + 1):(overhead + n)],
          Y2[(overhead + 1):(overhead + n)])
  return(as.matrix(Y))
}

## simulate a QVAR(1), WITH a spatial component in Y_1
rqvar2a <- function (n, th10, th120, th111, th121, th20, th211, th221,
                     myCop, overhead = 1000) 
{
  Y1 <- rep(0, n + overhead)
  Y2 <- rep(0, n + overhead)
  Y1[1:2] <- 0 
  Y2[1:2] <- 0 
  
  innov = rCopula(n + overhead, myCop)
  
  U01 <- innov[, 1]
  U02 <- innov[, 2]
  
  for (t in 3:(n + overhead)) {
    U1 <- U01[t]
    U2 <- U02[t]
    Y2[t] <- th20(U2) + th211(U2) * Y1[t - 1] + th221(U2) * Y2[t - 1]
    Y1[t] <- th10(U1) + th120(U1) * Y2[t]
    + th111(U1) * Y1[t - 1] + th121(U1) * Y2[t - 1]
  }
  Y <- cbind(Y1[(overhead + 1):(overhead + n)],
             Y2[(overhead + 1):(overhead + n)])
  return(as.matrix(Y))
}

## simulate a QVAR(1), WITH a spatial component in Y_2
rqvar2b <- function (n, th10, th111, th121, th20, th210, th211, th221,
                     myCop, overhead = 1000) 
{
  Y1 <- rep(0, n + overhead)
  Y2 <- rep(0, n + overhead)
  Y1[1:2] <- 0 
  Y2[1:2] <- 0 
  
  innov = rCopula(n + overhead, myCop)
  
  U01 <- innov[,1]
  U02 <- innov[,2]
  
  for (t in 3:(n + overhead)) {
    U1 <- U01[t]
    U2 <- U02[t]
    Y1[t] <- th10(U1) + th111(U1) * Y1[t - 1] + th121(U1) * Y2[t - 1]
    Y2[t] <- th20(U2) + th210(U2) * Y1[t]
    + th211(U2) * Y1[t - 1] + th221(U2) * Y2[t - 1]
  }
  Y <- cbind(Y1[(overhead + 1):(overhead + n)],
             Y2[(overhead + 1):(overhead + n)])
  return(as.matrix(Y))
}
