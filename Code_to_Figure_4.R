
###################################################################
## Code to replicate the results in the Barunik and Kley (2018)
###################################################################
##
## Code to Figure 4
##
## Authors: Jozef Barunik, and Tobias Kley (2018)
###################################################################

# NOTE: this code produces a pdf figure that is saved to the current 
#       working folder


# load the auxiliary files
source("quantile_coherency_replication_pack.R")

# set seed for RNG
set.seed(1234)

# these packages are only needed for this figure
library(np)
library(quantilogram)

# Load the data
ffdata <- read.csv(file("12_Industry_Portfolios_Daily1.csv"),
                   header = TRUE, sep = ";")
ffdatamkt <- read.csv(file("F-F_Research_Data_Factors_daily.csv"),
                      header = TRUE, sep = ";")

datY10 <- ffdata[, 2]
datY20 <- ffdatamkt[, 2]

# remove mean and heteroscedasticity
model <- ugarchspec(variance.model = list(garchOrder = c(1,1)),
                    mean.model = list(armaOrder = c(1,1), include.mean = T),
                    distribution.model = "norm")

fit1  <- ugarchfit(model, data = datY10)
fit2  <- ugarchfit(model, data = datY20)
datY1 <- fit1@fit$residuals / fit1@fit$sigma
datY2 <- fit2@fit$residuals / fit2@fit$sigma
datY1 <- datY1 - mean(datY1)
datY2 <- datY2 - mean(datY2)

Y <- matrix(c(datY1, datY2), ncol=2)
n <- dim(Y)[1]

## ComputeCross-quantilogram
# Here we use the adapted version of the original code provided in the 
# quantilogram package

##====================
## setup values
##====================
## the maximum lag orders 
Kmax <- 60  

## quantile ranges
vecA1  <- as.matrix(c(0.05, 0.5, 0.95))
scaA2  <- as.matrix(c(0.05, 0.5, 0.95))
A1size <- length(vecA1)


## the number of repetition for the stationary bootstrap
Bsize  <- 1000

## the significance level 
sigLev <- 0.05

## the triming parameter for the self-normalization
scaW   <- 0.1
## - critical values
data(CV.SelfN.trim.0.10)
vecCV.SN <- as.matrix(CV.SelfN.trim.0.10$sig0.05)

##=========================
## data 
##=========================

DATA <- Y

##========================================================================
## Cross-Q
##========================================================================

#==========================
# Estimation and Bootstrap
#==========================
matCRQ  <- matrix(0, Kmax, A1size)
matCI.L <- matrix(0, Kmax, A1size)
matCI.R <- matrix(0, Kmax, A1size)

## Cross Q
for (j in 1:A1size){
  print(j)
  vecA <- matrix( c(vecA1[j], scaA2[j]), 2, 1)
  
  for (k in 1:Kmax){
    
    RES          <- CrossQ.StatBoot.OptAve(DATA, vecA, k, Bsize, sigLev)
    matCI.L[k,j] <- RES$vecCV[1]
    matCI.R[k,j] <- RES$vecCV[2]
    matCRQ[k,j]  <- RES$vCRQ
  }
}



##==============
## FIGURE 4
##==============
## x-axis
vec.lag <- as.matrix(seq(1,Kmax))

## Cross-Quantilogram with the SB-resample critical values 

name <- c("a","b","c")
for (j in 1:A1size){
  
  pdf(file = paste("fig_04", name[j], ".pdf", sep = ""),
      width = 4, height = 5.3)
    par(mar = c(4, 2, 5.5, 0.5) + 0.1)
  
    vecCRQ  <- matCRQ[, j, drop = FALSE]
    vecCI.L <- matCI.L[, j, drop = FALSE]
    vecCI.R <- matCI.R[, j, drop = FALSE]
    
    plot(vec.lag, vecCRQ, type = "h",
         xlab = "Lag", ylab = "",
         ylim = c(-0.1,0.1))
     
    ## line for 0
    abline(h = 0)
    
    ## add lines for CI
    lines(vec.lag, vecCI.L, lwd = 1.5, lty = 2 , col = "black") 
    lines(vec.lag, vecCI.R, lwd = 1.5, lty = 2 , col = "black") 
  dev.off()
}
