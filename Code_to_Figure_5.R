###################################################################
## Code to replicate the results in the Barunik and Kley (2018)
###################################################################
##
## This script is to estimate a VAR(1) model  to then
## simulate the corresponding quantile coherencies.
##
## Code also replicates Figures 5
##
## Authors: Jozef Barunik, and Tobias Kley (2018)
###################################################################

# NOTE: this code produces nine pdf figures:
#
#         - fig_05a1.pdf
#         - fig_05a2.pdf
#         - fig_05a3.pdf
#         - fig_05b1.pdf
#         - fig_05b2.pdf
#         - fig_05b3.pdf
#         - fig_05c1.pdf
#         - fig_05c2.pdf
#         - fig_05c3.pdf
#
#       that are saved to the current working folder.

# NOTE: simulation of the corresponding quantities is lengthy, 

# load the auxiliary files
source("quantile_coherency_replication_pack.R")

# set seed for RNG
set.seed(1234)

R <- 1000

le <- 1024
quantile <- seq(0.05, 0.95, 0.05)
which_taus <- (2:98) / 100


# load data
ffdata <- read.csv(file("12_Industry_Portfolios_Daily1.csv"),
                   header = TRUE, sep = ";")
ffdatamkt <- read.csv(file("F-F_Research_Data_Factors_daily.csv"),
                      header = TRUE, sep = ";")

###############################################################################
## PREPROCESS THE DATA
###############################################################################

datY10 <- ffdata[,2]
datY20 <- ffdatamkt[,2]

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
colnames(Y) <- c("V1", "V2")
n <- dim(Y)[1]

###############################################################################
## FIT VAR(1)
###############################################################################

fit <- VAR(Y, p = 1, type = "const")

rho <- cor(fit$varresult$V1$residuals, fit$varresult$V2$residuals)
v1 <- var(fit$varresult$V1$residuals)
v2 <- var(fit$varresult$V2$residuals)

th10 <- function(u) {qnorm(u)}
th11 <- function(u) {fit$varresult$V1$coefficients[1]}
th12 <- function(u) {fit$varresult$V1$coefficients[2]}

th20 <- function(u) {qnorm(u)}
th21 <- function(u) {fit$varresult$V2$coefficients[1]}
th22 <- function(u) {fit$varresult$V2$coefficients[2]}

###############################################################################
##  Simulate the qSD when the innovations are Gaussian
###############################################################################

myCop <- normalCopula(rho, dim = 2)

model <- function (n) {
  rqvar1(n, th10, th11, th12, th20, th21, th22,
      myCop, overhead = 1000)
}

qSD <- quantileSD(le, ts = model, levels.1 = quantile, R = R)
save(qSD, file = "VAR1_Gaussian.Rdata")

###############################################################################
## Simulate the qSD when the innovations are Gumbel
###############################################################################

myCop <- gumbelCopula(2.7, dim = 2)
model <- function (n) {
  rqvar1(n, th10, th11, th12, th20, th21, th22,
      myCop, overhead = 1000)
}
qSD <- quantileSD(le, ts = model, levels.1 = quantile, R = R)
save(qSD, file = "VAR1_Gumbel.Rdata")

###############################################################################
## Simulate the qSD when the innovations are Clayton
###############################################################################

myCop <- claytonCopula(4, dim = 2)
model <- function (n) {
  rqvar1(n, th10, th11, th12, th20, th21, th22,
      myCop, overhead = 1000)
}
qSD <- quantileSD(le, ts = model, levels.1 = quantile, R = R)
save(qSD, file = "VAR1_Clayton.Rdata")



###############################################################################
## Generate Figure 5
###############################################################################

# clean global environment
rm(list = ls())

#set thickness
tht <- 1.5

which_taus <- (2:98)/100

# parameters for the pdfs
le <- 1024
quantile <- seq(0.05, 0.95, 0.05)
quantilelegend <- c("0.5 | 0.5", "0.05 | 0.05", "0.95 | 0.95")

lim1 <- 0
lim2 <- 1
lim3 <- -0.4
lim4 <- 0.2

###############################################################################
## VAR(1) Gaussian
###############################################################################

load("VAR1_Gaussian.Rdata")

# (ii) Quantile coherency (for 4 different combinations of parameters
#                              and all frequencies)

V <- getValues(qSD, levels.1 = c(0.5, 0.05, 0.95),
                    levels.2 = c(0.5, 0.05, 0.95))
V <- V[1:(dim(V)[1]/2),,,,]

frequencies <- (0:le)/(le+1)
freq <- frequencies[which(frequencies != 0)][1:(le/2)]

pdf("fig_05a1.pdf", width = 4, height = 3.5)
  par(mar = c(4, 2, 2, 0.5) + 0.1)
  
  plot(x = freq, xlim = c(0, 0.5), ylim = c(lim1, lim2), type = "l",
       xlab = expression(omega / 2 * pi), ylab = "")
  for (i1 in 1:3) {
    lines(x = freq, y = Re(V[,1,i1,2,i1]/sqrt(V[,1,i1,1,i1]*V[,2,i1,2,i1])),
          col = "black", lty = i1, lwd = tht)
  }
  axis(side = 3, at = c(1/5, 1/22, 1/250), labels = c("W", "M", "Y"))
  abline(v=c(1/5, 1/22, 1/250), col = "gray")
  legend("bottom", inset = .03, "center",
         c("0.5 | 0.5", "0.05 | 0.05", "0.95 | 0.95"),
         cex = 0.65, lwd = c(1, 1, 1), lty = c(1:3),
         horiz = TRUE, bg = "white")

dev.off()

pdf("fig_05a2.pdf", width = 4, height = 3.5)
  par(mar = c(4, 2, 2, 0.5) + 0.1)
  plot(x = freq, xlim = c(0, 0.5), ylim = c(-0.2, lim2), type = "l",
       xlab = expression(omega / 2 * pi), ylab = "")
  lines(x = freq, y = Re(V[,1,1,2,3]/sqrt(V[,1,1,1,1]*V[,2,3,2,3])),
        col = "black", lty = i1, lwd = tht) 
  axis(side = 3, at = c(1/5, 1/22, 1/250), labels = c("W", "M", "Y"))
  abline(v = c(1/5, 1/22, 1/250), col = "gray")
  legend("bottom", inset = .03, "center", c("0.05 | 0.95"),
         cex = 0.8, lwd = c(1, 1, 1), lty = c(1:3),
         horiz = TRUE, bg = "white")
dev.off()

# (iii) Quantile coherency (for 3 frequencies and all taus)

V <- getValues(qSD, frequencies = 2*pi * c(1/5, 1/22, 1/250))

pdf("fig_05a3.pdf", width = 4, height = 3.5)

  par(mar = c(4, 2, 2, 0.5) + 0.1)
  
  plot(x = freq, xlim = c(0, 1), ylim = c(lim4, lim2), type="l",
      xlab = expression(tau), ylab = "")
  
  for (i in 1:3) {   
    qcoh <- rep(0, length(quantile))
    for (j in 1:length(quantile)) {
      qcoh[j] <- Re(V[i, 1, j, 2, j]/sqrt(V[i, 1, j, 1, j]*V[i, 2, j, 2, j]))
    }
    lines(x = quantile, y = qcoh, col = "black", lty = i, lwd = tht)
  }
  legend("bottom", inset = .03, "center", c("W", "M", "Y"),
         cex = 0.65, lwd = c(1, 1, 1),
         lty = c(1:3), horiz = TRUE, bg = "white")

dev.off()

###############################################################################
## VAR(1) Gumbel
###############################################################################

load("VAR1_Gumbel.Rdata")

# (ii) Quantile coherency (for 4 different combinations of parameters
#                              and all frequencies)

V <- getValues(qSD, levels.1 = c(0.5, 0.05, 0.95),
                    levels.2 = c(0.5, 0.05, 0.95))
V <- V[1:(dim(V)[1]/2),,,,]

frequencies <- (0:le)/(le+1)
freq <- frequencies[which(frequencies != 0)][1:(le/2)]

pdf("fig_05b1.pdf", width = 4, height = 3.5)
  par(mar = c(4, 2, 2, 0.5) + 0.1)
  
  plot(x = freq, xlim = c(0, 0.5), ylim = c(lim1, lim2), type = "l",
       xlab = expression(omega / 2*pi), ylab = "")
  for (i1 in 1:3) {
    lines(x = freq, y = Re(V[,1,i1,2,i1]/sqrt(V[,1,i1,1,i1]*V[,2,i1,2,i1])),
          col = "black", lty = i1, lwd = tht)
  }
  axis(side = 3, at = c(1/5, 1/22, 1/250), labels = c("W", "M", "Y"))
  abline(v = c(1/5, 1/22, 1/250), col = "gray")
  legend("bottom", inset = .03, "center",
         c("0.5 | 0.5", "0.05 | 0.05", "0.95 | 0.95"),
         cex = 0.65, lwd = c(1, 1, 1), lty = c(1:3),
         horiz = TRUE, bg = "white")

dev.off()

pdf("fig_05b2.pdf", width = 4, height = 3.5)
  par(mar = c(4, 2, 2, 0.5) + 0.1)
  plot(x = freq, xlim = c(0, 0.5), ylim = c(-0.2, lim2), type="l",
       xlab = expression(omega / 2 * pi), ylab = "")
  lines(x = freq, y = Re(V[,1,1,2,3]/sqrt(V[,1,1,1,1]*V[,2,3,2,3])),
        col = "black", lty = i1, lwd = tht) 
  axis(side = 3, at = c(1/5, 1/22, 1/250), labels = c("W", "M", "Y"))
  abline(v = c(1/5, 1/22, 1/250), col = "gray")
  legend("bottom", inset = .03, "center", c("0.05 | 0.95"),
         cex = 0.8, lwd = c(1, 1, 1), lty = c(1:3),
         horiz = TRUE, bg = "white")
dev.off()

# (iii) Quantile coherency (for 3 frequencies and all taus)

V <- getValues(qSD, frequencies = 2*pi * c(1/5, 1/22, 1/250))

pdf("fig_05b3.pdf", width = 4, height = 3.5)

par(mar = c(4, 2, 2, 0.5) + 0.1)

  plot(x = freq, xlim = c(0, 1), ylim = c(lim4, lim2), type = "l",
      xlab = expression(tau), ylab = "")
  
  for (i in 1:3) {   
    qcoh <- rep(0, length(quantile))
    for (j in 1:length(quantile)) {
      qcoh[j] <- Re(V[i, 1, j, 2, j]/sqrt(V[i, 1, j, 1, j]*V[i, 2, j, 2, j]))
    }
    lines(x = quantile, y = qcoh, col = "black", lty = i, lwd = tht)
  }
  legend("bottom", inset = .03, "center", c("W", "M", "Y"),
         cex = 0.65, lwd = c(1, 1, 1),
         lty = c(1:3), horiz = TRUE, bg = "white")

dev.off()


###############################################################################
## VAR(1) Clayton
###############################################################################

load("VAR1_Clayton.Rdata")

# (ii) Quantile coherency (for 4 different combinations of parameters
#                              and all frequencies)

V <- getValues(qSD, levels.1 = c(0.5, 0.05, 0.95),
                    levels.2 = c(0.5, 0.05, 0.95))
V <- V[1:(dim(V)[1]/2),,,,]

frequencies <- (0:le)/(le+1)
freq <- frequencies[which(frequencies != 0)][1:(le/2)]

pdf("fig_05c1.pdf", width = 4, height = 3.5)
  par(mar = c(4, 2, 2, 0.5) + 0.1)
  
  plot(x = freq, xlim = c(0, 0.5), ylim = c(lim1, lim2), type = "l",
       xlab = expression(omega / 2*pi), ylab = "")
  for (i1 in 1:3) {
    lines(x = freq, y = Re(V[,1,i1,2,i1]/sqrt(V[,1,i1,1,i1]*V[,2,i1,2,i1])),
          col = "black", lty = i1, lwd = tht)
  }
  axis(side = 3, at = c(1/5, 1/22, 1/250), labels = c("W", "M", "Y"))
  abline(v = c(1/5, 1/22, 1/250), col = "gray")
  legend("bottom", inset = .03, "center",
         c("0.5 | 0.5", "0.05 | 0.05", "0.95 | 0.95"),
         cex = 0.65, lwd = c(1, 1, 1), lty = c(1:3),
         horiz = TRUE, bg = "white")

dev.off()

pdf("fig_05c2.pdf", width = 4, height = 3.5)
  par(mar = c(4, 2, 2, 0.5) + 0.1)
  plot(x = freq, xlim = c(0, 0.5), ylim = c(-0.2, lim2), type = "l",
       xlab = expression(omega / 2 * pi), ylab = "")
  lines(x = freq, y = Re(V[,1,1,2,3]/sqrt(V[,1,1,1,1]*V[,2,3,2,3])),
        col = "black", lty = i1, lwd = tht)
  axis(side = 3, at = c(1/5, 1/22, 1/250), labels = c("W", "M", "Y"))
  abline(v = c(1/5, 1/22, 1/250), col = "gray")
  legend("bottom", inset = .03, "center", c("0.05 | 0.95"),
         cex = 0.8, lwd =c (1, 1, 1), lty = c(1:3),
         horiz = TRUE, bg = "white")
dev.off()

# (iii) Quantile coherency (for 3 frequencies and all taus)

V <- getValues(qSD, frequencies = 2 * pi * c(1/5, 1/22, 1/250))

pdf("fig_05c3.pdf", width = 4, height = 3.5)

  par(mar = c(4, 2, 2, 0.5) + 0.1)
  
  plot(x = freq, xlim = c(0, 1), ylim = c(lim4, lim2), type = "l",
      xlab = expression(tau), ylab = "")
  
  for (i in 1:3) {   
    qcoh <- rep(0, length(quantile))
    for (j in 1:length(quantile)) {
      qcoh[j] <- Re(V[i, 1, j, 2, j]/sqrt(V[i, 1, j, 1, j]*V[i, 2, j, 2, j]))
    }
    lines(x = quantile, y = qcoh, col = "black", lty = i, lwd = tht)
  }
  legend("bottom", inset = .03, "center", c("W", "M", "Y"),
         cex = 0.65, lwd = c(1, 1, 1),
         lty = c(1:3), horiz = TRUE, bg="white")

dev.off()
