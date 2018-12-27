
###################################################################
## Code to replicate the results in the Barunik and Kley (2018)
###################################################################
##
## Code to the Supplementary Figure 4
##
## Authors: Jozef Barunik, and Tobias Kley (2018)
###################################################################

# NOTE: this code produces a pdf figure that is saved to the current 
#       working folder


# load the auxiliary files
source("quantile_coherency_replication_pack.R")

# set seed for RNG
set.seed(1234)

# Set length of the generated time series
le <- 2048

# Choose quantiles for which the quantities are computed
quantile <- c(0.5,0.05,0.95)
quantilelegend <- c("0.05 | 0.05",
                     "0.5 | 0.5",
                    "0.05 | 0.95",
                     "0.5 | 0.05")

# simulate the quantile cross-spectral densities
qSD <- quantileSD(le, ts = indepR, levels.1 = quantile, R = 1000)
V <- getValues(qSD)
V <- V[1:(dim(V)[1]/2),,,,]

# FIGURE S.4 (a)
pdf(file = paste("fig_S04a.pdf", sep = ''), width = 9, height = 3)
  par(mfrow = c(1, 3))
  par(mar = c(4, 4, 1, 0.5) + 0.1)
    
  plot(x = c(0, 0.5), y = c(0, 0), type = "l", ylim = c(-1, 1), lty = 3,
       ylab = "traditional coherency (Re)", xlab = expression(omega / 2 * pi))
    
  for (i in seq(-.9, .9, .1)) {
    lines(x = c(0, 0.5), y = c(i, i), lty = 3)  
  }
  lines(x = c(0, 0.5), y = c(0.6, 0.6), lty = 1, lwd = 2)
  
  plot(transl(seq(-1, 1, 0.01), 0.05, 0.05), type = "l",
       xaxt = "n", ylim = c(-1, 1), lty = 2,
       ylab = "quantile coherency (Re)", xlab = "traditional coherency (Re)")
  lines(transl(seq(-1, 1, 0.01), 0.5, 0.5), col = "black", lty = 1)
  lines(transl(seq(-1, 1, 0.01), 0.05, 0.95), col = "black", lty = 3)
  lines(transl(seq(-1, 1, 0.01), 0.5, 0.05), col = "black", lty = 4)
  axis(side = 1, at = c(1, 50, 100, 150, 200), labels = c(-1, -0.5, 0, 0.5, 1))
  abline(h = c(transl(0.6, 0.05, 0.05)), v = c(161), col="gray")
  abline(h = c(transl(0.6, 0.5, 0.5)), col = "gray")
  abline(h = c(transl(0.6, 0.05, 0.95)), col = "gray")
  abline(h = c(transl(0.6, 0.5, 0.05)), col = "gray")
  legend("topleft", inset = .03, quantilelegend,
         cex = 0.73, lty = c(2, 1, 3, 4), ncol = 2)
    
  # Plot Quantile Coherency
  frequencies  <-  (0:le)/(le+1)
  freq <- frequencies[which(frequencies != 0)][1:(le/2)]
    
  plot(x = freq, xlim = c(0, 0.5), ylim = c(-1, 1), type = "l",
       xlab = expression(omega / 2 * pi), ylab = "")
  lines(x = freq, y = Re(V[,1,2,2,2] / sqrt(V[,1,2,1,2] * V[,2,2,2,2])),
        col = "black", lty = 2)
  lines(x = freq, y = Re(V[,1,1,2,1] / sqrt(V[,1,1,1,1] * V[,2,1,2,1])),
        col = "black", lty = 1)
  lines(x = freq, y = Re(V[,1,2,2,3] / sqrt(V[,1,2,1,2] * V[,2,3,2,3])),
        col = "black", lty = 3)
  lines(x = freq, y = Re(V[,1,1,2,2] / sqrt(V[,1,1,1,1] * V[,2,2,2,2])),
        col = "black", lty = 4)
  legend("bottom", inset = .03, quantilelegend,
         cex = 0.73, lty = c(2, 1, 3, 4), ncol = 2)

dev.off()



# FIGURE S.4 (b)
#
# (1) Determine a set of VAR(1) models such that the
# 			process is stable.
# (2) Determine the corresponding traditional and quantile
# 			spectral densities and coherencies.

M <- matrix(0, ncol = 4)

S1 <- c(seq(-1.95, 1.95, 0.05))
S2 <- c(seq(-1.95, 1.95, 0.05))

for (x1 in c(0)) {
  for (x2 in S2) {
    for (x3 in c(x2)) {
      for (x4 in c(x1)) {
        
        is_stable <- FALSE
        
        # if quadratic
        detA <- x1 * x4 - x2 * x3
        if ( detA != 0) {
          D <- (x1 - x4)^2 + 4 * x2 * x3
          if ( D < 0 ) {
            is_stable <- TRUE
          } else {
            z01 <- 0.5 * (x1 + x4 + sqrt(D)) / detA
            z02 <- 0.5 * (x1 + x4 - sqrt(D)) / detA
            if ( (abs(round(z01, 8)) > 1) && (abs(round(z02, 8)) > 1) ) {
              is_stable <- TRUE
            }  
          }
          
          # if linear
        } else {
          if ( x1 + x4 == 0 ) {
            is_stable <- TRUE
          } else {
            if ( 1 / abs(round(x1 + x4, 8)) > 1 ) {
              is_stable <- TRUE
            }
          }
        }
        
        if (is_stable) {
          M <- rbind(M, round(c(x1, x2, x3, x4), 8))
        }       
      }
    }
  }
}
M <- M[-1,]


Sigma <- diag(c(1, 1))
levels <- c(0.05, 0.5, 0.95)

maxLag <- 511

D <- 2 # TODO
K <- length(levels)

tradCoh <- array(0, dim = c(nrow(M), maxLag + 1, D, D))
quantCoh <- array(0, dim = c(nrow(M), maxLag + 1, D, K, D, K))

for (i in 1:nrow(M)) {
  A <- matrix(M[i, 1:4], ncol = 2, byrow = TRUE)
  tqCoh <- getTradQuantCoherencyInVAR1(A, Sigma, levels, maxLag = maxLag)
  tradCoh[i,,,] <- tqCoh$tradCoh
  quantCoh[i,,,,,] <- tqCoh$quantCoh
  cat(i," ", nrow(M), "\n")
}

quantile <- c(0.5, 0.05, 0.95)
quantilelegend <- c("0.05 | 0.05",
                     "0.5 | 0.5",
                    "0.05 | 0.95",
                     "0.5 | 0.05")

pdf(file = paste("fig_S04b.pdf", sep = ""), width = 9, height = 3)
  par(mfrow = c(1, 3))
  par(mar = c(4, 4, 1, 0.5) + 0.1)
  
  
  plot(x = 0:256/512, y = tradCoh[1, 1:257, 1, 2], type = "l",
       ylim = c(-1, 1), lty = 3,
       ylab = "traditional coherency (Re)", xlab = expression(omega / 2 * pi))
  
  for (i in seq(2, 39, 1)) {
    lines(x = 0:256/512, y = tradCoh[i, 1:257, 1, 2], lty = 3)  
  }
  lines(x = 0:256/512, y = tradCoh[29, 1:257, 1, 2], lty = 1, lwd = 2)
  
  abline(h = Re(tradCoh[29, 53, 1, 2]), col = "gray")
  abline(v = 52 / 512, col = "gray")
  
  freq <- 53
  i1 <- 1
  i2 <- 2
  k1 <- 1
  k2 <- 1
  plot(x = Re(tradCoh[, freq, i1, i2]), y = Re(quantCoh[,freq,i1,k1,i2,k2]),
       type = "l", xlim = c(-1, 1), ylim = c(-1, 1),
       lty = 2, ylab = "quantile coherency (Re)",
       xlab = "traditional coherency (Re)",
       main = expression("frequency =" ~ 2 * pi * 52 / 512))
  abline(h = Re(quantCoh[29, 53, 1, k1, 2, k2]), col = "gray")
  
  k1 <- 2
  k2 <- 2
  lines(x = Re(tradCoh[, freq, i1, i2]),
        y = Re(quantCoh[, freq, i1, k1, i2, k2]), col = "black", lty = 1)
  abline(h = Re(quantCoh[29, 53, 1, k1, 2, k2]), col = "gray")
  
  k1 <- 1
  k2 <- 3
  lines(x = Re(tradCoh[, freq, i1, i2]),
        y = Re(quantCoh[, freq, i1, k1, i2, k2]), col = "black", lty = 3)
  abline(h = Re(quantCoh[29, 53, 1, k1, 2, k2]), col = "gray")
  
  k1 <- 2
  k2 <- 1
  lines(x = Re(tradCoh[, freq, i1, i2]),
        y=Re(quantCoh[, freq, i1, k1, i2, k2]), col = "black", lty = 4)
  abline(h = Re(quantCoh[29, 53, 1, k1, 2, k2]), col = "gray")
  
  abline(v = Re(tradCoh[29, 53, 1, 2]), col = "gray")
  
  legend("topleft", inset = .03, quantilelegend,
         cex = 0.73, lty = c(2, 1, 3, 4), ncol = 2)
  
  # Plot Quantile Coherency
  frequencies <- (0:le) / (le + 1)
  freq <- frequencies[which(frequencies != 0)][1:(le/2)]
  
  plot(x = 0:256/512, xlim = c(0, 0.5), ylim = c(-1, 1), type="l",
       xlab = expression(omega / 2 * pi), ylab = "")
  abline(v = 52 / 512, col = "gray")
  
  freq <- 53
  k1 <- 1
  k2 <- 1
  lines(x = 0:256/512, y = Re(quantCoh[29, 1:257, i1, k1, i2, k2]),
        col = "black", lty = 2)
  abline(h = Re(quantCoh[29, 53, 1, k1, 2, k2]), col = "gray")
  
  k1 <- 2
  k2 <- 2
  lines(x = 0:256/512, y = Re(quantCoh[29, 1:257, i1, k1, i2, k2]),
        col = "black", lty = 1)
  abline(h = Re(quantCoh[29, 53, 1, k1, 2, k2]), col = "gray")
  
  k1 <- 1
  k2 <- 3
  lines(x = 0:256/512, y = Re(quantCoh[29, 1:257, i1, k1, i2, k2]),
        col = "black", lty = 3)
  abline(h = Re(quantCoh[29, 53, 1, k1, 2, k2]), col = "gray")
  
  k1 <- 2
  k2 <- 1
  lines(x = 0:256/512, y = Re(quantCoh[29, 1:257, i1, k1, i2, k2]),
        col = "black", lty = 4)
  abline(h = Re(quantCoh[29, 53, 1, k1, 2, k2]), col = "gray")
  
  legend("bottomright", inset = .03, quantilelegend,
         cex = 0.73, lty = c(2, 1, 3, 4), ncol = 2)

dev.off()



# FIGURE S.4 (c)

M <- matrix(0, ncol=4)

S1 <- c(seq(-.99, -.95, 0.01),
        seq(-.9, -.05, 0.05),
        seq(-.04, 0.04, 0.01),
        seq(.05, 0.9, 0.05),
        seq(.95, .99, 0.01))
S2 <- S1


for (x1 in S1) {
  for (x2 in S2) {
    for (x3 in c(x2)) {
      for (x4 in c(x1)) {
        
        is_stable <- FALSE
        
        # if quadratic
        detA <- x1 * x4 - x2 * x3
        if ( detA != 0) {
          D <- (x1 - x4)^2 + 4 * x2 * x3
          if ( D < 0 ) {
            is_stable <- TRUE
          } else {
            z01 <- 0.5 * (x1 + x4 + sqrt(D)) / detA
            z02 <- 0.5 * (x1 + x4 - sqrt(D)) / detA
            if ( (abs(round(z01, 8)) > 1) && (abs(round(z02, 8)) > 1) ) {
              is_stable <- TRUE
            }  
          }
          
          # if linear
        } else {
          if ( x1 + x4 == 0 ) {
            is_stable <- TRUE
          } else {
            if ( 1 / abs(round(x1 + x4, 8)) > 1 ) {
              is_stable <- TRUE
            }
          }
        }
        
        if (is_stable) {
          M <- rbind(M, round(c(x1, x2, x3, x4), 8))
        }       
      }
    }
  }
}
M <- M[-1,]


Sigma <- diag(c(1, 1))
levels <- c(0.05, 0.5, 0.95) 

maxLag <- 1023

D <- 2
K <- length(levels)

tradCoh <- array(0, dim = c(nrow(M), maxLag + 1, D, D))
quantCoh <- array(0, dim = c(nrow(M), maxLag + 1, D, K, D, K))

for (i in 1:nrow(M)) {
  A <- matrix(M[i, 1:4], ncol = 2, byrow = TRUE)
  tqCoh <- getTradQuantCoherencyInVAR1(A, Sigma, levels, maxLag = maxLag)
  tradCoh[i,,,] <- tqCoh$tradCoh
  quantCoh[i,,,,,] <- tqCoh$quantCoh
  cat(i, " ", nrow(M), "\n")
}



quantile <- c(0.5, 0.05, 0.95)
quantilelegend <- c("0.05 | 0.05",
                    "0.05 | 0.95")

pdf(file=paste("fig_S04c.pdf", sep = ""), width = 10.5, height = 3.5)
  
  par(mfrow = c(1, 3))
  par(mar = c(4, 4, 1, 0.5) + 0.1)
  
  plot(x = c(0, 0.5), y = c(0,0), type = "l",
       ylim = c(-1, 1), lty = 3,
       ylab = "traditional coherency (Re)", xlab = expression(omega / 2 * pi))
  
  for (i in seq(1, 761, 20)) {
    lines(x = 0:512/1024, y = Re(tradCoh[i, 1:513, 1, 2]), lty = 3)  
  }
  
  miB <- c(287, 614, 1221)
  
  for(i in miB) {
    lines(x = 0:512/1024, y = Re(tradCoh[i, 1:513, 1, 2]), lty = 1, lwd = 2)  
  }
  
  abline(h = mean(Re(tradCoh[miB, 105, 1, 2])), col = "gray")
  abline(v = 104 / 1024, col = "gray")
  
  
  ## TRANSLATION PLOT
  
  freq <- 105
  i1 <- 1
  i2 <- 2
  
  ## end: determine points on the alpha convex hull
  
  plot(x = c(-1, 1), y = c(0, 0),
       type = "n", xlim = c(-1, 1), ylim = c(-1, 1),
       ylab = "quantile coherency (Re)", xlab = "traditional coherency (Re)",
       main = expression("frequency =" ~ 2 * pi * 52 / 512),
       pch = 16, cex = .5, col = "darkgray")
  
  
  k1 <- 1
  k2 <- 1
  X <- matrix(c(Re(tradCoh[, freq, i1, i2]),
                Re(quantCoh[, freq, i1, k1, i2, k2])), ncol = 2)
  
  aCH <- aConvHull(X, a = 0.1)
  polygon(aCH, density = 20, angle = -45, lty = 2)
  
  abline(v = mean(Re(tradCoh[miB,freq,1,2])), col = "gray")
  abline(h = Re(quantCoh[miB,freq,1,k1,2,k2]), col = "gray")
  
  k1 <- 1
  k2 <- 3
  
  X <- matrix(c(Re(tradCoh[, freq, i1, i2]),
                Re(quantCoh[, freq, i1, k1, i2, k2])), ncol = 2)
  aCH <- aConvHull(X)
  polygon(aCH, density = 20, angle = 45, lty = 3)
  
  abline(v = mean(Re(tradCoh[miB, freq, 1, 2])), col = "gray")
  abline(h = Re(quantCoh[miB, freq, 1, k1, 2, k2]), col = "gray")
  
  legend("topleft", inset = .03, quantilelegend,
         cex = 0.73, lty = c(2, 3), ncol = 2)
  
  # COHERENCY AT QUANTILES
  
  plot(x = 0:512/1024, xlim = c(0, 0.5), ylim = c(-1, 1), type = "l",
       xlab = expression(omega / 2 * pi), ylab = "")
  abline(v = 105 / 1024, col = "gray")
  
  freq <- 105
  
  for (i in miB) {
    
    k1 <- 1
    k2 <- 1
    lines(x = 0:512/1024, y = Re(quantCoh[i, 1:513, 1, k1, 2, k2]),
          col = "black", lty = 2)
    abline(h = Re(quantCoh[i, freq, 1, k1, 2, k2]), col = "gray")
    
    k1 <- 1
    k2 <- 3
    lines(x = 0:512/1024, y = Re(quantCoh[i, 1:513, 1, k1, 2, k2]),
          col = "black", lty = 3)
    abline(h = Re(quantCoh[i, freq, 1, k1, 2, k2]), col = "gray")
  }
  
  
  legend("bottomright", inset = .03, quantilelegend,
         cex = 0.73, lty = c(2, 3), ncol = 2)

dev.off()

