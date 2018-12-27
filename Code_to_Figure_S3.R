
###################################################################
## Code to replicate the results in the Barunik and Kley (2018)
###################################################################
##
## Code to the Supplementary Figure 3
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
quantile <- c(0.5, 0.05, 0.25, 0.95) 
quantilelegend <- c("0.05 | 0.05",
                    "0.95 | 0.95",
                    "0.25 | 0.25",
                     "0.5 | 0.5",
                     "0.5 | 0.95")

# FIGURE S.3
# ==========

# choose Data Generating Process
model <- quar33

# chose name of the output file
name <- "fig_S03.pdf"

# simulate the quantile cross-spectral densities
qSD <- quantileSD(le, ts = model, levels.1 = quantile, R = 1000)
V <- getValues(qSD)
V <- V[1:(dim(V)[1]/2),,,,]

pdf(file = name, width = 8, height = 6)
  par(mar = c(4, 2, 2, 0.5) + 0.1)
  par(mfrow = c(2, 3))
  
  frequencies <- (0:le) / (le+1)
  freq <- frequencies[which(frequencies != 0)][1:(le/2)]

  # Compute Traditional Coherency
  ccR <- matrix(0, 1000, le / 2)
  ccI <- matrix(0, 1000, le / 2)
  for(i in 1:1000){
    dog <- mvspec(model(le), spans = c(100, 100), plot = FALSE)
    ccR[i,] <- Re(dog$fxx[1,2,] / (sqrt(dog$fxx[1,1,]) * sqrt(dog$fxx[2,2,])))
    ccI[i,] <- Im(dog$fxx[1,2,] / (sqrt(dog$fxx[1,1,]) * sqrt(dog$fxx[2,2,])))
  }
  plot(x = freq, y = apply(ccR, 2, mean), xlim = c(0, 0.5), ylim = c(-1, 1),
       type = "l", xlab = expression(omega / 2 * pi), ylab = "")
  
  # Plot Quantile Coherency (Real part)
  plot(x = freq, xlim = c(0, 0.5), ylim = c(-1, 1), type = "l",
       xlab = expression(omega / 2 * pi), ylab = "")
  for (i1 in 1:length(quantile)) {
    lines(x = freq, y = Re(V[,1,i1,2,i1] / sqrt(V[,1,i1,1,i1] * V[,2,i1,2,i1])),
          col = "black", lty = i1)
  }
  i2 <- 1
  i3 <- 4
  lines(x = freq, y = Re(V[,1,i2,2,i3] / sqrt(V[,1,i2,1,i2] * V[,2,i3,2,i3])),
        lty = 5)
  legend("bottom", inset = .03, "center", quantilelegend,
         cex = 0.73, lty = c(2, 4, 3, 1, 5), ncol = 3)

  # Plot Quantile Co-spectrum (Real part)
  plot(x = freq, xlim = c(0, 0.5), ylim = c(-0.05, 0.05), type = "l",
       xlab = expression(omega / 2 * pi), ylab = "")
  for (i1 in 1:length(quantile)) {
    lines(x = freq, y = Re(V[,1,i1,2,i1]), lty = i1)
  }
  i2 <- 1
  i3 <- 4
  lines(x = freq, y = Re(V[,1,i2,2,i3]), lty = 5)
  legend("bottom", inset = .03, "center", quantilelegend,
         cex = 0.73, lty = c(2, 4, 3, 1, 5), ncol = 3)


  # Traditional Coherency (Imaginary part)
  plot(x = freq, y = apply(ccI, 2, mean), xlim = c(0, 0.5), ylim = c(-1, 1),
       type = "l", xlab = expression(omega / 2 * pi), ylab = "")
  
  # Plot Quantile Coherency (Imaginary part)
  plot(x = freq, xlim = c(0, 0.5), ylim = c(-1, 1), type = "l",
       xlab = expression(omega / 2 * pi), ylab = "")
  for (i1 in 1:length(quantile)) {
    lines(x = freq, y = Im(V[,1,i1,2,i1] / sqrt(V[,1,i1,1,i1] * V[,2,i1,2,i1])),
          col = "black", lty = i1)
  }
  i2 <- 1
  i3 <- 4
  lines(x = freq, y = Im(V[,1,i2,2,i3] / sqrt(V[,1,i2,1,i2] * V[,2,i3,2,i3])),
        lty = 5)
  legend("bottom", inset = .03, "center", quantilelegend,
         cex = 0.73, lty = c(2, 4, 3, 1, 5), ncol = 3)

  # Plot Quantile Co-spectrum (Imaginary part)
  plot(x = freq, xlim = c(0, 0.5), ylim = c(-0.05, 0.05), type = "l",
       xlab = expression(omega / 2 * pi), ylab = "")
  for (i1 in 1:length(quantile)) {
    lines(x = freq, y = -Im(V[,1,i1,2,i1]), lty = i1)
  }
  i2 <- 1
  i3 <- 4
  lines(x = freq, y = Im(V[,1,i2,2,i3]), lty = 5)
  legend("bottom", inset = .03, "center", quantilelegend,
         cex = 0.73, lty = c(2, 4, 3, 1, 5), ncol = 3)

dev.off()
