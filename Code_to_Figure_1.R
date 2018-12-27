
###################################################################
## Code to replicate the results in the Barunik and Kley (2018)
###################################################################
##
## Code to Figure 1
##
## Authors: Jozef Barunik, and Tobias Kley (2018)
###################################################################

# NOTE: this code produces three pdf figures fig_01a.pdf,
#       fig_01b.pdf, and fig_01c.pdf that are saved to the current 
#       working folder


# load the auxiliary files
source("quantile_coherency_replication_pack.R")

# set seed for RNG
set.seed(1234)

# Set length of the generated time series
le <- 2048

# Choose quantile levels for which the quantities are computed
quantile <- c(0.5, 0.05, 0.25, 0.75, 0.95)

# legend for the plots
quantilelegend <- c("0.05 | 0.05",
                    "0.95 | 0.95",
                    "0.25 | 0.25",
                    "0.75 | 0.75",
                     "0.5 | 0.5")

# FIGURE 1 (a)
# ============

# choose Data Generating Process
model <- indepq

# chose name of the output file
name <- "fig_01a.pdf"

# simulate the quantile cross-spectral densities
qSD <- quantileSD(le, ts = model, levels.1 = quantile, R = 1000)
V <- getValues(qSD)
V <- V[1:(dim(V)[1]/2),,,,]

pdf(file = name, width = 3, height = 3.5)
  par(mar = c(4, 2, 2, 0.5) + 0.1)
  
  # Plot Quantile Coherency
  frequencies <- (0:le)/(le+1)
  freq <- frequencies[which(frequencies != 0)][1:(le/2)]
  
  plot(x = freq, xlim = c(0, 0.5), ylim = c(-1, 1), type = "l",
       xlab = expression(omega / 2 * pi), ylab = "")
  for (i1 in 1:length(quantile)) {
    lines(x = freq, y = Re(V[,1,i1,2,i1] / sqrt(V[,1,i1,1,i1] * V[,2,i1,2,i1])),
          col = "black", lty = i1, lwd = 1.5)
  }
  legend("bottom", inset = .03, "center", quantilelegend, cex = 0.53,
         lty = c(2, 5, 3, 4, 1), ncol = 3)
dev.off()

# FIGURE 1 (b)
# ============

# choose Data Generating Process
model <- indepqshift

# chose name of the output file
name <- "fig_01b.pdf"

# simulate the quantile cross-spectral densities
qSD <- quantileSD(le, ts = model, levels.1 = quantile, R = 1000) 
V <- getValues(qSD)
V <- V[1:(dim(V)[1]/2),,,,]

pdf(file = name, width = 3, height = 3.5)
  par(mar = c(4, 2, 2, 0.5) + 0.1)
  
  # Plot Quantile Coherency
  frequencies <- (0:le)/(le+1)
  freq <- frequencies[which(frequencies != 0)][1:(le/2)]
  
  plot(x = freq, xlim = c(0, 0.5), ylim = c(-1, 1), type = "l",
       xlab = expression(omega / 2 * pi), ylab = "")
  for (i1 in 1:length(quantile)) {
    lines(x = freq, y = Re(V[,1,i1,2,i1] / sqrt(V[,1,i1,1,i1] * V[,2,i1,2,i1])),
          col = "black", lty = i1, lwd = 1.5)
  }
  legend("bottom", inset = .03, "center", quantilelegend, cex = 0.53,
         lty=c(2, 5, 3, 4, 1), ncol = 3)
dev.off()

# FIGURE 1 (c)
# ============

# choose Data Generating Process
model <- indep

# chose name of the output file
name <- "fig_01c.pdf"

# simulate the quantile cross-spectral densities
qSD <- quantileSD(le, ts = model, levels.1 = quantile, R = 1000)
V <- getValues(qSD)
V <- V[1:(dim(V)[1]/2),,,,]

pdf(file = name, width = 3, height = 3.5)
  par(mar = c(4, 2, 2, 0.5) + 0.1)
  
  # Plot Quantile Coherency
  frequencies <- (0:le)/(le+1)
  freq <- frequencies[which(frequencies != 0)][1:(le/2)]
  
  plot(x = freq, xlim = c(0, 0.5), ylim = c(-1, 1), type = "l",
       xlab = expression(omega / 2 * pi), ylab = "")
  for (i1 in 1:length(quantile)) {
    lines(x = freq, y = Re(V[,1,i1,2,i1]/sqrt(V[,1,i1,1,i1]*V[,2,i1,2,i1])),
          col = "black",lty = i1, lwd = 1.5)
  }
  legend("bottom", inset = .03, "center", quantilelegend, cex = 0.53,
         lty = c(2, 5, 3, 4, 1), ncol = 3)
dev.off()
