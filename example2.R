# TVAR(6) process 
# 
#
#   Model 12 in
#   Rosen, O., Stoffer, D., and Wood, S, (2009)
#     ``Local Spectral Analysis via a Bayesian Mixture of Smoothing Splines,''   
#       JASA,104(485), 249-262
#
rm(list = ls())

##### load library
library("fields")
#library("signal")

##### load R functions of Bayesian lattice filters
source("./Rfun/BayesLatticeLik.R")
source("./Rfun/BayesLattice.R")
source("./Rfun/dynamicLik.R")
source("./Rfun/tvar_spec.R")


##### Simulate a TVAR(6) process
N <- 1024
et <- rnorm(N)
xt <- rnorm(6)
true_coeffs <- matrix(0, nrow=6, ncol=N)

signal <- numeric(N)

for (tt in 1:N) {
  	theta1 <- 0.05 + tt*(0.1/(N-1))
		theta2 <- 0.25
		theta3 <- 0.45 - tt*(0.1/(N-1))
    
		a1_r <- 1/1.1*exp(2i*pi*theta1)
		a1_c <- Conj(a1_r)
	
		a2_r <- 1/1.12*exp(2i*pi*theta2)
		a2_c <- Conj(a2_r)
	
		a3_r <- 1/1.1*exp(2i*pi*theta3)
		a3_c <- Conj(a3_r)
	
		ar_root <- c(a1_r, a1_c, a2_r, a2_c, a3_r, a3_c)
		coeff <- Re(poly(ar_root))
		coeff <- -coeff[2:7]
		true_coeffs[,tt] <- coeff
	
		signal[tt] <- xt[1]*coeff[1] + xt[2]*coeff[2]+ xt[3]*coeff[3]+ xt[4]*coeff[4]+ xt[5]*coeff[5] + xt[6]*coeff[6] + et[tt]
		xt[2:6] <- xt[1:5]
		xt[1] <- signal[tt]
}

plot(1:N, signal, main="TVAR(6)", type="l", xlab="Time", ylab="")

##### Procedure 1: Search orders of TVAR models #####
disSys <- seq(0.8,1, by = 0.02)
disMea <- seq(0.8,1, by = 0.02)
P <- 10

para_combin <- BayesLatticeLik(signal, P, disSys, disMea)

# BLFscree plot
plot(1:P,para_combin[,4], type = "l",
     xlab = "Order", 
		 ylab = "log(likelihood)", 
		 main = "")
   

##### Procedure 2: Obtain time-varying coefficients, innovation variance, and PARCOR of TVAR(P) #####
sel_order <- 6 ## Remember to change it according to the result of order selection
Dfactor <- para_combin[1:sel_order,2:3]

tvarp <- BayesLattice(signal, Dfactor)

par(mfrow = c(3,1))
plot(1:N, tvarp$parcor[1,], type="l", main="PARCOR coefficients", xlab="", ylab="", col=1, ylim=c(-1,1))
for (ii in 2:6) lines(1:N, tvarp$parcor[ii,], col=ii)
abline(h=0, col="gray")

ylim <- range(tvarp$coefs)
plot(1:N, tvarp$coefs[1,], type="l", main="TVAR coefficients", xlab="", ylab="", col=1, ylim=ylim)
for (ii in 2:6) lines(1:N, tvarp$coefs[ii,], col=ii)

plot(1:N, tvarp$s2, type="l", main="Variance", xlab="", ylab="", col=1)


##### Procedure 3: Make spectrum plots #####
times <- 1:N
freqs <- seq(0, 0.5, by = 0.005)
spec <- tvar_spec(tvarp$coefs, tvarp$s2, times, freqs)
spec <- log(spec)

tspec <- tvar_spec(true_coeffs, rep(1,N), times, freqs)
tspec <- log(tspec)

set.panel() # reset plotting device
par(oma=c( 0,0,0,6)) # margin of 4 spaces width at right hand side
set.panel( 3,1) 

plot(times, signal, type ="l", xlab = "", ylab = "", main = "TVAR(6)")

image(times, freqs, t(tspec), main="True Spectrogram", xlab = "", ylab = "Frequency", col = tim.colors(64))
image.plot(tspec, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)

image(times, freqs, t(spec), main="Posterior mean", xlab = "", ylab = "Frequency", col = tim.colors(64))
image.plot(spec, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)

set.panel() 


