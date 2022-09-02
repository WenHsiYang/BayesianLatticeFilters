# TVAR(2) process 
# 
#            
#    X_t = a_t * X_{t-1} - 0.81X_{t-2} + e_t, for  1<= t <= 1024
#            
#        with a_t = 0.8(1-0.5cos(pi*t/1024)) and e_t ~ N(0,1).
#
#   Model 10 in
#   Rosen, O., Wood, S, and Stoffer, D., (2012)
#     ``AdaptSPEC: Adaptive Spectral Estimation for Nonstationary Time Series,''   
#       JASA,107(500), 1575-1589
#
#   Model 9 in
#   Rosen, O., Stoffer, D., and Wood, S, (2009)
#     ``Local Spectral Analysis via a Bayesian Mixture of Smoothing Splines,''   
#       JASA,104(485), 249-262
#
rm(list = ls())

##### load library
library("fields") # for plotting image plots

##### load R functions of Bayesian lattice filters
source("./Rfun/BayesLatticeLik.r")
source("./Rfun/BayesLattice.r")
source("./Rfun/dynamicLik.r")
source("./Rfun/tvar_spec.r")

##### Simulate a TVAR(2) process
N <- 1024
et <- rnorm(N)
xt1 <- rnorm(1)
xt2 <- rnorm(1)

true_coeffs <- matrix(NA, 2, N)
signal <- numeric(N)

for (tt in 1:N) {
  at <- 0.8*(1-0.5*cos(pi*tt/N))
  signal[tt] <- at*xt1 - 0.81*xt2 + et[tt]
	
	true_coeffs[1,tt] <- at
	true_coeffs[2,tt] <- -0.81
	
	xt2 <- xt1
  xt1 <- signal[tt]
}

plot(1:N, signal, main="TVAR(2)", type="l", xlab="Time", ylab="")

##### Procedure 1: Search orders of TVAR models #####
disSys <- seq(0.8,1, by = 0.02)
disMea <- seq(0.8,1, by = 0.02)
P <- 5

para_combin <- BayesLatticeLik(signal, P, disSys, disMea)

# BLFscree plot
plot(1:P,para_combin[,4], type = "l",
     xlab = "Order", 
		 ylab = "log(likelihood)", 
		 main = "")
   

##### Procedure 2: Obtain time-varying coefficients, innovation variance, and PARCOR of TVAR(P) #####
sel_order <- 2 ## Remember to change it according to the result of order selection
Dfactor <- para_combin[1:sel_order,2:3] 

tvarp <- BayesLattice(signal, Dfactor)

par(mfrow = c(3,1))
plot(1:N, tvarp$parcor[1,], type="l", main="PARCOR coefficients", xlab="", ylab="", col=1, ylim=c(-1,1))
lines(1:N, tvarp$parcor[2,], col=2)
abline(h=0, col="gray")

ylim <- range(tvarp$coefs)
plot(1:N, tvarp$coefs[1,], type="l", main="TVAR coefficients", xlab="", ylab="", col=1, ylim=ylim)
lines(1:N, tvarp$coefs[2,], col=2)

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

plot(times, signal, type ="l", xlab = "", ylab = "", main = "TVAR(2)")

image(times, freqs, t(tspec), main="True Spectrogram", xlab = "", ylab = "Frequency", col = tim.colors(64))
image.plot(tspec, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)

image(times, freqs, t(spec), main="Posterior mean", xlab = "", ylab = "Frequency", col = tim.colors(64))
image.plot(spec, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)

set.panel() 




    