# A chirp signal
# 
#
rm(list = ls())

##### load library
library("fields")
library("signal")

##### load R functions of Bayesian lattice filters
source("./Rfun/BayesLatticeLik.R")
source("./Rfun/BayesLattice.R")
source("./Rfun/dynamicLik.R")
source("./Rfun/tvar_spec.R")

## Simulate a chirp signal with quadratic instantaneous frequency deviation
## The chirp is sampled at 1 kHz for 2 seconds. The instantaneous frequency is 100 Hz at t = 0 and crosses 200 Hz at t = 1 second.
Fs <- 1000
times <- seq(0, 2, by=1/Fs)
signal <- chirp(times, 100, 1, 200, "quadratic")

N <- length(signal)

signal <- signal + rnorm(N, 0, 10^-5) ## add noises with small variance to obtain stabilized
                                         ## numerical analysis 

plot(times, signal, main="A chirp signal", type="l", xlab="Time", ylab="")

##### Procedure 1: Search orders of TVAR models #####
disSys <- seq(0.8,1, by = 0.02)
disMea <- seq(0.8,1, by = 0.02)
P <- 25

para_combin <- BayesLatticeLik(signal, P, disSys, disMea)

# BLFscree plot
plot(1:P,para_combin[,4], type = "l",
     xlab = "Order", 
		 ylab = "log(likelihood)", 
		 main = "")

which(((para_combin[2:P,4] - para_combin[1:(P-1),4])/abs(para_combin[1:(P-1),4]))*100 < 0.5)

##### Procedure 2: Obtain time-varying coefficients, innovation variance, and PARCOR of TVAR(P) #####
sel_order <- 12 ## Remember to change it according to the result of order selection   
Dfactor <- para_combin[1:sel_order,2:3] 

tvarp <- BayesLattice(signal, Dfactor)

par(mfrow = c(3,1))
plot(1:N, tvarp$parcor[1,], type="l", main="PARCOR coefficients", xlab="", ylab="", col=1, ylim=c(-1,1))
for (ii in 2:sel_order) lines(1:N, tvarp$parcor[ii,], col=ii)
abline(h=0, col="gray")

ylim <- range(tvarp$coefs)
plot(1:N, tvarp$coefs[1,], type="l", main="TVAR coefficients", xlab="", ylab="", col=1, ylim=ylim)
for (ii in 2:sel_order) lines(1:N, tvarp$coefs[ii,], col=ii)

plot(1:N, tvarp$s2, type="l", main="Variance", xlab="", ylab="", col=1)


##### Procedure 3: Make spectrum plots #####
freqs <- seq(0, 0.5, by = 0.005)
spec <- tvar_spec(tvarp$coefs, tvarp$s2, 1:N, freqs)
spec <- log(spec)

STFT <-  specgram(signal, Fs = Fs) # short time Fourier Transform
logS <- log(abs(STFT$S))

set.panel() # reset plotting device
par(oma=c( 0,0,0,6)) # margin of 4 spaces width at right hand side
set.panel( 3,1) 

plot(times, signal, type ="l", xlab = "", ylab = "", main = "TVAR(6)")

image(STFT$t, STFT$f, t(logS), main="STFT", xlab = "", ylab = "Frequency (Hz)", col = tim.colors(64))
image.plot(logS, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)

image(times, freqs*Fs, t(spec), main="Posterior mean", xlab = "", ylab = "Frequency (Hz)", col = tim.colors(64))
image.plot(spec, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)

set.panel() 


