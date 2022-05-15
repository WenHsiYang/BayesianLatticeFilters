# Estimate coefficients, innovation variances, and PAROCR coefficients
# of TVAR(P) 
#
#	Description:
# 
#	Usage:
#		BayesLattice(signal, Dfactor)
#
# Arguments:
#   signal :  1 x N vector dependent variable
#   Dfactor:  P x 2 matrix with columns: gamma and delta
#
# Values: 
#  coefs : P x N matrix of estimated coefficients of TVAR models 
#  s2  	 : 1 x N vector of estimated innovation variances  
#  parcor: P x N matrix of estimated PARCOR
#
#	Note:
#

BayesLattice = function(signal, Dfactor){

N <- length(signal)
P <- nrow(Dfactor)

fcoefs <- array(0, dim = c(P,P,N))
bcoefs <- array(0, dim = c(P,P,N))

## initial priors
m0 <- 0.00001 
C0 <- 1
s0 <- var(signal[1:10]) 
n0 <- 1 

finnov <- signal 
binnov <- signal
del <- Dfactor[1,]

cat("*** Order: 1\n", sep="")
flush.console()
	
# get the instantaneous forward PARCOR
resultf <- dynamicLik(finnov, c(0, binnov[1:(N-1)]), del, m0, C0, s0, n0)
#s2 <- resultf$s

# get the instantaneous back PARCOR
resultb <- dynamicLik(binnov,c(finnov[2:N], 0), del, m0, C0, s0, n0);

fcoefs[1,1,] <- resultf$m
bcoefs[1,1,] <- resultb$m 


if (P > 1) {
	finnov <- resultf$e
	binnov <- resultb$e
	
	for (kk in 2:P) {
  	cat("*** Order: ",kk,"\n", sep="")
		flush.console()
	
		del <- Dfactor[kk,]

		# get the instantaneous forward PARCOR
		resultf <- dynamicLik(finnov,c(numeric(kk), binnov[1:(N-kk)]), del, m0, C0, s0, n0)
		#s2 <- resultf$s
	
		# get the instantaneous back PARCOR
		resultb <- dynamicLik(binnov,c(finnov[(kk+1):N], numeric(kk)), del, m0, C0, s0, n0)
		
		fcoefs[kk,kk,] <- resultf$m
		bcoefs[kk,kk,] <- resultb$m
    
		ii <- kk - 1
		for (jj in 1:(kk-1)) {
			hh <- kk - jj
			fcoefs[kk,jj,] <- fcoefs[ii,jj,] - fcoefs[kk,kk,] * bcoefs[ii,hh,]
		}
    
		ii <- kk - 1
		for (jj in 1:(kk-1)) {
			hh <- kk - jj
			bcoefs[kk,jj,] <- bcoefs[ii,jj,] - bcoefs[kk,kk,] * fcoefs[ii,hh,]
		}
    
		finnov <- resultf$e
		binnov <- resultb$e
  }
}

s2 <- resultf$s

parcor <- matrix(0, nrow = P, ncol = N)
for (ii in 1:P) {
    parcor[ii,] <- fcoefs[ii,ii,]
}

coefs <- matrix(0, nrow = P, ncol = N)
for (ii in 1:P) {
    coefs[ii,] <- fcoefs[P,ii,]
}

return(list(coefs = coefs, parcor = parcor, s2 = s2))
}

############ subfunctions ####################
# Please see the details of dynamicLik
dynamicLik = function(y, x, del, m0, C0, s0, n0){

# set output
N <- length(y)
m <- numeric(N)  
CC <- numeric(N)  
s <- numeric(N)    
n <- numeric(N)  

d <- del[1] 
b <- del[2]
mt <- m0 
Ct <- C0 
st <- s0 
nt <- n0
llike <- 0
# forward filtering
for (ii in 1:N) {
  G <- x[ii]
  A <- Ct*G/d 
	Q <- G*A + st 
	A <- A/Q 
	e <- y[ii] - G * mt
	nt <- b * nt
	# non-standardized Student's t-distribution
  llike <- llike + lgamma((nt+1)/2) - lgamma(nt/2) - log(nt*Q)/2 - (nt+1)*log(1+e*e/(Q*nt))/2 
    
	mt <- mt + A*e 
	m[ii] <- mt
	r <- nt + e*e/Q 
	nt <- nt+1 
	r <- r/nt 
	st <- st*r 
	n[ii] <- nt 
	s[ii] <- st   
  Ct <- r * (Ct/d - A*A*Q)
	CC[ii] <- Ct 	
}

# backward smoothing
e <- numeric(N)  
for (ii in (N-1):1) {
  m[ii] <- (1-d)*m[ii] + d*m[ii+1]
  e[ii] <- y[ii]- m[ii]*x[ii]
     
	n[ii] <- (1-b)*n[ii] + b*n[ii+1]  
  st <- s[ii] 
	s[ii] <- 1 / ((1-b)/st + b/s[ii+1]) 
  CC[ii] <- s[ii] * ((1-d)*CC[ii]/st + d*d*CC[ii+1]/st) 
}
	return(list(m = m, CC = CC, s = s, e = e, n = n, llike = llike))
}
