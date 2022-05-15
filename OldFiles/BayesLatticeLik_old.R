# Search maximum log likelihood values at each lattice stage for TVAR(P) 
#
#	Description:
#
#	Usage:
# 	BayesLatticeLik(signal, P, disSys, disMea)
#
# Arguments:
#   signal: 1 x N dependent variable.
#   P     : maximum order of TVAR models.
#   disSys: a vector of discount factors for system errors (between 0 and 1).
#   disMea: a vector of discount factors for measurement errors (between 0 and 1).  
#
# Values: 
#  para_combin: a P x 4 matrix with columns: order, alpha, disMea, and loglikelihood.  
#
#	Note:
#

BayesLatticeLik = function(signal, P, disSys, disMea){

N <- length(signal)
likp <- array(0, dim = c(length(disSys),length(disMea),P))
para_combin <- matrix(0, nrow = P, ncol = 4)

# initial values at t_0
m0 <- 0.00001 
C0 <- 1
s0 <- var(signal[1:10])
n0 <- 1 

# start the searching
finnov <- signal 
binnov <- signal

for (kk in 1:P) {
	cat("*** Order: ",kk,"\n", sep="")
	flush.console()
	
	maxlike <- -1e300 # give a very small initial value for maxlike
	para_combin[kk,1] <- kk
	
	for (nSys in 1:length(disSys)) {
		for (nMea in  1 : length(disMea)) {
      del <- c(disSys[nSys],disMea[nMea])
			resultf <- dynamicLik(finnov,c(numeric(kk), binnov[1:(N-kk)]), del, m0, C0, s0, n0)
			
			likp[nSys,nMea,kk] <- resultf$llike
			if (resultf$llike > maxlike) {
				para_combin[kk,2] <- disSys[nSys]
				para_combin[kk,3] <- disMea[nMea]
				para_combin[kk,4] <- resultf$llike
				maxlike <- resultf$llike;
			}
		}
	}
	del <- c(para_combin[kk,2], para_combin[kk,3])
	resultf <- dynamicLik(finnov,c(numeric(kk), binnov[1:(N-kk)]), del, m0, C0, s0, n0)
  
	resultb <- dynamicLik(binnov,c(finnov[(kk+1):N], numeric(kk)), del, m0, C0, s0, n0)
		
	finnov <- resultf$e
	binnov <- resultb$e		
}
	return(para_combin)
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
