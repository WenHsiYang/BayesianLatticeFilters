# Compute log-likelihood values of dynamic linear models
#
# Description:
#		
# *** Model:
#       obs equation: f_t = phi_t * b_t + e_t,   e_t ~ N(0,sig_t^2)
#
#     	trans equation: phi_t = phi_{t-1} + w_t,   w_t ~ N(0,W_t)
#
# *** Assumption:
#     	evolution error variance: sig_t^2 = (sig_{t-1}^2 * delta)/eta_t,   
#                         
#                              eta_t ~ Beta(delta*k_{t-1}/2, (1-delta)*k_{t-1}/2)
#
#     	evolution sys variance: W_t = C_{t-1}(1-alpha)/alpha
# 
# *** Initial prior:
#     	phi_0 ~ N(m_0, C_0)
#     	sig_0^{-2} ~ Gamma(k_0/2, d_0/2)
#
#	Usage:
#		dynamicLik(y, x, del, m0, C0, s0, n0) 
#
# Arguments:
#   y	 : 1 x N dependent variable
#   x  : 1 x N independent variable 
# 	del: Discount factors: del[1] = alpha and del[2] = delta  
#  	m0 : initial prior mean for state
#  	C0 :  initial variance for state
#  	s0 :  initial prior estimate of obs var
#  	n0 :  initial prior df
#
# Values: 
#  m    : 1 x N post mean 
#  s    : 1 x N post obs var estimates
#  e    : 1 x N estimated innovations 
#  llike: log-likelihood
#
# Note:
#		

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
