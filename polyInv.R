# from polynomial roots back to polynomial coefficients

x <- polyroot(c(-27,-72,-6,1)) # z0 + z1*x + z2*x^2 + z3*x^3

ee <- x

nn <- length(ee)

pp <- rep(0,nn+1)


pp[1] <- 1

for (ii in 1:nn) {

  pp[2:(ii+1)] <- pp[2:(ii+1)] - ee[ii]*pp[1:ii]

}

Re(pp[length(pp):1])

