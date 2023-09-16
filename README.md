# Bayesian Lattice Filters (BLF)

Estimate the parameters of time-varying autoregression models and reveal the time-frequency representations of signals.

## Major functions

 - tvar_spec: Compute the time-varying representation of a TVAR(p) model.
 - dynamicLik: Compute log-likelihood values of dynamic linear models.
 - BayesLatticeLik: Compute log-likelihood values of dynamic linear models.
 - BayesLattice: Estimate coefficients, innovation variances, and PAROCR coefficients of a TVAR(p) model.
 
## Experience the efficiency and robustness of BLF
  
  - Example 1: TVAR(2)
  - Example 2: TVAR(6)
  - Example 3: Chirp signals

## Further reading
* Yang, W. H., Holan, S. H., & Wikle, C. K. (2016). Bayesian lattice filters for time-varying autoregression and time–frequency analysis. Bayesian Analysis, 11(4), 977-1003. [link](https://projecteuclid.org/journals/bayesian-analysis/volume-11/issue-4/Bayesian-Lattice-Filters-for-Time-Varying-Autoregression-and-TimeFrequency-Analysis/10.1214/15-BA978.full?tab=ArticleLink)

  * Note: in the paper, Equation (8) contains a typo; that is, $d_{t,k}^{(m)}$ on the right side of the equation should be replaced with $d_{t,m}^{(m)}$.