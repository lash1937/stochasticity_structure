# All population and community model functions

# deterministic, population model
# Equation 1 from text
# R = growth rate
# alpha = intraspecific competition coefficient
# N = population size at time t
# N_next = population size at time t+1
do.pop.det.BH <- function(R, alpha, N) {
  N_next <- R*N/(1+alpha*N)
  return(N_next)
}

# demographic stochasticity, population model
# Equation 2 from text
# R = growth rate
# alpha = intraspecific competition coefficient
# N = population size at time t
# N_next = population size at time t+1
do.pop.dem.BH <- function(R, alpha, N) {
  N_next <- rpois(1, R*N/(1+alpha*N))
  return(N_next)
}

# environmental stochasticity, population model
# Equation 3 from text
# R = growth rate
# alpha = intraspecific competition coefficient
# N = population size at time t
# N_next = population size at time t+1
# zeta = magnitude of the effect of env. stochasticity
# sigma = environmental condition at time t
do.pop.env.BH <- function(R, alpha, N, zeta, sigma) {
  N_next <- R*N/(1+alpha*N) + N*zeta*sigma
  return(N_next)
}

# environmental and demographic stochasticity, population model
# Equation 4 from text
# R = growth rate
# alpha = intraspecific competition coefficient
# N = population size at time t
# N_next = population size at time t+1
# zeta = magnitude of the effect of env. stochasticity
# sigma = environmental condition at time t
do.pop.both.BH <- function(R, alpha, N, zeta, sigma) {
  N_next <- rpois(1,R*N/(1+alpha*N) + N*zeta*sigma)
  return(N_next)
}

# demographic stochasticity, community model
# Equation 6 from text
# R = growth rate
# alphas = competition coefficients
# Nall = vector of population sizes of all species
# N = focal species' population size at time t
# N_next = population size at time t+1
do.com.dem.BH <- function(R, alphas, N, Nall) {
  N_next <- rpois(1, R*N/(1+sum(alphas*Nall)))
  return(N_next)
}

# environmental stochasticity, community model
# Equation 7 from text
# R = growth rate
# alphas = competition coefficients
# Nall = vector of population sizes of all species
# N = focal species' population size at time t
# N_next = population size at time t+1
# zeta = magnitude of the effect of env. stochasticity
# sigma = environmental condition at time t
do.com.env.BH <- function(R, alphas, N, Nall, zeta, sigma) {
  N_next <- R*N/(1+sum(alphas*Nall)) + N*zeta*sigma
  return(N_next)
}

# environmental and demographic stochasticity, community model
# Equation 8 from text
# R = growth rate
# alphas = competition coefficients
# Nall = vector of population sizes of all species
# N = focal species' population size at time t
# N_next = population size at time t+1
# zeta = magnitude of the effect of env. stochasticity
# sigma = environmental condition at time t
do.com.both.BH <- function(R, alphas, N, Nall, zeta, sigma) {
  N_next <- rpois(1, max(R*N/(1+sum(alphas*Nall)) + N*zeta*sigma,0))
  return(N_next)
}

# community model with seedbank and environmental stochasticity
# Equation 10 from text
# R = growth rate
# alphas = competition coefficients
# Nall = vector of population sizes of all species
# N = focal species' population size at time t
# N_next = population size at time t+1
# zeta = magnitude of the effect of env. stochasticity
# sigma = environmental condition at time t
# s = seedbank survival
# g = fraction of seeds that germinate under environmental conditions at t=0
do.com.seedbank <- function(R, alphas, N, Nall, zeta, sigma, s, g) {
  germ <- g+zeta*sigma
  if(germ < 0) {germ <- 0}
  if(germ>1) {germ <- 1}
  N_next <- N*s*(1-germ) +germ*R*N/(1+sum(alphas*Nall))
  return(N_next)
}

