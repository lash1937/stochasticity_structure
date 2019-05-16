# code to create figure 4,
# diving into the effects of environmental stochasticity on population and community structure

# source the functions to run each model
source("model_functions.R")

# ----------------------------------------------------------------------------------------------------------
# single run population results (panel B)

# time to run the model
time <- 100

# compare no autocorrelation (alpha=0) to positive autocorrelation (alpha=.75)
alpha <- c(0, .75)
beta <- (1-alpha^2)^.5

# create timeseries of environmental conditions
sigma_no <- rep(NA, time)
sigma_no[1] <- 0
sigma_pos <- sigma_no

for (t in 1:(time-1)) {
  sigma_no[t+1] <- alpha[1]*sigma_no[t]+beta[1]*rnorm(1,0,1)
  sigma_pos[t+1] <- alpha[2]*sigma_pos[t]+beta[2]*rnorm(1,0,1)
}

# density independent growth rates
R1 <- 1.5     # speices 1
R2 <- 1.5     # speices 1
# intraspecific competition coefficients
alpha1 <- .05
alpha2 <- .05
# magnitude of the effect of env. stochasticity
zeta <- .25

# set up vector to hold results
results_no <- results_pos <- rep(NA, time)
results_no[1] <- results_pos[1] <- 20

# run model
for (t in 1:(time-1)) {
  results_no[t+1] <- do.pop.env.BH(R=R1, alpha=alpha1, N=results_no[t], zeta=zeta, sigma=sigma_no[t])
  results_pos[t+1] <- do.pop.env.BH(R=R2, alpha=alpha2, N=results_pos[t], zeta=zeta, sigma=sigma_pos[t])
  # if population sizes are below 1 individual, set to 0
  if(results_no[t+1] < 1) { results_no[t+1] <- 0}
  if(results_pos[t+1] < 1) { results_pos[t+1] <- 0}
}

# Create figure (panel B)
quartz(width=5, height=5)
par(mar=c(3,3,2,2))

plot(results_no,  col="darkgoldenrod2", type="l", lwd=2, lty=6, ylab="", yaxt="n", xaxt="n",
     xlab="", cex.axis=1.25, xlim=c(0, 100), ylim=c(0,30), main="")
#abline(v=75, lty=2, col="black")
axis(side=1, at=c(0, 50, 100), labels=c(0, 50, 100), cex.axis=1.25)
mtext("Time", side=1, line=2, outer=FALSE, col="black", cex=1.25)
mtext("Abundance", side=2, line=2, outer=FALSE, col="black", cex=1.25)
axis(side=2, at=c(0, 15, 30), labels=c(0, 15, 30), cex.axis=1.25)
lines(results_pos, col="firebrick4", lwd=2)

# ----------------------------------------------------------------------------------------------------------
# across runs population distributions (panel C)

# number of runs for creating the distributions
runs  <- 1000

# set up matrix to hold timeseries results
dist_no <- dist_pos <- matrix(NA, nrow=runs, ncol=time)
dist_no[,1] <- dist_pos[,1] <- 20

# set up vector to hold results of correlations in each timeseries
# between time step t and t-1
auto_no <- auto_pos <- rep(NA, runs)

# run model
for (counter in 1:runs) {
  # create new timeseries of environmental conditions in each run
  sigma_no <- rep(NA, time)
  sigma_no[1] <- 0
  sigma_pos <- sigma_no
  
  for (t in 1:(time-1)) {
    sigma_no[t+1] <- alpha[1]*sigma_no[t]+beta[1]*rnorm(1,0,1)
    sigma_pos[t+1] <- alpha[2]*sigma_pos[t]+beta[2]*rnorm(1,0,1)
    
    # determine population dynamics
    dist_no[counter,t+1] <- do.pop.env.BH(R=R1, alpha=alpha1, 
                                          N=dist_no[counter,t], zeta=zeta, sigma=sigma_no[t])
    dist_pos[counter,t+1] <- do.pop.env.BH(R=R2, alpha=alpha2, 
                                           N=dist_pos[counter,t], zeta=zeta, sigma=sigma_pos[t])
    if(dist_no[counter, t+1] < 1) { dist_no[counter, t+1] <- 0}
    if(dist_pos[counter, t+1] < 1) { dist_pos[counter, t+1] <- 0}
  }
  # calculate correlation in population dynamics
  auto_no[counter] <- cor(dist_no[counter,1:(time-1)], dist_no[counter,2:time])
  auto_pos[counter] <- cor(dist_pos[counter,1:(time-1)], dist_pos[counter,2:time])
}

# create distributions of expected correlations
density_results1 <- density(auto_no)
density_results2 <- density(auto_pos)

# Create figure (panel C)
quartz(width=5, height=5)
par(mar=c(3,3,2,2))

plot(density_results1,  col="darkgoldenrod2", lwd=2, ylab="", lty=6, yaxt="n", xaxt="n",
     xlab="", cex.axis=1.25, xlim=c(0, 1), ylim=c(0,14), main="")
axis(side=1, at=c(0, .25, .5, .75, 1), labels=TRUE, cex.axis=1.25)
mtext("Autocorrelation", side=1, line=2, outer=FALSE, col="black", cex=1.25)
mtext("Probability Density", side=2, line=2, outer=FALSE, col="black", cex=1.25)
axis(side=2, at=c(0, 7, 14), labels=TRUE, cex.axis=1.25)
lines(density_results2, col="firebrick4", lwd=2)

# ----------------------------------------------------------------------------------------------------------
# single run community results (panel D)

# number of species
species <- 20

# set parameters for model run
# magnitude of the effect of env. stochasticity
zeta <- .2
# survival in the seedbank
surv <- .8
# fraction of seeds that germinate under environmental conditions at t=0
g <- .5
# density indepdendent growth rates
R_stabilizing <- runif(species, 1.1, 1.5)
# interspecific competition coefficients
alphas_stabilizing <- matrix(runif(species*species, .001, .005), ncol=species, nrow=species)
# intraspecific competition coefficients
diag(alphas_stabilizing) <- .007

# set up matrix to hold results
results_no <- matrix(NA, nrow=species, ncol=time)
results_no[,1] <- 20
results_pos <- results_no

# create timeseries of environmental conditions per each species
# i.e. each species responds independently (rather than synchronously)
sigma_no <- matrix(NA, species, time)
sigma_no[,1] <- 0
sigma_pos <- sigma_no

# run model
for (t in 1:(time-1)) {
  for (s in 1:species) {
    sigma_no[s, t+1] <- alpha[1]*sigma_no[s, t]+beta[1]*rnorm(1,0,1)
    sigma_pos[s, t+1] <- alpha[2]*sigma_pos[s, t]+beta[2]*rnorm(1,0,1)
    
    results_no[s,t+1] <- do.com.seedbank(R=R_stabilizing[s], alphas=alphas_stabilizing[s,], 
                                         N=results_no[s,t], Nall=results_no[,t],
                                         zeta=zeta, sigma=sigma_no[s,t], s=surv, g=g)
    results_pos[s,t+1] <- do.com.seedbank(R=R_stabilizing[s], alphas=alphas_stabilizing[s,], 
                                          N=results_pos[s,t], Nall=results_pos[,t],
                                          zeta=zeta, sigma=sigma_pos[s,t], s=surv, g=g)
    if(results_pos[s, t+1] < 1) { results_pos[s, t+1] <- 0}
    if(results_no[s, t+1] < 1) { results_no[s, t+1] <- 0}
  }
}

# calculate diversity through time
diversity_no<- diversity_pos <- rep(NA, time)

for (t in 1:time) {
  diversity_no[t] <- sum(results_no[,t] > 0)
  diversity_pos[t] <- sum(results_pos[,t] > 0)
}

# Create figure (panel D)
quartz(width=5, height=5)
par(mar=c(3,3,2,2))

plot(diversity_no,  col="plum2", type="l", lwd=2, lty=6, ylab="", yaxt="n", xaxt="n",
     xlab="", cex.axis=1.25, xlim=c(0, 100), ylim=c(0,20), main="")
abline(v=75, lty=2, col="black")
axis(side=1, at=c(0, 50, 100), labels=TRUE, cex.axis=1.25)
mtext("Time", side=1, line=2, outer=FALSE, col="black", cex=1.25)
mtext("Diversity", side=2, line=2, outer=FALSE, col="black", cex=1.25)
axis(side=2, at=c(0, 10, 20), labels=c(0, 10, 20), cex.axis=1.25)
lines(diversity_pos, col="darkorchid3", lwd=2)

# ----------------------------------------------------------------------------------------------------------
# distribution for community results (panel E)

# time point to create the distribution
t_star <- 75

# set up arrays to hold results
results_no <- array(NA, c(runs, species, time))
results_no[,,1] <- 20
results_pos <- results_no

for (counter in 1:runs) {
  # redraw parameters for each run
  R_stabilizing <- runif(species, 1.1, 1.5)
  alphas_stabilizing <- matrix(runif(species*species, .001, .005), ncol=species, nrow=species)
  diag(alphas_stabilizing) <- .007
  
  # redraw enviromental conditions for each run
  sigma_no <- matrix(NA, species, time)
  sigma_no[,1] <- 0
  sigma_pos <- sigma_no
  
  for (t in 1:(time-1)) {
    
    for (s in 1:species) {
      sigma_no[s, t+1] <- alpha[1]*sigma_no[s, t]+beta[1]*rnorm(1,0,1)
      sigma_pos[s, t+1] <- alpha[2]*sigma_pos[s, t]+beta[2]*rnorm(1,0,1)
      
      results_no[counter,s,t+1] <- do.com.seedbank(R=R_stabilizing[s], alphas=alphas_stabilizing[s,], 
                                                   N=results_no[counter,s,t], Nall=results_no[counter,,t],
                                                   zeta=zeta, sigma=sigma_no[s,t], s=surv, g=g)
      results_pos[counter,s,t+1] <- do.com.seedbank(R=R_stabilizing[s], alphas=alphas_stabilizing[s,], 
                                                    N=results_pos[counter,s,t], Nall=results_pos[counter,,t],
                                                    zeta=zeta, sigma=sigma_pos[s,t], s=surv, g=g)
      
      if(results_pos[counter,s, t+1] < 1) { results_pos[counter,s, t+1] <- 0}
      if(results_no[counter,s, t+1] < 1) { results_no[counter,s, t+1] <- 0}
    }
  }
}

# extract results at t=t_star
t_star_pos <- results_pos[,,t_star]
t_star_no <- results_no[,,t_star]

# determine diversity at t_star for each run
diversity_no <- diversity_pos <- rep(NA, runs)

for(counter in 1:runs) {
  diversity_no[counter] <- sum(t_star_no[counter,]>0)
  diversity_pos[counter] <- sum(t_star_pos[counter,]>0)
}

# create distributions of expected diversity
density_pos <- density(diversity_pos, from=0, adjust=1.5)
density_no <- density(diversity_no, from=0, adjust=1.5)

# Create figure (panel E)
quartz(width=5, height=5)
par(mar=c(3,3,2,2))

plot(density_pos,  col="darkorchid3", lwd=2, ylab="", lty=1, yaxt="n", xaxt="n",
     xlab="", cex.axis=1.25, xlim=c(0, 20), ylim=c(0,.3), main="")
axis(side=1, at=c(0, 10, 20), labels=c(0, 10, 20), cex.axis=1.25)
mtext("Diversity (time=t*)", side=1, line=2, outer=FALSE, col="black", cex=1.25)
mtext("Probability", side=2, line=2, outer=FALSE, col="black", cex=1.25)
axis(side=2, at=c(0, .15, .3), labels=TRUE, cex.axis=1.25)
lines(density_no, col="plum2", lwd=2, lty=6)

