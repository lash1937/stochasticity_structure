# code to create appendix S3, figure 1
# comparing community alpha diversity with: 
# demographic stochasticity only, environmental stochasticity only, or both forms of stochasticty.

# source the functions to run each model
source("model_functions.R")

# time to run the model
time <- 50
# number of species
species <- 20

# autocorrelation of environmental noise
# alpha = 0 is uncorrelated (white) noise
alpha <- 0
# scale beta such that sigma has the same variance across all alpha values
beta <- (1-alpha^2)^.5

# create timeseries of environmental conditions
sigma <- rep(NA, time)
sigma[1] <- 0

for (t in 1:(time-1)) {
  sigma[t+1] <- alpha*sigma[t]+beta*rnorm(1,0,1)
}

# set parameters for model run
R <- runif(species, 2, 2.5)
alphas <- matrix(runif(species*species, .005, .01), ncol=species, nrow=species)
diag(alphas) <- .03
zeta <- .35

# single model run (e.g. panel C)
# set up matrix to hold results
results_dem <- results_env <- results_both <- matrix(NA, nrow=species, ncol=time)
results_dem[,1] <- results_env[,1] <- results_both[,1] <- 20

# run models
for (t in 1:(time-1)) {
  for(s in 1:species) {
    results_dem[s,t+1] <- do.com.dem.BH(R=R[s], alpha=alphas[s,], N=results_dem[s,t], Nall=results_dem[,t])
    results_env[s,t+1] <- do.com.env.BH(R=R[s], alpha=alphas[s,], N=results_env[s,t], 
                                        Nall=results_env[,t], zeta=zeta, sigma=sigma[t])
    results_both[s,t+1] <- do.com.both.BH(R=R[s], alpha=alphas[s,], N=results_both[s,t], 
                                          Nall=results_both[,t], zeta=zeta, sigma=sigma[t])
    # set results to 0 if < 1
    # do not need to do so when including demographic stochasticity, as the Poisson distribution
    # will only produce integer values
    if(results_env[s,t+1] < 1) { results_env[s,t+1] <- 0}
  }
}

# determine diversity through time
div_dem <- div_env <- div_both <- rep(NA, time)

for (t in 1:time) {
  div_dem[t] <- sum(results_dem[,t]>1)
  div_env[t] <- sum(results_env[,t]>1)
  div_both[t] <- sum(results_both[,t]>1)
}

# ----------------------------------------------------------------------------------------------------------
# Panel d: creating distributions across multiple runs

# time point to create the distribution
t_star <- 40
# number of runs for creating the distributions
runs <- 1000

# create arrays to hold results
dist_dem <- dist_env <- dist_both <- array(NA, c(runs, species, t_star))
dist_dem[,,1] <- dist_env[,,1] <- dist_both[,,1] <- 20

for (counter in 1:runs) {
  # redraw parameters and recreate environmental timeseries for each run
  R <- runif(species, 2, 2.5)
  alphas <- matrix(runif(species*species, .005, .01), ncol=species, nrow=species)
  diag(alphas) <- .03
  sigma <- rep(NA, time)
  sigma[1] <- 0
  
  for (t in 1:(t_star-1)) {
    sigma[t+1] <- alpha*sigma[t]+beta*rnorm(1,0,1)
    
    for(s in 1:species) {
      dist_dem[counter,s,t+1] <- do.com.dem.BH(R=R[s], alpha=alphas[s,], N=dist_dem[counter,s,t], 
                                               Nall=dist_dem[counter,,t])
      dist_env[counter,s,t+1] <- do.com.env.BH(R=R[s], alpha=alphas[s,], N=dist_env[counter,s,t], 
                                               Nall=dist_env[counter,,t], zeta=zeta, sigma=sigma[t])
      dist_both[counter,s,t+1] <- do.com.both.BH(R=R[s], alpha=alphas[s,], N=dist_both[counter,s,t], 
                                                 Nall=dist_both[counter,,t], zeta=zeta, sigma=sigma[t])
      if(dist_env[counter,s,t+1] < 1) { dist_env[counter,s,t+1] <- 0}
    }
  }
}

# determine diversity at time=t_star for each run
div_dist_dem <- div_dist_env <- div_dist_both <- rep(NA, runs)

for (counter in 1:runs) {
  div_dist_dem[counter] <- sum(dist_dem[counter,,t_star] > 0)
  div_dist_env[counter] <- sum(dist_env[counter,,t_star] > 0)
  div_dist_both[counter] <- sum(dist_both[counter,,t_star] > 0)
}

# create distributions of expected diversity 
density_dem <- density(div_dist_dem, from=0, adjust=1.5)
density_env <- density(div_dist_env, from=0, adjust=1.5)
density_both <- density(div_dist_both, from=0, adjust=1.5)

# ----------------------------------------------------------------------------------------------------------
# plot results
quartz(width=5, height=1.5)
par(mfrow=c(1,3), ps=8, mar=c(3.75,2,1,2), oma=c(1,1,0,5), family="Times")


plot(div_dem,  col="blue", type="l", lwd=2, lty=1, ylab="", yaxt="n", xaxt="n",
     xlab="", cex.axis=1.25, xlim=c(0, 50), ylim=c(0,20), main="")
abline(v=40, lty=2, col="black")
axis(side=1, at=c(0, 25, 50), labels=c(0, 25, 50), cex.axis=1.5)
mtext("Time", side=1, line=2, outer=FALSE, col="black", cex=1.25)
mtext("Diversity", side=2, line=2, outer=FALSE, col="black", cex=1.25)
axis(side=2, at=c(0, 10, 20), labels=c(0, 10, 20), cex.axis=1.5)
lines(div_env, col="darkgreen", lwd=2)
lines(div_both, col="black", lwd=2)

plot(density_dem,  col="blue", lwd=2, ylab="", lty=1, yaxt="n", xaxt="n",
     xlab="", cex.axis=1.25, xlim=c(0, 20), ylim=c(0,.4), main="")
axis(side=1, at=c(0, 10, 20), labels=c(0, 10, 20), cex.axis=1.5)
mtext("Diversity (time=t*)", side=1, line=2, outer=FALSE, col="black", cex=1.25)
mtext("Probability", side=2, line=2, outer=FALSE, col="black", cex=1.25)
axis(side=2, at=c(0, .2, .4), labels=c(0, .2, .4), cex.axis=1.5)
lines(density_env, col="darkgreen", lwd=2, lty=1)
lines(density_both, col="black", lwd=2, lty=1)

plot.new()
legend("center", c("dem", "env", "both"), col=c("blue", "darkgreen", "black"), lwd=3, cex=1.5, bty="n")



