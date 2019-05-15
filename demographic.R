# code to create figure 3,
# diving into the effects of demographic stochasticity on population and community structure

# source the functions to run each model
source("model_functions.R")

# ----------------------------------------------------------------------------------------------------------
# single run population results (panel B)

# time to run the model
time <- 50
# density independent growth rates
R1 <- 1.6     # speices 1
R2 <- 1.6     # species 2
# intraspecific competition coefficients
alpha1 <- .02
alpha2 <- .1

# set up vector to hold results
results1 <- results2 <- rep(NA, time)
results1[1] <- results2[1] <- 20

# run model
for (t in 1:(time-1)) {
  results1[t+1] <- do.pop.dem.BH(R=R1, alpha=alpha1, N=results1[t])
  results2[t+1] <- do.pop.dem.BH(R=R2, alpha=alpha2, N=results2[t])
}

# Create figure (panel B)
quartz(width=5, height=5)
par(mar=c(3,3,2,2))

plot(results2,  col="darkgoldenrod2", type="l", lwd=2, lty=1, ylab="", yaxt="n", xaxt="n",
     xlab="", cex.axis=1.25, xlim=c(0, 50), ylim=c(0,50), main="")
abline(v=40, lty=2, col="black")
axis(side=1, at=c(0, 25, 50), labels=c(0, 25, 50), cex.axis=1.25)
mtext("Time", side=1, line=2, outer=FALSE, col="black", cex=1.25)
mtext("Abundance", side=2, line=2, outer=FALSE, col="black", cex=1.25)
axis(side=2, at=c(0, 25, 50), labels=c(0, 25, 50), cex.axis=1.25)
lines(results1, col="firebrick4", lwd=2)

# ----------------------------------------------------------------------------------------------------------
# across runs population distributions (panel C)

# time point to create the distribution
t_star <- 40
# number of runs for creating the distributions
runs  <- 1000

# create matrix to hold results
results1 <- results2 <- matrix(NA, nrow=runs, ncol=time)
results1[,1] <- results2[,1] <- 20

# run model
for (counter in 1:runs) {
  for (t in 1:(time-1)) {
    results1[counter,t+1] <- do.pop.dem.BH(R=R1, alpha=alpha1, N=results1[counter,t])
    results2[counter,t+1] <- do.pop.dem.BH(R=R2, alpha=alpha2, N=results2[counter,t])
  }
}

# extract abundance at time=t_star for each run
dist_results1 <- results1[,t_star]
dist_results2 <- results2[,t_star]

# create distributions of expected diversity 
density_results1 <- density(dist_results1, from=0)
density_results2 <- density(dist_results2, from=0)

# Create figure (panel C)
quartz(width=5, height=5)
par(mar=c(3,3,2,2))

plot(density_results2,  col="darkgoldenrod2", lwd=2, ylab="", lty=1, yaxt="n", xaxt="n",
     xlab="", cex.axis=1.25, xlim=c(0, 50), ylim=c(0,.4), main="")
axis(side=1, at=c(0, 25, 50), labels=c(0, 25, 50), cex.axis=1.25)
mtext("Abundance (time=t*)", side=1, line=2, outer=FALSE, col="black", cex=1.25)
mtext("Probability", side=2, line=2, outer=FALSE, col="black", cex=1.25)
axis(side=2, at=c(0, .2, .4), labels=c(0, .2, .4), cex.axis=1.25)
lines(density_results1, col="firebrick4", lwd=2)

# ----------------------------------------------------------------------------------------------------------
# single run community results (panel D)

# number of species
species <- 20

# set parameters for model run
# parameters for the stabilizing case
R_stabilizing <- runif(species, 2, 2.5)
alphas_stabilizing <- matrix(runif(species*species, .002, .005), ncol=species, nrow=species)
diag(alphas_stabilizing) <- .03

# parameters for the neutral case
R_neutral <- rep(mean(R_stabilizing), species)
alphas_neutral <- matrix(rep(mean(alphas_stabilizing), species*species), ncol=species, nrow=species)


# set up matrix to hold results
results_stabilizing <- results_neutral <- matrix(NA, nrow=species, ncol=time)
results_stabilizing[,1] <- results_neutral[,1] <- 20

# run model
for (t in 1:(time-1)) {
  for (s in 1:species) {
    results_stabilizing[s,t+1] <- do.com.dem.BH(R=R_stabilizing[s], alphas=alphas_stabilizing[s,], 
                                                N=results_stabilizing[s,t], Nall=results_stabilizing[,t])
    results_neutral[s,t+1] <- do.com.dem.BH(R=R_neutral[s], alphas=alphas_neutral[s,], 
                                            N=results_neutral[s,t], Nall=results_neutral[,t])
  }
}

# determine diversity through time
diversity_stabilizing <- diversity_neutral <- rep(NA, time)

for (t in 1:time) {
  diversity_stabilizing[t] <- sum(results_stabilizing[,t] > 0)
  diversity_neutral[t] <- sum(results_neutral[,t] > 0)
}

# Create figure (panel D)
quartz(width=5, height=5)
par(mar=c(3,3,2,2))

plot(diversity_stabilizing,  col="darkorchid3", type="l", lwd=2, lty=1, ylab="", yaxt="n", xaxt="n",
     xlab="", cex.axis=1.25, xlim=c(0, 50), ylim=c(0,20), main="")
abline(v=40, lty=2, col="black")
axis(side=1, at=c(0, 25, 50), labels=c(0, 25, 50), cex.axis=1.25)
mtext("Time", side=1, line=2, outer=FALSE, col="black", cex=1.25)
mtext("Diversity", side=2, line=2, outer=FALSE, col="black", cex=1.25)
axis(side=2, at=c(0, 10, 20), labels=c(0, 10, 20), cex.axis=1.25)
lines(diversity_neutral, col="plum2", lwd=2)

# ----------------------------------------------------------------------------------------------------------
# distribution for community results (panel E)

# set up arrays to hold results
dist_stabilizing <- dist_neutral <- array(NA, c(runs, species, time))
dist_stabilizing[,,1] <- dist_neutral[,,1] <- 20

# model run
for (counter in 1:runs) {
  # redraw parameters for each run
  R_stabilizing <- runif(species, 2, 2.5)
  alphas_stabilizing <- matrix(runif(species*species, .002, .005), ncol=species, nrow=species)
  diag(alphas_stabilizing) <- .03
  
  R_neutral <- rep(mean(R_stabilizing), species)
  alphas_neutral <- matrix(rep(mean(alphas_stabilizing), species*species), ncol=species, nrow=species)
  
  for (t in 1:(time-1)) {
    for (s in 1:species) {
      dist_stabilizing[counter,s,t+1] <- do.com.dem.BH(R=R_stabilizing[s], alphas=alphas_stabilizing[s,], 
                                                       N=dist_stabilizing[counter,s,t], 
                                                       Nall=dist_stabilizing[counter,,t])
      
      dist_neutral[counter,s,t+1] <- do.com.dem.BH(R=R_neutral[s], alphas=alphas_neutral[s,], 
                                                   N=dist_neutral[counter,s,t], 
                                                   Nall=dist_neutral[counter,,t])
    }
  }
}

# extract results for t=t_star
t_star_stabilizing <- dist_stabilizing[,,t_star]
t_star_neutral <- dist_neutral[,,t_star]

# determine diversity at t=t_star for each run
diversity_stabilizing <- diversity_neutral <- rep(NA, runs)

for(counter in 1:runs) {
  diversity_stabilizing[counter] <- sum(t_star_stabilizing[counter,]>0)
  diversity_neutral[counter] <- sum(t_star_neutral[counter,]>0)
}

# create distributions of expected diversity
density_stabilizing <- density(diversity_stabilizing, from=0, adjust=1.5)
density_neutral <- density(diversity_neutral, from=0, adjust=1.5)

# Create figure (panel E)
quartz(width=5, height=5)
par(mar=c(3,3,2,2))

plot(density_stabilizing,  col="darkorchid3", lwd=2, ylab="", lty=1, yaxt="n", xaxt="n",
     xlab="", cex.axis=1.25, xlim=c(0, 20), ylim=c(0,.3), main="")
axis(side=1, at=c(0, 10, 20), labels=c(0, 10, 20), cex.axis=1.25)
mtext("Diversity (time=t*)", side=1, line=2, outer=FALSE, col="black", cex=1.25)
mtext("Probability", side=2, line=2, outer=FALSE, col="black", cex=1.25)
axis(side=2, at=c(0, .15, .3), labels=c(0, .15, .3), cex.axis=1.25)
lines(density_neutral, col="plum2", lwd=2)


