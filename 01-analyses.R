# run model in nimble using parallel processing
load("./data_5s.Rdata")
library(parallel)
library(coda)
this_cluster <- makeCluster(3)
set.seed(01151929)

# create function to run model
run <- function(seed, dat, const){
  library(nimble)
  code <- nimbleCode( # model code 
    {
      dist.beta[1] ~ dunif(0, 3000)
      for (z in 2:4){ dist.beta[z] ~ dnorm(0, sd=100) }
      dint.beta[1] <- logit(p.dint.beta)
      p.dint.beta ~ dunif(0,1)
      for (j in 2:9){ dint.beta[j] ~ dnorm(0, sd=100) }
      turb.sigma ~ T(dnorm(0,sd=10), 0, )
      phi ~ dunif(0,1)
      psi1 ~ dunif(0,1)
      
      # State model
      for(i in 1:ntracks) {
        z[i,1] ~ dbern(psi1)
        for(t in 2:last[i]) {
          z[i,t] ~ dbern(Ez[i,t-1])
          Ez[i,t-1] <- gamma[i,t-1]*(1-z[i,t-1]) + phi*z[i,t-1]
          gamma[i, t-1] <- dint[i,t-1]*g[i,t-1]
          g[i,t-1] <- exp(-dist[i,t-1]*dist[i,t-1]/(2*dist.sigma[i,t-1]*dist.sigma[i,t-1])) # half-normal distance function
          
          log(dist.sigma[i,t-1]) <- log(dist.beta[1]) + 
            dist.beta[2]*spd[i,t-1] + 
            dist.beta[3]*app[i,t-1] +
            dist.beta[4]*spd[i,t-1]*app[i,t-1] 
          logit(dint[i,t-1]) <-     eps.turb[ turbine[i,t-1] ] + 
            dint.beta[1] + 
            dint.beta[2]*spd[i,t-1] + 
            dint.beta[3]*app[i,t-1] +
            dint.beta[4]*N[i,t-1] +
            dint.beta[5]*E[i,t-1] +
            dint.beta[6]*spd[i,t-1]*app[i,t-1] +
            dint.beta[7]*N[i,t-1]*E[i,t-1] +
            dint.beta[8]*ht[i,t-1] + dint.beta[9]*ht[i,t-1]^2 
        } # time t
      } # track i
      for (k in 1:nturb){ eps.turb[k]~dnorm(0, sd=turb.sigma) } # ntracks
    } 
  ) # end model code
  
  inits <- function(){ list(
    p.dint.beta = runif(1),
    dist.beta = c(runif(1, 50, 500), runif(3, -0.2, 0.2)),
    dint.beta = c(runif(1), runif(8, -0.2, 0.2)),
    phi = runif(1),
    psi1 = runif(1),
    turb.sigma = runif(1)
  ) }
  
  params <- c("p.dint.beta", "dint.beta", "dist.beta",
              "phi", "psi1",
              "turb.sigma", 
              "eps.turb") 
  
  mod<- nimbleModel(code, calculate=T, constants = const, 
                    data = dat, inits = inits())
  
  cmod <- compileNimble(mod)
  mcmc1 <- buildMCMC(cmod)
  cmcmc <- compileNimble(mcmc1)
  chain.out <- runMCMC(cmcmc, 
                       niter = 50000,
                       nburnin = 25000,
                       thin = 25,
                       setSeed = seed)
  
  return(chain.out)
}
# run model in parallel
post <- parLapply(cl = this_cluster, 
                  X = 1:3, 
                  fun = run, 
                  dat=datl, 
                  const=constl)

stopCluster(this_cluster)