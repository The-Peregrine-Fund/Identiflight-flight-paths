# Run models in nimble
##################
# run global model
##################
load("./data_10s.Rdata")
library(nimble)

code <- nimbleCode( # start model code
  {
    dist.beta[1] ~ dunif(0, 3000)
    dint.beta[1] <- logit(p.dint.beta)
    p.dint.beta ~ dunif(0,1)
    for (j in 2:13){ dint.beta[j] ~ dnorm(0, sd=100) }
    for (z in 2:11){ dist.beta[z] ~ dnorm(0, sd=100) }
    turb.sigma ~ dunif(0,20)
    month.sigma ~ dunif(0,20)
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
          dist.beta[2]*Sn[i,t-1] + dist.beta[3]*RotSp[i,t-1] + dist.beta[4]*perp[i,t-1] + dist.beta[5]*app[i,t-1] +
          dist.beta[6]*Sn[i,t-1]*RotSp[i,t-1] + dist.beta[7]*Sn[i,t-1]*perp[i,t-1]+ dist.beta[8]*Sn[i,t-1]*app[i,t-1] +
          dist.beta[9]*RotSp[i,t-1]*perp[i,t-1] + dist.beta[10]*RotSp[i,t-1]*app[i,t-1] + 
          dist.beta[11]*perp[i,t-1]*app[i,t-1]
        logit(dint[i,t-1]) <- dint.beta[1] + eps.turb[ turbine[i,t-1] ] + eps.month[ month[i] ] +
          dint.beta[2]*Sn[i,t-1] + dint.beta[3]*RotSp[i,t-1] + dint.beta[4]*perp[i,t-1] + dint.beta[5]*app[i,t-1] +
          dint.beta[6]*Sn[i,t-1]*RotSp[i,t-1] + dint.beta[7]*Sn[i,t-1]*perp[i,t-1]+ dint.beta[8]*Sn[i,t-1]*app[i,t-1] +
          dint.beta[9]*RotSp[i,t-1]*perp[i,t-1] + dint.beta[10]*RotSp[i,t-1]*app[i,t-1] + 
          dint.beta[11]*perp[i,t-1]*app[i,t-1] +
          dint.beta[12]*ht[i,t-1] + dint.beta[13]*ht2[i,t-1] 
      } # time t
    } # track i
    for (k in 1:nturb){ eps.turb[k]~dnorm(0, sd=turb.sigma) } # ntracks
    for (m in 1:nmonth){ eps.month[m]~dnorm(0, sd=month.sigma) }
  } 
) # end model code

z.init <- datl$z
z.init[is.na(z.init)] <- as.numeric(rbinom(prob=0.5, size=1, n=sum(is.na(z.init)) )) 
z.init[!is.na(datl$z)] <- NA

constl$turbine[] <- as.numeric(factor(constl$turbine))

inits <- function(){ list(
  z=z.init,
  p.dint.beta = runif(1),
  dist.beta = c(runif(1, 50, 500), runif(10, -0.2, 0.2)),
  dint.beta = c(runif(1), runif(12, -0.2, 0.2)),
  eps.turb = runif(constl$nturb, -0.1, 0.1),
  eps.month = runif(constl$nmonth, -0.1, 0.1),
  month.sigma = runif(1),
  phi = runif(1),
  psi1 = runif(1),
  turb.sigma = runif(1),
  dist.sigma = array(runif(constl$ntracks*(constl$ntime-1), 100, 1000) , dim=c(constl$ntracks, (constl$ntime-1))),
  dint = array(runif(constl$ntracks*(constl$ntime-1)) , dim=c(constl$ntracks, (constl$ntime-1))),
  Ez = array(runif(constl$ntracks*(constl$ntime-1)) , dim=c(constl$ntracks, (constl$ntime-1))),
  gamma = array(runif(constl$ntracks*(constl$ntime-1)) , dim=c(constl$ntracks, constl$ntime-1)),
  g= array(runif(constl$ntracks*(constl$ntime-1)) , dim=c(constl$ntracks, (constl$ntime-1)))
) }

params <- c("p.dint.beta", "dint.beta", "dist.beta",
            "phi", "psi1",
            "turb.sigma", "month.sigma",
            "eps.turb", "eps.month")

ni <- 50000; nc <- 3;  nb <- 25000; nt <- 25; nt2 <- 20

mod<- nimbleModel(code, calculate=T, constants = constl, 
                  data = datl, inits = inits())

out <- nimbleMCMC(  model=mod,
                    code = code,
                    monitors = params,
                    nchains=nc,
                    thin = nt,
                    niter = ni,
                    nburnin = nb,
                    progressBar=T,
                    summary=T,
                    WAIC=F,
                    samplesAsCodaMCMC = T,
                    samples=T )                                                      nburnin = nb))

save(mod, out, 
     file="./results_Bayes_10s.Rdata")

##################
# cross validation
##################
library(nimble)
load("./data_10s_cv.Rdata")
for (xx in 1:10){
  datl <- alldat[[xx]]
  constl <- allconst[[xx]]
  datl$turbine <- array(as.numeric(factor(constl$turbine)), dim=dim(constl$turbine))
  datl$month <- as.numeric(factor(constl$month))
  constl <- constl[-c(6,7)]
  
  code <- nimbleCode(
    {
      dist.beta[1] ~ dunif(0, 3000)
      dint.beta[1] <- logit(p.dint.beta)
      p.dint.beta ~ dunif(0,1)
      for (j in 2:8){ dint.beta[j] ~ dnorm(0, sd=100) }
      for (z in 2:8){ dist.beta[z] ~ dnorm(0, sd=100) }
      turb.sigma ~ dunif(0,20)
      month.sigma ~ dunif(0,20)
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
            dist.beta[2]*Sn[i,t-1] + dist.beta[3]*RotSp[i,t-1] + dist.beta[4]*perp[i,t-1] + dist.beta[5]*app[i,t-1] +
            dist.beta[6]*Sn[i,t-1]*RotSp[i,t-1] + 
            dist.beta[7]*Sn[i,t-1]*perp[i,t-1] + 
            dist.beta[8]*RotSp[i,t-1]*app[i,t-1]
          logit(dint[i,t-1]) <- dint.beta[1] + eps.turb[ turbine[i,t-1] ] + eps.month[ month[i] ] +
            dint.beta[2]*Sn[i,t-1] + dint.beta[3]*RotSp[i,t-1] + dint.beta[4]*app[i,t-1] +
            dint.beta[5]*Sn[i,t-1]*RotSp[i,t-1] + 
            dint.beta[6]*RotSp[i,t-1]*app[i,t-1] + 
            dint.beta[7]*ht[i,t-1] + dint.beta[8]*ht2[i,t-1] 
        } # time t
      } # track i
      for (k in 1:nturb){ eps.turb[k]~dnorm(0, sd=turb.sigma) } # ntracks
      for (m in 1:nmonth){ eps.month[m]~dnorm(0, sd=month.sigma) }
    } # code 
  ) 
  
  z.init <- datl$z
  z.init[is.na(z.init)] <- as.numeric(rbinom(prob=0.5, size=1, n=sum(is.na(z.init)) )) 
  z.init[!is.na(datl$z)] <- NA
  
  ##############################
  ## end model code
  ###############################
  inits <- function(){ list(
    z=z.init,
    p.dint.beta = runif(1),
    dist.beta = c(runif(1, 50, 500), runif(7, -0.2, 0.2)),
    dint.beta = runif(8, -0.2, 0.2),
    eps.turb = runif(constl$nturb, -0.1, 0.1),
    eps.month = runif(constl$nmonth, -0.1, 0.1),
    month.sigma = runif(1),
    phi = runif(1),
    psi1 = runif(1),
    turb.sigma = runif(1),
    dist.sigma = array(runif(constl$ntracks*(constl$ntime-1), 100, 1000) , dim=c(constl$ntracks, (constl$ntime-1))),
    dint = array(runif(constl$ntracks*(constl$ntime-1)) , dim=c(constl$ntracks, (constl$ntime-1))),
    Ez = array(runif(constl$ntracks*(constl$ntime-1)) , dim=c(constl$ntracks, (constl$ntime-1))),
    gamma = array(runif(constl$ntracks*(constl$ntime-1)) , dim=c(constl$ntracks, constl$ntime-1)),
    g= array(runif(constl$ntracks*(constl$ntime-1)) , dim=c(constl$ntracks, (constl$ntime-1)))
  ) }
  
  params <- c("p.dint.beta", "dint.beta", "dist.beta",
              "phi", "psi1",
              "turb.sigma", "month.sigma",
              "eps.turb", "eps.month")
  
  ni <- 50000; nc <- 3;  nb <- 25000; nt <- 25
  # ni <- 2000; nc <- 3;  nb <- 500; nt <- 1; nt2 <- 1
  
  mod<- nimbleModel(code, calculate=T, constants = constl, 
                    data = datl, inits = inits())
  cmod <- compileNimble(mod)
  
  mc <- configureMCMC(cmod, nchains=nc,
                      monitors= params, thin=nt, nburnin = nb)
  #monitors2=params2, thin2=nt2)
  mcmc <- buildMCMC(mc)
  cmcmc <- compileNimble(mcmc, project=cmod)
  
  try(samps <- runMCMC(  mcmc = cmcmc,
                         nchains=nc,
                         thin = nt,
                         #thin2 = nt2,
                         niter = ni,
                         nburnin = nb,
                         progressBar=T,
                         summary=F,
                         WAIC=F,
                         samplesAsCodaMCMC = T,
                         samples=T ))
  
  flnm2 <- paste("./results_Bayes_10s_cv_", xx, ".Rdata", sep="")       
  save(samps, file=flnm2)
  rm(list=ls())
} # xx