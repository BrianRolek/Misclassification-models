## ---- basic --------
m <- "misclass_basic"
library (nimble)
load(".\\data\\final-data.Rdata")

code <- nimbleCode(
  {
    ## PRIORS 
    mean.b ~ dunif(0, 1)
    b.b[1] <- logit(mean.b)
    mean.p10 ~dunif(0, 1)
    p10.b[1] <- logit(mean.p10)
    mean.p11 ~ dunif(0, 1)
    p11.b[1] <- logit(mean.p11)
    mean.phi ~ dunif(0, 1)
    phi.alpha[1] <- logit(mean.phi)
    mean.gamma ~ dunif(0, 1)
    gam.alpha[1] <- logit(mean.gamma)
    mean.psi ~ dunif(0, 1)
    psi.b <- logit(mean.psi)
    
    ## LIKELIHOOD
    ## first year
    for(i in 1:nsite){
      logit(psi1[i,1]) <- psi.b
      z[i,1] ~ dbern(psi1[i,1])
      for (j in 1:nvisit){
        pi[i,j,1,1] <- z[i,1]*(1-p11[i,j,1]) + (1-z[i,1])*(1-p10[i,j,1])      
        pi[i,j,1,2] <- z[i,1]*(1-b[i,j,1])*p11[i,j,1] + (1-z[i,1])*p10[i,j,1]
        pi[i,j,1,3] <- z[i,1]*b[i,j,1]*p11[i,j,1]   
        Y[i,j,1] ~ dcat(pi[i,j,1,1:3])
      } #j
      
      # subsequent years
      for (t in 2:nyear){
        for (j in 1:nvisit){ 
          pi[i,j,t,1] <- z[i,t]*(1-p11[i,j,t]) + (1-z[i,t])*(1-p10[i,j,t])      
          pi[i,j,t,2] <- z[i,t]*(1-b[i,j,t])*p11[i,j,t] + (1-z[i,t])*p10[i,j,t]
          pi[i,j,t,3] <- z[i,t]*b[i,j,t]*p11[i,j,t]
          Y[i,j,t] ~ dcat(pi[i,j,t,1:3])
        } # j
        # dynamics
        z[i,t] ~ dbern(muZ[i,t])
        muZ[i,t]<- z[i,t-1]*phi[i,t-1] + (1-z[i,t-1])*gamma[i,t-1]
        logit(phi[i,t-1]) <- phi.alpha[1]
        logit(gamma[i,t-1]) <- gam.alpha[1]
      } # t nyear 
      
      for (t in 1:nyear){
        for (j in 1:nvisit){
          # detection models
          logit(b[i,j,t]) <- b.b[1] 
          logit(p11[i,j,t]) <- p11.b[1] 
          logit(p10[i,j,t]) <- p10.b[1]  
        } } }# t j i
    
    # Predicted values for certainty, detection, and misclassification
    psi[1] <- mean.psi
    n.occ[1] <- sum(z[1:nsite,1])
    for (t in 2:nyear){
      n.occ[t] <- sum(z[1:nsite,t])
      logit(phi.est[t-1]) <- phi.alpha[1] 
      logit(gam.est[t-1]) <- gam.alpha[1] 
      psi[t] <- psi[t-1]*phi.est[t-1] + (1-psi[t-1])*gam.est[t-1]
      growthr[t-1] <- psi[t]/psi[t-1] 
      turnover[t-1] <- (1-psi[t-1]) * gam.est[t-1]/psi[t]
    } # t
  }
)

datl <- list(
  Y=dat.conv$Y,
  nsite=dim(dat.conv$Y)[[1]],
  nvisit=dim(dat.conv$Y)[[2]],
  nyear=dim(dat.conv$Y)[[3]]
)

params<-c("mean.p11", "p11.b", 
          "mean.b", "b.b",  
          "mean.p10", "p10.b", 
          "psi.b", 
          "mean.phi", "phi.alpha", "phi.est",
          "mean.gamma",  "gam.alpha", "gam.est",
          "psi", "n.occ", "growthr", "turnover",
          "z")

# Inits
maxstates <- apply(datl$Y, c(1, 3), max, na.rm=T )
z.inits <- ifelse(maxstates>1, 1, 0)
z.inits[is.na(z.inits)] <- 0
Y.inits <- array(NA, dim=dim(dat.conv$Y))
Y.inits[is.na(dat.conv$Y)] <- 1
occ1 <- apply(datl$Y, c(1,3), sum)

inits <- function()list ( 
  Y = Y.inits,
  z = z.inits,
  muZ= ifelse(occ1>0, 0.1, 0.8), 
  p10 = array(runif(datl$nsite*datl$nvisit*datl$nyear, 0.001, 0.1), dim=c(datl$nsite,datl$nvisit,datl$nyear)),
  p11 = array(runif(datl$nsite*datl$nvisit*datl$nyear, 0.3, 0.7), dim=c(datl$nsite,datl$nvisit,datl$nyear)),
  b = array(runif(datl$nsite*datl$nvisit*datl$nyear, 0.5, 0.99), dim=c(datl$nsite,datl$nvisit,datl$nyear)),
  mean.psi=runif(1),
  mean.p11=runif(1),
  mean.p10=runif(1),
  mean.b=runif(1),
  psi.b= runif(1, -5, 5),
  b.b= runif(1, -5, 5),
  p11.b= runif(1, -5, 5),
  p10.b= runif(1, -5, 5),
  mean.gamma= runif(1),
  mean.phi=runif(1),
  phi.alpha= runif(1, -5, 5),
  gam.alpha= runif(1, -5, 5)
) 
n.chains=3; n.thin=200; n.iter=600000; n.burnin=400000

mod <- list()
mod<- nimbleModel(code, calculate=T, constants = datl[-1], 
                  data = list(Y=datl$Y), inits = inits())

out <- nimbleMCMC(
  model = mod,
  code = code,
  monitors = params,
  nchains = n.chains,
  thin = n.thin,
  niter = n.iter,
  nburnin = n.burnin,
  progressBar = T,
  summary = T,
  WAIC = F,
  samplesAsCodaMCMC = T,
  samples=T
)

#flnm <- paste(".\\",m, "_", Sys.Date(), ".Rdata", sep="")
#save(out, mod, file=flnm)

