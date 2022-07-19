## ---- BLISS --------
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
    psi.b[1] <- logit(mean.psi)
    for (xx in 2:4) { b.b[xx] ~ dnorm(0, sd=100) }
    for (xx in 2:3) { p10.b[xx] ~ dnorm(0, sd=100) }
    for (xx in 2) { p11.b[xx] ~ dnorm(0, sd=100) }
    
    for (t in 1:(nyear-1)){ 
      eps.phi[t] ~ dnorm(0, sd=sig.phi) 
      eps.gam[t] ~ dnorm(0, sd=sig.gam) 
    } # t
    for (t in 1:nyear){ 
      eps.p10[t] ~ dnorm(0, sd=sig.p10) 
      eps.p11[t] ~ dnorm(0, sd=sig.p11)
      eps.b[t] ~ dnorm(0, sd=sig.b)} # t 
    sig.p10 ~ T(dnorm(0,sd=10),0, )
    sig.p11 ~ T(dnorm(0,sd=10),0, )
    sig.b ~ T(dnorm(0,sd=10),0, )
    bliss[1] ~ dcat( priors1[1:12] )
    bliss[3] ~ dcat( priors3[1:12] )
    for (jj in 1:12){
    priors1[jj] <-  0.0833
    priors3[jj] <-  0.0833
    } # jj
    bliss[2] ~ dcat( priors2[1:2] )
    bliss[4] ~ dcat( priors4[1:2] )
    for (kk in 1:2){
      priors2[kk] <- 0.5
      priors4[kk] <- 0.5
    } # kk
    ## LIKELIHOOD
    ## first year
    for(i in 1:nsite){
      logit(psi1[i,1]) <- psi.b[1]
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
        logit(phi[i,t-1]) <- 
          phi.alpha0 + eps.phi[t-1] +
          equals(bliss[1],1) * phi.alpha[1] *  YSF.std[i,t, 1 ]  +  
          equals(bliss[1],1) * phi.alpha[2]*  YSF.std[i,t, 1 ]^2  +
          equals(bliss[1],2) * phi.alpha[1] *  YSF.std[i,t, 2 ]  +  
          equals(bliss[1],2) * phi.alpha[2]*  YSF.std[i,t, 2 ]^2  +
          equals(bliss[1],3) * phi.alpha[1] *  YSF.std[i,t, 3 ]  +  
          equals(bliss[1],3) * phi.alpha[2]*  YSF.std[i,t, 3 ]^2  +
          equals(bliss[1],4) * phi.alpha[1] *  YSF.std[i,t, 4 ]  +  
          equals(bliss[1],4) * phi.alpha[2]*  YSF.std[i,t, 4 ]^2  +
          equals(bliss[1],5) * phi.alpha[1] *  YSF.std[i,t, 5 ]  +  
          equals(bliss[1],5) * phi.alpha[2]*  YSF.std[i,t, 5 ]^2  +
          equals(bliss[1],6) * phi.alpha[1] *  YSF.std[i,t, 6 ]  +  
          equals(bliss[1],6) * phi.alpha[2]*  YSF.std[i,t, 6 ]^2 +
          equals(bliss[1],7) * phi.alpha[1] *  YSF.std[i,t, 1 ]  +  
          equals(bliss[1],8) * phi.alpha[1] *  YSF.std[i,t, 2 ]  +  
          equals(bliss[1],9) * phi.alpha[1] *  YSF.std[i,t, 3 ]  +  
          equals(bliss[1],10) * phi.alpha[1] *  YSF.std[i,t, 4 ]  +  
          equals(bliss[1],11) * phi.alpha[1] *  YSF.std[i,t, 5 ]  +  
          equals(bliss[1],12) * phi.alpha[1] *  YSF.std[i,t, 6 ]  +  
          equals(bliss[2],1) * phi.alpha[3] * sin(SEAS[i,t, 1 ]*2*3.1416) + equals(bliss[2],1) * phi.alpha[4] * cos(SEAS[i,t, 1 ]*2*3.1416) +
          equals(bliss[2],2) * phi.alpha[3] * sin(SEAS[i,t, 2 ]*2*3.1416) + equals(bliss[2],2) * phi.alpha[4] * cos(SEAS[i,t, 2 ]*2*3.1416)
        
        logit(gamma[i,t-1]) <- 
          gamma.alpha0 + eps.gam[t-1] +
          equals(bliss[3],1) * gam.alpha[1] *  YSF.std[i,t, 1 ]  +  
          equals(bliss[3],1) * gam.alpha[2]*  YSF.std[i,t, 1 ]^2  +
          equals(bliss[3],2) * gam.alpha[1] *  YSF.std[i,t, 2 ]  +  
          equals(bliss[3],2) * gam.alpha[2]*  YSF.std[i,t, 2 ]^2  +
          equals(bliss[3],3) * gam.alpha[1] *  YSF.std[i,t, 3 ]  +  
          equals(bliss[3],3) * gam.alpha[2]*  YSF.std[i,t, 3 ]^2  +
          equals(bliss[3],4) * gam.alpha[1] *  YSF.std[i,t, 4 ]  +  
          equals(bliss[3],4) * gam.alpha[2]*  YSF.std[i,t, 4 ]^2  +
          equals(bliss[3],5) * gam.alpha[1] *  YSF.std[i,t, 5 ]  +  
          equals(bliss[3],5) * gam.alpha[2]*  YSF.std[i,t, 5 ]^2  +
          equals(bliss[3],6) * gam.alpha[1] *  YSF.std[i,t, 6 ]  +  
          equals(bliss[3],6) * gam.alpha[2]*  YSF.std[i,t, 6 ]^2 +
          equals(bliss[3],7) * gam.alpha[1] *  YSF.std[i,t, 1 ]  +  
          equals(bliss[3],8) * gam.alpha[1] *  YSF.std[i,t, 2 ]  +  
          equals(bliss[3],9) * gam.alpha[1] *  YSF.std[i,t, 3 ]  +  
          equals(bliss[3],10) * gam.alpha[1] *  YSF.std[i,t, 4 ]  +  
          equals(bliss[3],11) * gam.alpha[1] *  YSF.std[i,t, 5 ]  +  
          equals(bliss[3],12) * gam.alpha[1] *  YSF.std[i,t, 6 ]  +  
          equals(bliss[4],1) * gam.alpha[3] * sin(SEAS[i,t, 1 ]*2*3.1416) + equals(bliss[4],1) * gam.alpha[4] * cos(SEAS[i,t, 1 ]*2*3.1416) +
          equals(bliss[4],2) * gam.alpha[3] * sin(SEAS[i,t, 2 ]*2*3.1416) + equals(bliss[4],2) * gam.alpha[4] * cos(SEAS[i,t, 2 ]*2*3.1416)
         } # t nyear 
      
        for (t in 1:nyear){
          for (j in 1:nvisit){
          # detection models
          logit(b[i,j,t]) <- b.b[1] + b.b[2]*date[i,j,t] + b.b[3]*date[i,j,t]*2+ b.b[4]*date[i,j,t]^3 + eps.b[t]
          logit(p11[i,j,t]) <- p11.b[1] + p11.b[2]*hr[i,j,t] + eps.p11[t]
          logit(p10[i,j,t]) <- p10.b[1]  + p10.b[2]*date[i,j,t] + p10.b[3]*date[i,j,t]^2 + eps.p10[t]
        } } }# t j i

    psi[1] <- mean.psi
    n.occ[1] <- sum(z[1:nsite,1])
    for (t in 2:nyear){
      n.occ[t] <- sum(z[1:nsite,t])
      logit(phi.est[t-1]) <- phi.alpha0 + eps.phi[t-1]
      logit(gam.est[t-1]) <- gamma.alpha0 + eps.gam[ t-1]
      psi[t] <- psi[t-1]*phi.est[t-1] + (1-psi[t-1])*gam.est[t-1]
      growthr[t-1] <- psi[t]/psi[t-1] 
      turnover[t-1] <- (1-psi[t-1]) * gam.est[t-1]/psi[t]
      } # t
  }
)

YSF.std <- YSF2.std <- array(NA, dim=dim(YSF))
for (i in 1:6){
  YSF.std[,,i] <- (YSF[,,i]-mean(YSF[,,i], na.rm=T)) / sd(YSF[,,i])
}

datl <- list(
  Y=dat.conv$Y,
  date= dat.conv$date,
  hr= dat.conv$hr,
  nsite=dim(dat.conv$Y)[[1]],
  nvisit=dim(dat.conv$Y)[[2]],
  nyear=dim(dat.conv$Y)[[3]],
  YSF.std=YSF.std,
  SEAS=SEAS,
  nYSF= dim(YSF)[[3]],
  nSEAS= dim(SEAS)[[3]]
)

params<-c("mean.p11", "p11.b", "eps.p11", "sig.p11",
          "mean.b", "b.b", "eps.b", "sig.b",  
          "mean.p10", "p10.b", "eps.p10", "sig.p10",
          "psi.b", 
          "mean.phi", "phi.alpha0", "phi.alpha", "eps.phi", "sig.phi", "phi.est",
          "mean.gamma",  "gamma.alpha0", "gam.alpha", "eps.gam", "sig.gam", "gam.est",
          "psi", "n.occ", "growthr", "turnover", 
          "bliss",
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
  mean.p10=runif(1, 0.01, 0.1),
  mean.b=runif(1),
  psi.b= runif(1, -5, 5),
  b.b= runif(4, -5, 5),
  p11.b= runif(2, -5, 5),
  p10.b= runif(3, -5, 5),
  mean.gamma= runif(1),
  mean.phi=runif(1),
  phi.alpha0= runif(1, -5, 0),
  gamma.alpha0= runif(1, 0, 5),
  phi.alpha= runif(4, -2, 2),
  gam.alpha= runif(4, -2, 2),
  sig.p10=runif(1),
  sig.p11=runif(1),
  sig.b=runif(1),
  sig.phi=runif(1),
  sig.gam=runif(1),
  eps.phi = runif((datl$nyear-1), -0.1, 0.1),
  eps.gam = runif((datl$nyear-1), -0.1, 0.1),
  eps.p10 = runif((datl$nyear), -0.1, 0.1),
  eps.p11 = runif((datl$nyear), -0.1, 0.1),
  eps.b = runif((datl$nyear), -0.1, 0.1),
  bliss = c(sample(1:12,1), sample(c(1,2),1), 
           sample(1:12,1), sample(c(1,2),1))
) 

n.chains=3; n.thin=200; n.iter=600000; n.burnin=400000
mod<- nimbleModel(code, calculate=T, constants = datl[-1], 
                  data = list(Y=datl$Y), inits = inits())

out <- nimbleMCMC(
  model=mod,
  code = code,
  monitors = params,
  nchains=n.chains,
  thin = n.thin,
  niter = n.iter,
  nburnin = n.burnin,
  progressBar=T,
  summary=T,
  WAIC=F,
  samplesAsCodaMCMC = T,
  samples=T
)

# flnm <- paste(".\\outputs\\misclass-BLISS_", Sys.Date(), ".Rdata", sep="")
# save(out, file=flnm)

