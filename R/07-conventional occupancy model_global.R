## GVS step 2
## ---- global model --------
m <- "07-conventional occupancy model_global"
library (nimble)
load(".\\data\\final-data.Rdata")

code <- nimbleCode(
  {
    ## PRIORS 
    mean.p11 ~ dunif(0, 1)
    p11.b[1] <- logit(mean.p11)
    mean.phi ~ dunif(0, 1)
    phi.alpha[1] <- logit(mean.phi)
    mean.gamma ~ dunif(0, 1)
    gam.alpha[1] <- logit(mean.gamma)
    mean.psi ~ dunif(0, 1)
    psi.b <- logit(mean.psi)
    
    for (aa in 2:6){
      phi.alpha[aa] ~ dnorm(0,10)
    }
    for (aa in 2:9){
      gamma.alpha[aa] ~ dnorm(0,10)
    }

    p11.b[2] ~ dnorm(0, sd=100) 
    for (t in 1:(nyear-1)){ 
      eps.phi[t] ~ dnorm(0, sd=sig.phi) 
      eps.gam[t] ~ dnorm(0, sd=sig.gam) 
    } # t
    for (t in 1:nyear){ 
      eps.p11[t] ~ dnorm(0, sd=sig.p11)
    }
      sig.phi ~ T(dnorm(0,sd=10),0, )
      sig.gam ~ T(dnorm(0,sd=10),0, )
      sig.p11 ~ T(dnorm(0,sd=10),0, )
    # priors for the w model inclusion terms, gam.alpha
    # this ensures that each of the 8 model combos has equal probability: Pr(m)= 1/8

    ## LIKELIHOOD
    ## first year
    for(i in 1:nsite){
      logit(psi1[i,1]) <- psi.b
      z[i,1] ~ dbern(psi1[i,1])
        for (t in 2:nyear){
          z[i,t] ~ dbern(muZ[i,t])
          muZ[i,t]<- z[i,t-1]*phi[i,t-1] + (1-z[i,t-1])*gamma[i,t-1]
         } } # i t
      
    for(i in 1:nsite){
      for (t in 1:(nyear-1)){
        logit(phi[i,t]) <- 
          phi.alpha[1] + 
          phi.alpha[2] * YSF.std[i,t, 6 ] + 
          phi.alpha[3] * sin(SEAS[i,t, 1 ]*2*3.1416) + 
          phi.alpha[4] * cos(SEAS[i,t, 1 ]*2*3.1416) +
          phi.alpha[5] * sin(SEAS[i,t, 1 ]*2*3.1416) * YSF.std[i,t, 6 ] + 
          phi.alpha[6] * cos(SEAS[i,t, 1 ]*2*3.1416) * YSF.std[i,t, 6 ] +
          eps.phi[t]
        logit(gamma[i,t]) <- 
          gam.alpha[1] + 
          gam.alpha[2] * YSF.std[i,t, 5 ] +gam.alpha[3] * YSF.std[i,t, 5 ]^2 +
          gam.alpha[4] * sin(SEAS[i,t, 1 ]*2*3.1416) + 
          gam.alpha[5] * cos(SEAS[i,t, 1 ]*2*3.1416) +
          gam.alpha[6] * sin(SEAS[i,t, 1 ]*2*3.1416) * YSF.std[i,t, 5 ] +
          gam.alpha[7] * cos(SEAS[i,t, 1 ]*2*3.1416) * YSF.std[i,t, 5 ] +
          gam.alpha[8] * sin(SEAS[i,t, 1 ]*2*3.1416) * YSF.std[i,t, 5 ] * YSF.std[i,t, 5 ]^2 +
          gam.alpha[9] * cos(SEAS[i,t, 1 ]*2*3.1416) * YSF.std[i,t, 5 ] * YSF.std[i,t, 5 ]^2 +
          eps.gam[t]
      }} # t nyear 
      
      for (i in 1:nsite){
        for (j in 1:nvisit){
          for (t in 1:nyear){
          # detection models
          logit(p11[i,j,t]) <- p11.b[1] + p11.b[2]*hr[i,j,t] + eps.p11[t]
          Y[i,j,t] ~ dbern(z[i,t] * p11[i,j,t] )
          } } }# t j i
    
    # Predicted occupancy values 
    psi[1] <- mean.psi
    n.occ[1] <- sum(z[1:nsite,1])
    for (t in 2:nyear){
      n.occ[t] <- sum(z[1:nsite,t])
      phi.est[t-1] <- mean(phi[1:nsite,t-1])
      gam.est[t-1] <- mean(gamma[1:nsite,t-1])
      psi[t] <- psi[t-1]*phi.est[t-1] + (1-psi[t-1])*gam.est[t-1]
      growthr[t-1] <- psi[t]/psi[t-1] 
      turnover[t-1] <- (1-psi[t-1]) * gam.est[t-1]/psi[t]
    } # t
  }
)

#  scale and center fire data
YSF.std <- array(NA, dim=dim(YSF))
for (j in 1:6){
  YSF.std[,,j] <- (YSF[,,j]-mean(YSF[,,j], na.rm=T)) / sd(YSF[,,j])
}

datl <- list(
  Y=ifelse(dat.conv$Y>1, 1, 0),
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

params <- c("mean.p11", "p11.b", "eps.p11", "sig.p11",
          "psi.b", 
          "mean.phi", "phi.alpha", "eps.phi", "sig.phi", 
          "mean.gamma",  "gam.alpha", "eps.gam", "sig.gam", 
          "psi", "n.occ", 
          "z")

# Inits
maxstates <- apply(datl$Y, c(1, 3), max, na.rm=T )
z.inits <- ifelse(maxstates>0, 1, 0)
z.inits[is.na(z.inits)] <- 0
Y.inits <- datl$Y
Y.inits[is.na(dat.conv$Y)] <- 0

inits <- function()list ( 
  Y = Y.inits,
  z = z.inits,
  #muZ= ifelse(z.inits>0, 0.1, 0.8),
  p11 = array(runif(datl$nsite*datl$nvisit*datl$nyear, 0.3, 0.7), dim=c(datl$nsite,datl$nvisit,datl$nyear)),
  mean.psi=runif(1),
  mean.p11=runif(1),
  p11.b=runif(2, -5, 5),
  mean.gamma= runif(1),
  mean.phi=runif(1),
  phi.alpha= runif(6, -5, 5),
  gam.alpha= runif(9, -5, 5),
  sig.p11=runif(1),
  sig.phi=runif(1),
  sig.gam=runif(1),
  eps.phi = runif((datl$nyear-1), -0.1, 0.1),
  eps.gam = runif((datl$nyear-1), -0.1, 0.1),
  eps.p11 = runif((datl$nyear), -0.1, 0.1),
  gamma= array(runif(datl$nsite*(datl$nyear-1), 0, 1), dim=c(datl$nsite,datl$nyear-1)),
  phi= array(runif(datl$nsite*(datl$nyear-1), 0, 1), dim=c(datl$nsite,datl$nyear-1)),
  psi= runif(datl$nyear, 0, 1)
) 
n.chains=3; n.thin=50; n.iter=100000; n.burnin=50000
#n.chains=3; n.thin=1; n.iter=5000; n.burnin=100 # trial runs

mod <- nimbleModel(code, calculate=T, constants = datl[-1], 
                  data = list(Y=datl$Y), inits = inits())

#mod$simulate("Y")

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

flnm <- paste(".\\outputs\\",m, "_", Sys.Date(), ".Rdata", sep="")
save(out, mod, file=flnm)

flnm <- paste("C:\\Users\\rolek.brian\\Documents\\Projects\\FGSP Misclassification\\Results\\",m, "_", Sys.Date(), ".Rdata", sep="")
save(outc=out, mod, file=flnm)
