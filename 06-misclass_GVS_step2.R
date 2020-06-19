m <- "06-misclass_GVS_step2"
library (nimble)
load("/scratch/brolek/fgsp_misclass/data/final-data.Rdata")
load("/scratch/brolek/fgsp_misclass/outputs/05-misclass_GVS_step1_2020-06-15.Rdata")
outg <- out
rm(list="out")

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
    sig.p10 ~ T(dnorm(0,sd=10),0, )
    sig.p11 ~ T(dnorm(0,sd=10),0, )
    bp.sig <- 100
    bg.sig <- 100
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
    
    # priors for the w model inclusion terms, phi.alpha
    # this ensures that each of the 8 model combos has equal probability: Pr(m)= 1/8
    wp[7] ~ dbern(0.5)
    wp[6] ~ dbern(0.23)
    wp[5] ~ dbern(0.23) 
    p.wp4 <- (1-wp[6])*0.5 + wp[6] 
    wp[4] ~ dbern(p.wp4)
    p.wp3 <- (1-wp[5])*0.5 + wp[5] 
    wp[3] ~ dbern(p.wp3)
    p.wp2 <- equals(wp[5]+wp[6],0)*0.5 + (1-equals(wp[5]+wp[6],0)) 
    wp[2] ~ dbern(p.wp2)
    wp[1] ~ dbern(1)
    wptemp[1] <- wp[1]
    
    # priors for the w model inclusion terms, gam.alpha
    # this ensures that each of the 8 model combos has equal probability: Pr(m)= 1/8
    wg[7] ~ dbern(0.5)
    wg[6] ~ dbern(0.23)
    wg[5] ~ dbern(0.23) 
    p.wg4 <- (1-wg[6])*0.5 + wg[6] 
    wg[4] ~ dbern(p.wg4)
    p.wg3 <- (1-wg[5])*0.5 + wg[5] 
    wg[3] ~ dbern(p.wg3)
    p.wg2 <- equals(wg[5]+wg[6],0)*0.5 + (1-equals(wg[5]+wg[6],0)) 
    wg[2] ~ dbern(p.wg2)
    wg[1] ~ dbern(1)
    wgtemp[1] <- wg[1]
    
    # set up the vectors/matrices for beta estimation, abundance
    for(b1 in 2:n.betasp){
      wptemp[b1] <- wp[posp[b1]]                # this uses GVS
      # wptemp[b1] <- 1                          # this forces you to fit the full model (all terms included)
      mean.bp[b1] <- post.bp[b1]*(1-wptemp[b1])  # prior is either 0 or full-model posterior mean
      for(b2 in 2:n.betasp){                   # set up the precision matrix (inverse variance) # allows for betas to be multivariate, if desired
        sig.bp[b1,b2] <- equals(b1,b2)*sd.bp[b1]*(1-wptemp[b1]) + (wptemp[b1]*bp.sig)
      } # b2
      phi.alpha[b1] ~ T(dnorm(mean.bp[b1], sd=sig.bp[b1,b1]), -30, 30)   # all beta coefficients
    } # b1
    #wptemp[7] <- 1 
    wptemp[7] <- wp[7]
    mean.sigphi <- post.bp[7]*(1-wptemp[7]) 
    sd.sigphi <- sd.bp[7]*(1-wptemp[7]) + (wptemp[7]*bp.sigphi)
    sig.phi ~ T(dnorm(mean.sigphi, sd=sd.sigphi), 0, )
    bp.sigphi <- 10
    # set up the vectors/matrices for beta estimation, abundance
    for(b1 in 2:n.betasg){
      wgtemp[b1] <- wg[posg[b1]]                # this uses GVS
      # wgtemp[b1] <- 1                          # this forces you to fit the full model (all terms included)
      mean.bg[b1] <- post.bg[b1]*(1-wgtemp[b1])  # prior is either 0 or full-model posterior mean
      for(b2 in 2:n.betasg){                   # set up the precision matrix (inverse variance) # allows for betas to be multivariate, if desired
        sig.bg[b1,b2] <- equals(b1,b2)*sd.bg[b1]*(1-wgtemp[b1]) + (wgtemp[b1]*bg.sig)
      } # b2
      gam.alpha[b1] ~ T(dnorm(mean.bg[b1], sd=sig.bg[b1,b1]), -30, 30)   # all beta coefficients
    } # b1
    #wgtemp[7] <- 1 
    wgtemp[8] <- wg[7]
    mean.siggam <- post.bg[8]*(1-wgtemp[8]) 
    sd.siggam <- sd.bg[8]*(1-wgtemp[8]) + (wgtemp[8]*bg.siggam)
    sig.gam ~ T(dnorm(mean.siggam, sd=sd.siggam), 0, )
    bg.siggam <- 10
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
        logit(phi[i,t-1]) <- 
          wptemp[1] * phi.alpha[1] + 
          wptemp[2] * phi.alpha[2] * YSF.std[i,t, 6 ] + 
          wptemp[3] * phi.alpha[3] * sin(SEAS[i,t, 2 ]*2*3.1416) + 
          wptemp[4] * phi.alpha[4] * cos(SEAS[i,t, 2 ]*2*3.1416) +
          wptemp[5] * phi.alpha[5] * sin(SEAS[i,t, 2 ]*2*3.1416) * YSF.std[i,t, 6 ] + 
          wptemp[6] * phi.alpha[6] * cos(SEAS[i,t, 2 ]*2*3.1416) * YSF.std[i,t, 6 ] +
          wptemp[7] * eps.phi[t-1]
        logit(gamma[i,t-1]) <- 
          wgtemp[1] * gam.alpha[1] + 
          wgtemp[2] * gam.alpha[2] * YSF.std[i,t, 6 ] + wgtemp[3] * gam.alpha[3] * YSF.std[i,t, 6 ]^2 +
          wgtemp[4] * gam.alpha[4] * sin(SEAS[i,t, 1 ]*2*3.1416) + 
          wgtemp[5] * gam.alpha[5] * cos(SEAS[i,t, 1 ]*2*3.1416) +
          wgtemp[6] * gam.alpha[6] * sin(SEAS[i,t, 1 ]*2*3.1416) * YSF.std[i,t, 6 ] * YSF.std[i,t, 6 ]^2 +
          wgtemp[7] * gam.alpha[7] * cos(SEAS[i,t, 1 ]*2*3.1416) * YSF.std[i,t, 6 ] * YSF.std[i,t, 6 ]^2 +
          wgtemp[8] * eps.gam[t-1]
        } # t nyear 
      
      for (t in 1:nyear){
        for (j in 1:nvisit){
          # detection models
          logit(b[i,j,t]) <-  b.b[1] + b.b[2]*date[i,j,t] + b.b[3]*date[i,j,t]^2+ b.b[4]*date[i,j,t]^3 + eps.b[t]
          logit(p11[i,j,t]) <- p11.b[1] + p11.b[2]*hr[i,j,t] + eps.p11[t]
          logit(p10[i,j,t]) <- p10.b[1]  + p10.b[2]*date[i,j,t] + p10.b[3]*date[i,j,t]^2 + eps.p10[t]
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

# function to extract samples from nimble output  
extr <- function (mod, var){
  var2 <- paste("^", var, sep="")
  col.ind <- grep(var2, colnames(outg$samples$chain1))
  mat1 <- as.matrix(outg$samples$chain1[,col.ind])
  mat2 <- as.matrix(outg$samples$chain2[,col.ind])
  mat3 <- as.matrix(outg$samples$chain3[,col.ind])
  all <- rbind(mat1, mat2, mat3) 
  print(colnames(all))
  return(all)
}

# extract samples
pa <- extr(outg, "phi.alpha")
pa.sig <- extr(outg, "sig.phi")
pg <- extr(outg, "gam.alpha")
pg.sig <- extr(outg, "sig.gam")

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
  nSEAS= dim(SEAS)[[3]], 
  n.betasp = 6,
  n.betasg = 6,
  posp = c(1, 2, 3, 4, 5, 6), 
  posg = c(1, 2, 2, 3, 4, 5, 6),
  post.bp = c(apply(pa, 2, mean), mean(pa.sig)),
  sd.bp = c(apply(pa, 2, sd), sd(pa.sig)),
  post.bg = c(apply(pg, 2, mean), mean(pg.sig)),
  sd.bg = c(apply(pg, 2, sd), sd(pg.sig)) 
)

params<-c("mean.p11", "p11.b", "eps.p11", "sig.p11",
          "mean.b", "b.b", "eps.b", "sig.b",
          "mean.p10", "p10.b", "eps.p10", "sig.p10",
          "psi.b", 
          "mean.phi", "phi.alpha", "eps.phi", "sig.phi", "phi.est",
          "mean.gamma",  "gam.alpha", "eps.gam", "sig.gam", "gam.est",
          "psi", "n.occ", "growthr", "turnover", 
          "wg", "wp", "wgtemp", "wptemp",
          "z")

# Inits
maxstates <- apply(datl$Y, c(1, 3), max, na.rm=T )
z.inits <- ifelse(maxstates>1, 1, 0)
z.inits[is.na(z.inits)] <- 0
Y.inits <- array(NA, dim=dim(dat.conv$Y))
Y.inits[is.na(dat.conv$Y)] <- 1

inits <- function()list ( 
  Y = Y.inits,
  z = z.inits,
  muZ= ifelse(z.inits>0, 0.1, 0.8), 
  p10 = array(runif(datl$nsite*datl$nvisit*datl$nyear, 0.001, 0.1), dim=c(datl$nsite,datl$nvisit,datl$nyear)),
  p11 = array(runif(datl$nsite*datl$nvisit*datl$nyear, 0.3, 0.7), dim=c(datl$nsite,datl$nvisit,datl$nyear)),
  b = array(runif(datl$nsite*datl$nvisit*datl$nyear, 0.5, 0.99), dim=c(datl$nsite,datl$nvisit,datl$nyear)),
  mean.psi=runif(1),
  mean.p11=runif(1),
  mean.p10=runif(1, 0.01, 0.1),
  mean.b=runif(1),
  b.b=c(runif(4, -5, 5)),
  p10.b=c(runif(3, -5, 5)),
  p11.b=runif(2, -5, 5),
  mean.gamma= runif(1),
  mean.phi=runif(1),
  phi.alpha= runif(6, -5, 5),
  gam.alpha= runif(7, -5, 5),
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
  wp= rbinom(n=7, size=1, prob=1),
  wg= rbinom(n=7, size=1, prob=1),
  bp.sig=100,
  bg.sig=100,
  gamma= array(runif(datl$nsite*(datl$nyear-1), 0, 1), dim=c(datl$nsite,datl$nyear-1)),
  phi= array(runif(datl$nsite*(datl$nyear-1), 0, 1), dim=c(datl$nsite,datl$nyear-1)),
  psi= runif(datl$nyear, 0, 1),
  gam.est=runif(datl$nyear-1, 0, 1),
  phi.est=runif(datl$nyear-1, 0, 1),
  growthr=runif(datl$nyear-1, 0.95, 1.05),
  turnover=runif(datl$nyear-1, 0, 1), 
  wgtemp = rep(1, 8),
  wptemp = rep(1, 7)
) 
n.chains=6; n.thin=200; n.iter=600000; n.burnin=400000
#n.chains=3; n.thin=1; n.iter=5000; n.burnin=100 # trial runs

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
  
flnm <- paste("/scratch/brolek/fgsp_misclass/outputs/",m, "_", Sys.Date(), ".Rdata", sep="")
save(out, mod, file=flnm)

