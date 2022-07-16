## ---- prior probs --------
library(MuMIn)
load( ".\\data\\final-data.Rdata")

vars <- c("eps.phi", "YSF.std","sin.SEAS", "cos.SEAS")
tf <- c(TRUE, FALSE)
eg <- expand.grid(tf, tf, tf, tf)
colnames(eg) <- vars

test1 <- model.matrix(~ YSF.std + sin.SEAS + cos.SEAS, eg)

dat <- matrix(rnorm(100*length(vars)),ncol=length(vars))
dat <- as.data.frame(dat)
names(dat) <- vars
dat$y <- rnorm(nrow(dat))

global <- lm(y ~ YSF.std + sin.SEAS + cos.SEAS +
               YSF.std:sin.SEAS + YSF.std:cos.SEAS, 
             data = dat, na.action = "na.fail")

# dummy dat

combos <- dredge(global, eval=F,
                 subset = 
                   dc(YSF.std, YSF.std:sin.SEAS) &&
                   dc(YSF.std, YSF.std:cos.SEAS) &&
                   dc(sin.SEAS, YSF.std:sin.SEAS) &&
                   dc(cos.SEAS, YSF.std:cos.SEAS) 
)

n.mods <- length(combos)
dat1 <- dat[1,]*0 + 1
terms <- colnames(model.matrix(formula(combos[[n.mods]]),dat1))
n.terms <- length(terms)

combo.mat <- matrix(0,nrow=n.mods,ncol=n.terms)
colnames(combo.mat) <- terms

for (i in 1:n.mods){
  combo.mat[i,match(colnames(model.matrix(formula(combos[[i]]),dat1)),
                    terms)] <- 1
}

combo.mat

# total model probs
apply(combo.mat,2,mean)

# highest order terms
mean(combo.mat[,"cos.SEAS:YSF.std"])
mean(combo.mat[,"sin.SEAS:YSF.std"])

# probabilities without higher order terms
# ba:sf
mean(combo.mat[which(combo.mat[,"cos.SEAS:YSF.std"]==0 & combo.mat[,"sin.SEAS:YSF.std"]==0),"YSF.std"])

##################
# COLONIZATION
#################
global2 <- lm(y ~ YSF.std + I(YSF.std^2) + 
                sin.SEAS + cos.SEAS +
                YSF.std:sin.SEAS + YSF.std:cos.SEAS + 
                I(YSF.std^2):sin.SEAS + I(YSF.std^2):cos.SEAS,
              data = dat, na.action = "na.fail")

combos2 <- dredge(global2, eval=F,
                  subset = 
                    dc(YSF.std, sin.SEAS:YSF.std, {sin.SEAS:I(YSF.std^2)}) &&
                    dc(YSF.std, cos.SEAS:YSF.std, {cos.SEAS:I(YSF.std^2)}) &&
                    dc(sin.SEAS, sin.SEAS:YSF.std) &&
                    dc(cos.SEAS, cos.SEAS:YSF.std) &&
                    dc(YSF.std, {I(YSF.std^2)}, {sin.SEAS:I(YSF.std^2)}) &&
                    dc(YSF.std, {I(YSF.std^2)}, {cos.SEAS:I(YSF.std^2)}) &&
                    dc({I(YSF.std^2)}, YSF.std ) &&
                    dc({sin.SEAS:I(YSF.std^2)}, sin.SEAS:YSF.std) #&&
                    # dc({cos.SEAS:I(YSF.std^2)}, cos.SEAS:YSF.std)
                  )

n.mods2 <- length(combos2)
dat2 <- dat[1,]*0 + 1
terms2 <- colnames(model.matrix(formula(combos2[[n.mods2]]),dat2))
n.terms2 <- length(terms2)

combo.mat2 <- matrix(0, nrow=n.mods2, ncol=n.terms2)
colnames(combo.mat2) <- terms2

for (i in 1:n.mods2){
  combo.mat2[i,match(colnames(model.matrix(formula(combos2[[i]]),dat2)),
                    terms2)] <- 1
}

combo.mat2

# total model probs
sub2 <- apply(combo.mat2[,6:9]==0, 1, all)
apply(combo.mat2[sub2,],2,mean)[1:5]
sub3 <- apply(combo.mat2[,8:9]==0, 1, all)
apply(combo.mat2[sub3,],2,mean)[6:7]
apply(combo.mat2,2,mean)[8:9]


vars <- c("eps.gam", "YSF","sin.SEAS", "cos.SEAS",
          "YSF2", "YSF:sin.SEAS", "YSF:cos.SEAS",
          "YSF2:sin.SEAS", "YSF2:cos.SEAS")
tf <- c(TRUE, FALSE)
tfl <- list()
for (i in 1:length(vars)){ tfl[[i]] <- tf }
eg <- expand.grid( tfl ) 
colnames(eg) <- vars
# subset
ind <- 
  !(eg$YSF==F & eg$YSF2==T) & !(eg$YSF==T & eg$YSF2==F) &
  !(eg$YSF==F & (eg$'YSF:sin.SEAS'==T | eg$'YSF2:sin.SEAS'==T)) &
  !(eg$YSF==F & (eg$'YSF:cos.SEAS'==T | eg$'YSF2:cos.SEAS'==T)) &
  !(eg$sin.SEAS==F & (eg$'YSF:sin.SEAS'==T | eg$'YSF2:sin.SEAS'==T)) &
  !(eg$cos.SEAS==F & (eg$'YSF:cos.SEAS'==T | eg$'YSF2:cos.SEAS'==T)) &
  !(eg$'YSF:sin.SEAS'==F & eg$'YSF2:sin.SEAS'==T) & !(eg$'YSF:sin.SEAS'==T & eg$'YSF2:sin.SEAS'==F) &
  !(eg$'YSF:cos.SEAS'==F & eg$'YSF2:cos.SEAS'==T) & !(eg$'YSF:cos.SEAS'==T & eg$'YSF2:cos.SEAS'==F)
eg <- eg[ind, ]

wg.priors <- list()
wg.priors[[1]] <- 1 
sub1 <- !apply(eg[,6:9], 1, any)
wg.priors[[2]] <-  apply(eg[sub1, c("YSF", "YSF2")], 2, mean)
sub2 <- !apply(eg[,c(6,8)], 1, any)
wg.priors[[3]] <- mean(eg[sub2, c("sin.SEAS")])
sub3 <- !apply(eg[,c(7,9)], 1, any)
wg.priors[[4]] <- mean(eg[sub3, c("cos.SEAS")])
wg.priors[[5]] <- apply(eg[, c("YSF:sin.SEAS", "YSF2:sin.SEAS")], 2, mean)
wg.priors[[6]] <- apply(eg[, c("YSF:cos.SEAS", "YSF2:cos.SEAS")], 2, mean)
wg.priors[[7]] <- mean(eg[sub3, c("eps.gam")])
wg.priors
