library(MuMIn)
load( "C:\\Users\\rolek.brian\\Documents\\Projects\\FGSP Misclassification\\Data\\final-data.Rdata")

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
dat2 <- dat[1,]*0 + 1

combos <- dredge(global,eval=F,
                 subset = 
                   dc(YSF.std, YSF.std:sin.SEAS) &&
                   dc(YSF.std, YSF.std:cos.SEAS) &&
                   dc(sin.SEAS, YSF.std:sin.SEAS) &&
                   dc(cos.SEAS, YSF.std:cos.SEAS) 
)

n.mods <- length(combos)
terms <- colnames(model.matrix(formula(combos[[n.mods]]),dat2))
n.terms <- length(terms)

combo.mat <- matrix(0,nrow=n.mods,ncol=n.terms)
colnames(combo.mat) <- terms

for (i in 1:n.mods){
  combo.mat[i,match(colnames(model.matrix(formula(combos[[i]]),dat2)),
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
