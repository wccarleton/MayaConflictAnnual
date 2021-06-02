poisTSCode <- nimbleCode({
   if(K > 1){
      A[1:K,1:K] <- inverse(diag(x=rep(10,K)))
      multimu[1:K] <- rep(0,K)
      B[1:K] ~ dmnorm(mean=multimu[1:K], prec=A[1:K,1:K])
      mu[1] <- inprod(X[1,1:K], B[1:K])
   }else{
      B ~ dnorm(mean=0,sd=10)
      mu[1] <- X[1]*B
   }
   sigma ~ dexp(0.7)
   lambda0 ~ dnorm(mean=0,sd=10)
   rho ~ dnorm(mean=0,sd=10)
   lambda[1] ~ dnorm(rho * lambda0,sd=sigma)
   for (j in 2:T){
      if(K > 1){
         mu[j] <- inprod(X[j,1:K], B[1:K])
      }else{
         mu[j] <- X[j] * B
      }
      lambda[j] ~ dnorm(rho * lambda[j-1],sd=sigma)
   }
   for (j in 1:T){
      Y[j] ~ dpois(XBase[j]*exp(mu[j]+lambda[j]))
      y[j] ~ dpois(XBase[j]*exp(mu[j]+lambda[j]))
   }
})
conflict_df <- read.table("Data/Maya_clim_conflict.csv",head=T)
Y <- conflict_df$Conflict
T <- length(Y)
XBase <- conflict_df$Monuments_adj
SSTr = conflict_df$Cariaco_Rub_SST-mean(conflict_df$Cariaco_Rub_SST)
d18O = conflict_df$YokCave_d18O-mean(conflict_df$YokCave_d18O)
X <- cbind(rep(1,T),SSTr,d18O)
K <- ncol(X)

poisTSData <- list(Y=Y,
                  X=X,
                  XBase=XBase)

poisTSConsts <- list(T=T,
                     K=K)

poisTSInits <- list(lambda0=0,
                     sigma=1,
                     rho=0,
                     B=rep(0,K))

poisTSModel <- nimbleModel(code=poisTSCode,
                        data=poisTSData,
                        inits=poisTSInits,
                        constants=poisTSConsts)

#compile nimble model to C++ code—much faster runtime
C_poisTSModel <- compileNimble(poisTSModel, showCompilerOutput = FALSE)

#configure the MCMC
poisTSModel_conf <- configureMCMC(poisTSModel)

#select the variables that we want to monitor in the MCMC chain
poisTSModel_conf$addMonitors2("y")

poisTSModel_conf$removeSamplers(c("B","lambda","lambda0"))
poisTSModel_conf$addSampler(target="B",type="AF_slice")
poisTSModel_conf$addSampler(target=c("lambda0","lambda"),type="AF_slice")

#build MCMC
poisTSModelMCMC <- buildMCMC(poisTSModel_conf,thin=1,enableWAIC=F)

#compile MCMC to C++—much faster
C_poisTSModelMCMC <- compileNimble(poisTSModelMCMC,project=poisTSModel)

#number of MCMC iterations
niter = 100000
niter90 = niter * 0.9

#set seed for replicability
set.seed(1)

#save the MCMC chain (monitored variables) as a matrix
samples <- runMCMC(C_poisTSModelMCMC, niter=niter)

save(samples,file="./mcmc_samples.RData")

#diags
ncol_samples <- ncol(samples)
coda_mcmc <- mcmc(samples[1000:niter,],thin=1)
mcmc_geweke <- geweke.diag(coda_mcmc)
write.csv(t(mcmc_geweke$z),file="./mcmc_geweke.csv")

alarm()
