
##
## Artemis bathch file run: Date 22/02/2022
##

library(pbapply)

source("source_code.R")

## 
## historical theta (i.e., soc) = 0.15, 0.25, 0.35
## delta = seq(0.10,0.15,by=0.01)
## true theta = seq(0.15,0.35,by=0.05)
## nMax = c(400)
## allocation of arms (1:1)
##

nSim <- 5000
soc_theta <- seq(0.15,0.35,by=0.05)
true_theta <- seq(0.15,0.35,by=0.05)
delta <- seq(0.10,0.15,by=0.01)
soc_arm <- c(1)
nMax <- c(400)
for(nn in 1:length(nMax)){
  for(j in 1:length(soc_theta)){
    for(i in 1:length(true_theta)){
      for(l in 1:length(delta)){
        for(arm in 1:length(soc_arm)){
          res <- trial_interim_simulator_2(nSim=nSim,n_sample=seq(100,nMax[nn],50),itr=1000,nburn=0,
                                           stanrun=FALSE, arm2_output=TRUE,
                                           delta=delta[l],true_theta=true_theta[i],soc_theta=soc_theta[j],
                                           soc_arm=soc_arm[arm],true_arm=1)
          save(res,file=paste0("res_arm",soc_arm[arm],"1_delta",delta[l],"_theta",true_theta[i],"_soc",soc_theta[j],"_nMax",nMax[nn],".RData"))
          rm(res);
        }
      }
    }
  }  
}
##
##