

########################################################################################################################
## based on two-arm CT
########################################################################################################################
gen_data_2 <- function(N_true,N_soc,true_theta,soc_theta,arm2_output=TRUE){
  y1 <- rbinom(n=N_true, size=1, prob=true_theta)
  y2 <- rbinom(n=N_soc, size=1, prob=soc_theta)
  if(isTRUE(arm2_output)){
    return(list(y_soc=y2,y_true=y1))
  }
  else{
    y <- c(y1,y2)
    y
  }
}
##
trial_2 <- function(y,
                    true_theta=0.05,
                    soc_theta=0.05,
                    delta=0.01,
                    a=1,b=1,
                    chains=2,
                    itr=1000,
                    nburn=500,
                    arm2_output=TRUE,
                    stanrun=FALSE,
                    BF=FALSE){
  ##
  ## true_theta = to simulate data
  ## soc_theta = historical value for the SOC
  ## delta = tolerance 
  ## N = sample size
  ## a,b = hyper-parameter for the beta distribution
  ##
  if(isTRUE(stanrun)){
    stop("not possible for two-arm trial")
  }
  else{
    if(isTRUE(arm2_output)){
      y_soc <- y$y_soc
      y_true <- y$y_true
      theta_draws_soc <- rbeta(itr, a + sum(y_soc), b + sum(!y_soc))
      theta_draws_true <- rbeta(itr, a + sum(y_true), b + sum(!y_true))
      theta_post_mean_soc <- mean(theta_draws_soc)
      theta_post_mean_true <- mean(theta_draws_true)
      decision_id <- theta_draws_soc + delta
      decision_pr <- sum(theta_draws_true<decision_id)/length(theta_draws_true)
      N <- length(c(y_soc,y_true))
      if(isTRUE(BF)){
        x <- theta_draws_true*N
        fnc_bf10 <- function(x,n,th0,a0,b0,a1,b1){
          1/hypotestBAF5x(x=x,n=n,th0=th0,a0=a0,b0=b0,a1=a1,b1=b1)$BaFa01
        }
        decision_id <- theta_post_mean_soc + delta
        bf10 <- sapply(x,fnc_bf10,n=N,th0=decision_id,a0=1,b0=1,a1=1,b1=1)
        decision_bf10_v3 <- sum(bf10>3)/length(bf10)
        decision_bf10_v10 <- sum(bf10>10)/length(bf10)
        decision_bf10_v30 <- sum(bf10>30)/length(bf10)
        decision_bf10_v100 <- sum(bf10>100)/length(bf10)
      }
      theta_post_mean <- theta_post_mean_true
    }
    else{
      N <- length(y)
      theta_draws <- rbeta(itr, a + sum(y), b + sum(!y))
      ##
      decision_id <- soc_theta + delta
      decision_pr <- sum(theta_draws<decision_id)/length(theta_draws)
      theta_post_mean <- mean(theta_draws)
      if(isTRUE(BF)){
        x <- theta_draws*N
        fnc_bf10 <- function(x,n,th0,a0,b0,a1,b1){
          1/hypotestBAF5x(x=x,n=n,th0=th0,a0=a0,b0=b0,a1=a1,b1=b1)$BaFa01
        }
        bf10 <- sapply(x,fnc_bf10,n=N,th0=decision_id,a0=1,b0=1,a1=1,b1=1)
        decision_bf10_v3 <- sum(bf10>3)/length(bf10)
        decision_bf10_v10 <- sum(bf10>10)/length(bf10)
        decision_bf10_v30 <- sum(bf10>30)/length(bf10)
        decision_bf10_v100 <- sum(bf10>100)/length(bf10)
      }
    }
  }
  if(isTRUE(BF)){
    return(list(sample.size=N,decision.pr=decision_pr,
                theta_post_mean=theta_post_mean,
                decision_bf10_v3=decision_bf10_v3,
                decision_bf10_v10=decision_bf10_v10,
                decision_bf10_v30=decision_bf10_v30,
                decision_bf10_v100=decision_bf10_v100))
  }
  else{
    return(list(sample.size=N,decision.pr=decision_pr,
                theta_post_mean=theta_post_mean))
  }
}
##
trial_interim_2 <- function(n_sample,
                            true_theta=0.05,
                            soc_theta=0.05,
                            soc_arm=1,
                            true_arm=1,
                            delta=0.01,
                            a=1,b=1,
                            chains=2,
                            itr=1000,
                            nburn=500,
                            arm2_output=TRUE,
                            stanrun=FALSE,
                            BF=FALSE){
  ##
  n_j <- length(n_sample)
  y_init <- list()
  for(j in 1:n_j){
    if(j%in%1){
      N <- n_sample[j]
      N_true <- true_arm*ceiling(N/(soc_arm+true_arm))
      N_soc <- soc_arm*ceiling(N/(soc_arm+true_arm))
      y_init[[j]] <- gen_data_2(N_true = N_true,
                                true_theta = true_theta,
                                N_soc = N_soc,
                                soc_theta = soc_theta)
    }
    else{
      N <- n_sample[j]-n_sample[j-1]
      N_true <- true_arm*ceiling(N/(soc_arm+true_arm))
      N_soc <- soc_arm*ceiling(N/(soc_arm+true_arm))
      y_init[[j]] <- mapply(c,y_init[[j-1]],gen_data_2(N_true = N_true, 
                                                       true_theta = true_theta,
                                                       N_soc = N_soc,
                                                       soc_theta = soc_theta),SIMPLIFY=FALSE)
    }
  }
  ##
  out <- lapply(y_init, function(x) trial_2(y = x,
                                            true_theta=true_theta,
                                            soc_theta=soc_theta,
                                            delta=delta,
                                            a=a,b=b,
                                            chains=chains,
                                            itr=itr,
                                            nburn=nburn,
                                            arm2_output=arm2_output,
                                            stanrun=stanrun,
                                            BF=BF))
  ##
  decision_pr <- sapply(out, function(x) x$decision.pr)
  theta_post_mean <- sapply(out, function(x) x$theta_post_mean)
  ##
  if(isTRUE(BF)){
    decision_bf10_v3 <- sapply(out, function(x) x$decision_bf10_v3)
    decision_bf10_v10 <- sapply(out, function(x) x$decision_bf10_v10)
    decision_bf10_v30 <- sapply(out, function(x) x$decision_bf10_v30)
    decision_bf10_v100 <- sapply(out, function(x) x$decision_bf10_v100)
    return(list(decision_pr = decision_pr, 
                theta_post_mean=theta_post_mean, 
                decision_bf10_v3=decision_bf10_v3,
                decision_bf10_v10=decision_bf10_v10,
                decision_bf10_v30=decision_bf10_v30,
                decision_bf10_v100=decision_bf10_v100))
  }
  else{
    return(list(decision_pr = decision_pr, 
                theta_post_mean=theta_post_mean))
  }
}
##
trial_interim_simulator_2 <- function(n_sample=seq(100,500,50),
                                      nSim=10, 
                                      true_theta=0.05,
                                      soc_theta=0.05,
                                      soc_arm=1,
                                      true_arm=1,
                                      delta=0.01,
                                      a=1,b=1,
                                      chains=1,
                                      itr=500,
                                      nburn=100,
                                      arm2_output=TRUE,
                                      stanrun=FALSE,
                                      BF=FALSE){
  ##
  out <- pblapply(1:nSim, 
                  function(x) 
                    trial_interim_2(n_sample=n_sample, 
                                    true_theta=true_theta,
                                    soc_theta=soc_theta,
                                    soc_arm=soc_arm,
                                    true_arm=true_arm,
                                    delta=delta,
                                    a=a,b=b,
                                    chains=chains,
                                    itr=itr,
                                    nburn=nburn,
                                    stanrun=stanrun,
                                    BF=BF))
  ##
  sample_size <- n_sample
  decision_pr <- sapply(out, function(x) x$decision_pr)
  dimnames(decision_pr)[[1]] <- paste0("ss",sample_size)
  theta_post_mean <- sapply(out, function(x) x$theta_post_mean)
  dimnames(theta_post_mean)[[1]] <- paste0("ss",sample_size)
  ##
  if(isTRUE(BF)){
    decision_bf10_v3 <- sapply(out, function(x) x$decision_bf10_v3)
    dimnames(decision_bf10_v3)[[1]] <- paste0("ss",sample_size)
    decision_bf10_v10 <- sapply(out, function(x) x$decision_bf10_v10)
    dimnames(decision_bf10_v10)[[1]] <- paste0("ss",sample_size)
    decision_bf10_v30 <- sapply(out, function(x) x$decision_bf10_v30)
    dimnames(decision_bf10_v30)[[1]] <- paste0("ss",sample_size)
    decision_bf10_v100 <- sapply(out, function(x) x$decision_bf10_v100)
    dimnames(decision_bf10_v100)[[1]] <- paste0("ss",sample_size)
    ##
    return(list(sample_size=sample_size,
                decision_pr=t(decision_pr),
                theta_post_mean=t(theta_post_mean), 
                decision_bf10_v3=t(decision_bf10_v3),
                decision_bf10_v10=t(decision_bf10_v10),
                decision_bf10_v30=t(decision_bf10_v30),
                decision_bf10_v100=t(decision_bf10_v100),
                n_sample=n_sample, nSim=nSim,
                true_theta=true_theta, soc_theta=soc_theta,
                delta=delta,a=a,b=b,chains=chains,itr=itr,nburn=nburn))
    
  }
  else{
    ##
    return(list(sample_size=sample_size,
                decision_pr=t(decision_pr),
                theta_post_mean=t(theta_post_mean), 
                n_sample=n_sample, nSim=nSim,
                true_theta=true_theta, soc_theta=soc_theta,
                delta=delta,a=a,b=b,chains=chains,itr=itr,nburn=nburn))
  }
  ##
}

## decision function based on probability

decision_fnc <- function(pr,D1,D2,nMax){
  ##
  pr <- pr[,1:which(dimnames(pr)[[2]]%in%paste0("ss",nMax))]
  ##
  unsafe <- c()
  safe <- c()
  for(i in 1:nrow(pr)){
    safe[i] <- names(which(pr[i,]>D2)[1]) # safe
    unsafe[i] <- names(which(pr[i,]<D1)[1]) # unsafe
  }
  Pr_safe <- length(na.omit(safe))/length(safe)
  Pr_unsafe <- length(na.omit(unsafe))/length(unsafe)
  sample_size <- c(unsafe[complete.cases(unsafe)])
  sample_size <- as.numeric(gsub("ss","",sample_size))
  sample_size <- c(sample_size, rep(nMax, nrow(pr) - length(sample_size)))
  Pr_stop_early <- mean(sample_size < nMax)
  dat <- data.frame(nMax=nMax, D1=D1, D2=D2,
                    Pr_safe=Pr_safe,Pr_unsafe=Pr_unsafe,Pr_stop_early=Pr_stop_early,
                    Mean_sample_size=mean(sample_size),SD_sample_size=sd(sample_size))
  round(dat,2)
}

########################################################################################################################
########################################################################################################################