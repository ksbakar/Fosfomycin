

########################################################################################################################
## based on two-arm CT
########################################################################################################################

gen_data_arm <- function(n=c(300),rand_pr=NULL,n_arm=NULL){
  ##
  ## n is the sample in one interim 
  ##
  dat <- data.table::CJ(
    id = 1:n,
    id_arm = NA_integer_
  )
  ##
  if(is.null(rand_pr)){
    stop("provide value for rand_pr")
  }
  ##
  if(is.null(n_arm)){
    n_arm <- length(rand_pr)
    dat$id_arm <- sample(1:n_arm, size=nrow(dat), replace=TRUE, prob=rand_pr)
  }
  else{
    dat$id_arm <- sample(n_arm, size=nrow(dat), replace=TRUE, prob=rand_pr)
  }
  ##  
  dat
}
input <- gen_data_arm(n=100,rand_pr=c(0.5,0.5))
##
gen_data_fosfo <- function(input,theta){
  ## theta => soc and treat; theta=c(soc, treat)
  n <- table(input$id_arm)
  input[1:sum(n), evt := rbinom(n=sum(n), size=1, prob=theta[id_arm])]
  #for(i in 1:length(n)){
  #  input[input$id_arm==i,"evt"] <- rbinom(n=n[i], size=1, prob=theta[i])   
  #}
  return(data.frame(input))
}
input <- gen_data_fosfo(input=gen_data_arm(n=100,rand_pr=c(0.5,0.5)),theta=c(0.15,0.25))
##
trail_fosfo <- function(input,
                        delta=0.01,
                        a=1,b=1, # priors 
                        itr=1000){
  ##
  n <- table(input$id_arm)
  theta_draws <- matrix(NA,itr,length(n))
  for(i in 1:length(n)){ # number of arms
    theta_draws[,i] <- rbeta(itr, a + sum(input[input$id_arm==i,"evt"]), b + sum(!input[input$id_arm==i,"evt"]))
  }
  ##
  theta_post <- cbind(colMeans(theta_draws),apply(theta_draws,2,sd))
  dimnames(theta_post)[[2]] <- c("mean","sd")
  dimnames(theta_post)[[1]] <- paste0("arm",1:length(n))
  decision_id <- theta_draws[,1] + delta
  decision_pr <- NULL
  for(i in 2:length(n)){ # number of arms 
    decision_pr <- cbind(decision_pr,sum(theta_draws[,i]<decision_id)/nrow(theta_draws))
  }
  ##
  return(list(theta_post=theta_post,decision_pr=decision_pr))
  ##
}
trail_fosfo(input)
##
trial_interim_fosfo <- function(n_sample=c(100,50,50,50,50),
                                rand_assignment=c(1,1),
                                theta=c(0.15,0.25),delta=0.01,a=1,b=1,itr=1000,
                                D0=0.10,D1=0.97){
  ##
  n_j <- length(n_sample)
  K <- length(theta)
  n_interimSS <- matrix(NA,K,n_j)
  n_rand_pr <- matrix(NA,K,n_j)
  decision_pr <- matrix(NA,K-1,n_j)
  y_init <- list()
  out <- list()
  for(j in 1:n_j){
    if(j==1){
      ##
      N <- n_sample[j]
      rand_pr <- rand_assignment/sum(rand_assignment)
      y_init[[j]] = gen_data_fosfo(gen_data_arm(n=N,rand_pr=rand_pr),theta=theta)
      n_interimSS[,j] <- table(y_init[[j]]$id_arm)
      n_rand_pr[,j] <- rand_pr
      out[[j]] = trail_fosfo(y_init[[j]],delta=delta,a=a,b=b,itr=itr)
      ##
    }
    else if(j==n_j){
      ##
      N <- n_sample[j]
      ##
      ## following code only works for 2 arms 
      if(isTRUE(out[[j-1]]$decision_pr < D0) | isTRUE(out[[j-1]]$decision_pr > D1)){
      #if(isTRUE(out[[j-1]]$decision_pr < D0)){
        rand_pr[2] <- 0
        rand_pr <- rand_pr/sum(rand_pr)
      }
      else{
        #rand_pr <- rand_assignment/sum(rand_assignment)
        rand_pr <- rand_pr
      }
      ##
      ##
      y_init[[j]] = data.frame(mapply(c,y_init[[j-1]],gen_data_fosfo(gen_data_arm(n=N,rand_pr=rand_pr),theta=theta),SIMPLIFY=FALSE))
      n_interimSS[,j] <- table(y_init[[j]]$id_arm)
      n_rand_pr[,j] <- rand_pr
      out[[j]] = trail_fosfo(y_init[[j]],delta=delta,a=a,b=b,itr=itr)
      ##
    }
    else{
      ##
      N <- n_sample[j]
      ##
      ## following code only works for 2 arms 
      if(isTRUE(out[[j-1]]$decision_pr < D0) | isTRUE(out[[j-1]]$decision_pr > D1)){
      #if(isTRUE(out[[j-1]]$decision_pr < D0)){        
        rand_pr[2] <- 0
        rand_pr <- rand_pr/sum(rand_pr)
      }
      else{
        #rand_pr <- rand_assignment/sum(rand_assignment)
        rand_pr <- rand_pr
      }
      ##
      ##
      y_init[[j]] = data.frame(mapply(c,y_init[[j-1]],gen_data_fosfo(gen_data_arm(n=N,rand_pr=rand_pr),theta=theta),SIMPLIFY=FALSE))
      n_interimSS[,j] <- table(y_init[[j]]$id_arm)
      n_rand_pr[,j] <- rand_pr
      out[[j]] = trail_fosfo(y_init[[j]],delta=delta,a=a,b=b,itr=itr)
      ##
    }
  }
  ##
  decision_pr[,] <- sapply(out, function(x) x$decision_pr)
  theta_post_mn <- sapply(out, function(x) x$theta_post[,1])
  theta_post_sd <- sapply(out, function(x) x$theta_post[,2])
  ##
  dimnames(n_interimSS)[[1]] = paste0("arm",1:length(rand_assignment))
  dimnames(n_rand_pr)[[1]] = paste0("arm",1:length(rand_assignment))
  dimnames(decision_pr)[[1]] = paste0("arm",2:length(rand_assignment))
  ##
  return(list(theta_post_mn=theta_post_mn,theta_post_sd=theta_post_sd,decision_pr=decision_pr,
              n_interimSS=n_interimSS,n_rand_pr=n_rand_pr))
}
trial_interim_fosfo()
##
trial_interim_simulator_fosfo <- function(nSim,n_sample=c(100,50,50,50,50),
                                          rand_assignment=c(1,1),
                                          theta=c(0.15,0.25),delta=0.01,a=1,b=1,itr=1000,
                                          D0=0.1,D1=0.97){
  ##
  library(pbapply)
  out <- pblapply(1:nSim, 
                  function(x) 
                    trial_interim_fosfo(n_sample=n_sample,rand_assignment=rand_assignment,
                                        theta=theta,delta=delta,a=a,b=b,itr=itr,D0=D0,D1=D1))
  ##
  decision_pr <- t(sapply(out, function(x) x$decision_pr))
  dimnames(decision_pr)[[2]] = paste0("interim",1:length(n_sample))
  ##
  sample_size <- sapply(out, function(x) x$n_interimSS)
  sample_size <- array(sample_size,dim=c(length(rand_assignment),length(n_sample),nSim))
  dimnames(sample_size)[[1]] = c("soc",paste0("treat",1:(length(rand_assignment)-1)))
  dimnames(sample_size)[[2]] = paste0("interim",1:length(n_sample))
  ##
  rand_pr <- sapply(out, function(x) x$n_rand_pr)
  rand_pr <- array(rand_pr,dim=c(length(rand_assignment),length(n_sample),nSim))
  dimnames(rand_pr)[[1]] = c("soc",paste0("treat",1:(length(rand_assignment)-1)))
  dimnames(rand_pr)[[2]] = paste0("interim",1:length(n_sample))
  ##
  theta_mn <- sapply(out, function(x) x$theta_post_mn)
  theta_mn <- array(theta_mn,dim=c(length(rand_assignment),length(n_sample),nSim))
  dimnames(theta_mn)[[1]] = c("soc",paste0("treat",1:(length(rand_assignment)-1)))
  dimnames(theta_mn)[[2]] = paste0("interim",1:length(n_sample))
  ##
  theta_sd <- sapply(out, function(x) x$theta_post_sd)
  theta_sd <- array(theta_sd,dim=c(length(rand_assignment),length(n_sample),nSim))
  dimnames(theta_sd)[[1]] = c("soc",paste0("treat",1:(length(rand_assignment)-1)))
  dimnames(theta_sd)[[2]] = paste0("interim",1:length(n_sample))
  ##
  return(list(sample_size = sample_size,
              decision_pr = decision_pr,
              rand_pr = rand_pr,
              theta_mn = theta_mn,
              theta_sd = theta_sd,
              nSim=nSim, theta=theta, delta=delta,a=a,b=b,itr=itr))
  ##
}
##
res <- trial_interim_simulator_fosfo(nSim=50)
##

## simulation 

theta = cbind(seq(0.05,0.35,0.05),seq(0.05,0.35,0.05))
nSim = 1000
n_sample = c(100,50,50,50,50)
delta = 0.10
rand_assignment = c(1,1)
D0 = 0.1; D1 = 0.97
store_res <- NULL
id_scenario <- 0
for(i in 1:nrow(theta)){
  for(j in 1:nrow(theta)){
    res <- trial_interim_simulator_fosfo(nSim=nSim,n_sample=n_sample,rand_assignment=rand_assignment,
                               theta=c(theta[i,1],theta[j,2]),delta=delta,a=1,b=1,itr=1000,
                               D0=D0,D1=D1)
    ##
    ## summary stat
    ##
    dat <- data.table::CJ(
      id_scenario = NA_integer_,
      interim = 1:length(n_sample),
      arm = 1:length(rand_assignment),
      delta = delta,
      true_theta = NA_real_,
      est_theta = NA_real_,
      est_theta_sd = NA_real_,
      est_sample =  NA_real_,
      est_sd =  NA_real_
    )
    id_scenario <- id_scenario + 1
    dat$id_scenario <-  id_scenario
    dat$true_theta = rep(c(theta[i,1],theta[j,2]),length(n_sample))
    for(jj in 1:length(n_sample)){
      ##
      dat[dat$interim==jj,"est_theta"] <- rowMeans(res$theta_mn[,jj,])
      dat[dat$interim==jj,"est_theta_sd"] <- rowMeans(res$theta_sd[,jj,])
      ##
      dat[dat$interim==jj,"est_sample"] <- rowMeans(res$sample_size[,jj,])
      dat[dat$interim==jj,"est_sd"] <- apply(res$sample_size[,jj,],1,sd)
      ##
    }
    ##
    D0=0.1; D1=0.97
    pr_decision <- NULL
    for(j in 1:length(n_sample)){
      pr_decision <- rbind(pr_decision,cbind(j,t(rbind(mean(res$decision_pr[,j]<D0),mean(res$decision_pr[,j]>D1)))))
    }
    dimnames(pr_decision)[[2]] <- c("interim","inferior_D0_0.1","noninferior_D1_0.97")
    pr_decision <- data.frame(pr_decision)
    pr_decision$arm <- 2
    ##
    dat <- merge(dat,pr_decision,by=c("interim","arm"),all.x=TRUE)
    #dat
    ##
    store_res <- rbind(store_res,dat)
    write.csv(store_res,file=paste0("result_fosfo_D0.10_D1.97_",Sys.Date(),".csv"),row.names=FALSE)
  }
}

###



##
## summary stat
##
i = 1
dat <- data.table::CJ(
  id_scenario = NA_integer_,
  interim = 1:length(n_sample),
  arm = 1:length(rand_assignment),
  true_theta = NA_real_,
  est_theta = NA_real_,
  est_theta_sd = NA_real_,
  est_sample =  NA_real_,
  est_sd =  NA_real_
)
dat
dat$id_scenario <- i 
dat$true_theta = rep(theta,length(n_sample))
for(j in 1:length(n_sample)){
  ##
  dat[dat$interim==j,"est_theta"] <- rowMeans(res$theta_mn[,j,])
  dat[dat$interim==j,"est_theta_sd"] <- rowMeans(res$theta_sd[,j,])
  ##
  dat[dat$interim==j,"est_sample"] <- rowMeans(res$sample_size[,j,])
  dat[dat$interim==j,"est_sd"] <- apply(res$sample_size[,j,],1,sd)
  ##
}
##
D0=0.1; D1=0.97
pr_decision <- NULL
for(j in 1:length(n_sample)){
  pr_decision <- rbind(pr_decision,cbind(j,t(rbind(mean(res$decision_pr[,j]<D0),mean(res$decision_pr[,j]>D1)))))
}
dimnames(pr_decision)[[2]] <- c("interim","inferior_D0_0.1","noninferior_D1_0.97")
pr_decision <- data.frame(pr_decision)
pr_decision$arm <- 2
##
dat <- merge(dat,pr_decision,by=c("interim","arm"),all.x=TRUE)
dat
##


########################################################################################################################
## The End
########################################################################################################################
