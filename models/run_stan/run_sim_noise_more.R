library(cmdstanr)
library(parallelly)
library(dplyr)

check_cmdstan_toolchain()

ncores <- ifelse(parallelly::availableCores()>4,8,parallelly::availableCores())

# set up stan model##### 
stanfile <- "multivariate_normal_2mix_semisupervised_testing"
file <- file.path(getwd(),paste0("models/stan_files/",stanfile,".stan"))
mod <- cmdstan_model(file)

# run stan model #####
sim_noise <- readRDS(paste0(getwd(),"/data/sim_noise_more.RDS"))

#IgG
for(i in 1:4){
  set.seed(1)
  whichsamp = i
  sampall <- sample(1:4,size = length(sim_noise$infected),replace = T)
  sim_noise_train <- sim_noise[sampall!=whichsamp,]
  sim_noise_train <- sim_noise_train %>% arrange(infected)
  sim_noise_test <- sim_noise[sampall==whichsamp,]

  N=length(sim_noise_train$IgM_t2)
  Ntest=length(log10(sim_noise_test$IgM_t2+1))
  Nc2=sum(sim_noise_train$infected)
  Nc1= N - Nc2
  nP = 1

  y1 = array(as.vector(c(log10(sim_noise_train$IgG_t2+1))), dim=c(N,nP))
  y2 = array(as.vector(c(log10(sim_noise_train$IgG_t3+1))), dim=c(N,nP))

  y1_test = array(as.vector(c((log10(sim_noise_test$IgG_t2+1)))), dim=c(Ntest,nP))
  y2_test = array(as.vector(c((log10(sim_noise_test$IgG_t3+1)))), dim=c(Ntest,nP))

  mixture_data=list(N=N,
                    Ntest=Ntest,
                    Nc1=Nc1,
                    Nc2=Nc2,
                    y1=y1,
                    y2=y2,
                    y1_test=y1_test,
                    y2_test=y2_test,
                    C = 2,
                    K = 2,
                    nP = nP)

  fit <- mod$sample(
    data = mixture_data,
    seed = 1,
    chains = 4,
    parallel_chains = ncores,
    adapt_delta=0.9,
    iter_warmup = 2000,
    iter_sampling = 1000,
    refresh = 500 # print update every 500 iters
  )

  fit$save_object(file = paste0(getwd(),"/stan_output/fit_",stanfile,"_sim_noise_more_nonLog_IgG_",whichsamp,".RDS"))

}

#HAI IgM
for(i in 1:4){
  set.seed(1)
  whichsamp = i
  sampall <- sample(1:4,size = length(sim_noise$infected),replace = T)
  sim_noise_train <- sim_noise[sampall!=whichsamp,]
  sim_noise_train <- sim_noise_train %>% arrange(infected)
  sim_noise_test <- sim_noise[sampall==whichsamp,]

  N=length(sim_noise_train$IgM_t2)
  Ntest=length(log10(sim_noise_test$IgM_t2+1))
  Nc2=sum(sim_noise_train$infected)
  Nc1= N - Nc2
  nP = 2

  y1 = array(as.vector(c(log10(sim_noise_train$IgM_t2+1),(sim_noise_train$HAI_t2))), dim=c(N,nP))
  y2 = array(as.vector(c(log10(sim_noise_train$IgM_t3+1),(sim_noise_train$HAI_t3))), dim=c(N,nP))

  y1_test = array(as.vector(c(log10(sim_noise_test$IgM_t2+1),(sim_noise_test$HAI_t2))), dim=c(Ntest,nP))
  y2_test = array(as.vector(c(log10(sim_noise_test$IgM_t3+1),(sim_noise_test$HAI_t3))), dim=c(Ntest,nP))

  mixture_data=list(N=N,
                    Ntest=Ntest,
                    Nc1=Nc1,
                    Nc2=Nc2,
                    y1=y1,
                    y2=y2,
                    y1_test=y1_test,
                    y2_test=y2_test,
                    C = 2,
                    K = 2,
                    nP = nP)

  fit <- mod$sample(
    data = mixture_data,
    seed = 1,
    chains = 4,
    parallel_chains = ncores,
    adapt_delta=0.9,
    iter_warmup = 2000,
    iter_sampling = 1000,
    refresh = 500 # print update every 500 iters
  )

  fit$save_object(file = paste0(getwd(),"/stan_output/fit_",stanfile,"_sim_noise_more_nonLog_HAI_IgM_",whichsamp,".RDS"))

}

#IgM IgG
for(i in 1:4){
  set.seed(1)
  whichsamp = i
  sampall <- sample(1:4,size = length(sim_noise$infected),replace = T)
  sim_noise_train <- sim_noise[sampall!=whichsamp,]
  sim_noise_train <- sim_noise_train %>% arrange(infected)
  sim_noise_test <- sim_noise[sampall==whichsamp,]

  N=length(sim_noise_train$IgM_t2)
  Ntest=length(log10(sim_noise_test$IgM_t2+1))
  Nc2=sum(sim_noise_train$infected)
  Nc1= N - Nc2
  nP = 2

  y1 = array(as.vector(c(log10(sim_noise_train$IgM_t2+1),log10(sim_noise_train$IgG_t2+1))), dim=c(N,nP))
  y2 = array(as.vector(c(log10(sim_noise_train$IgM_t3+1),log10(sim_noise_train$IgG_t3+1))), dim=c(N,nP))

  y1_test = array(as.vector(c(log10(sim_noise_test$IgM_t2+1),(log10(sim_noise_test$IgG_t2+1)))), dim=c(Ntest,nP))
  y2_test = array(as.vector(c(log10(sim_noise_test$IgM_t3+1),(log10(sim_noise_test$IgG_t3+1)))), dim=c(Ntest,nP))

  mixture_data=list(N=N,
                    Ntest=Ntest,
                    Nc1=Nc1,
                    Nc2=Nc2,
                    y1=y1,
                    y2=y2,
                    y1_test=y1_test,
                    y2_test=y2_test,
                    C = 2,
                    K = 2,
                    nP = nP)

  fit <- mod$sample(
    data = mixture_data,
    seed = 1,
    chains = 4,
    parallel_chains = ncores,
    adapt_delta=0.9,
    iter_warmup = 2000,
    iter_sampling = 1000,
    refresh = 500 # print update every 500 iters
  )

  fit$save_object(file = paste0(getwd(),"/stan_output/fit_",stanfile,"_sim_noise_more_nonLog_IgM_IgG_",whichsamp,".RDS"))

}

#HAI IgG
for(i in 1:4){
  set.seed(1)
  whichsamp = i
  sampall <- sample(1:4,size = length(sim_noise$infected),replace = T)
  sim_noise_train <- sim_noise[sampall!=whichsamp,]
  sim_noise_train <- sim_noise_train %>% arrange(infected)
  sim_noise_test <- sim_noise[sampall==whichsamp,]

  N=length(sim_noise_train$IgM_t2)
  Ntest=length(log10(sim_noise_test$IgM_t2+1))
  Nc2=sum(sim_noise_train$infected)
  Nc1= N - Nc2
  nP = 2

  y1 = array(as.vector(c(log10(sim_noise_train$IgG_t2+1),(sim_noise_train$HAI_t2))), dim=c(N,nP))
  y2 = array(as.vector(c(log10(sim_noise_train$IgG_t3+1),(sim_noise_train$HAI_t3))), dim=c(N,nP))

  y1_test = array(as.vector(c(log10(sim_noise_test$IgG_t2+1),(sim_noise_test$HAI_t2))), dim=c(Ntest,nP))
  y2_test = array(as.vector(c(log10(sim_noise_test$IgG_t3+1),(sim_noise_test$HAI_t3))), dim=c(Ntest,nP))

  mixture_data=list(N=N,
                    Ntest=Ntest,
                    Nc1=Nc1,
                    Nc2=Nc2,
                    y1=y1,
                    y2=y2,
                    y1_test=y1_test,
                    y2_test=y2_test,
                    C = 2,
                    K = 2,
                    nP = nP)

  fit <- mod$sample(
    data = mixture_data,
    seed = 1,
    chains = 4,
    parallel_chains = ncores,
    adapt_delta=0.9,
    iter_warmup = 2000,
    iter_sampling = 1000,
    refresh = 500 # print update every 500 iters
  )

  fit$save_object(file = paste0(getwd(),"/stan_output/fit_",stanfile,"_sim_noise_more_nonLog_HAI_IgG_",whichsamp,".RDS"))

}

# all
for(i in 1:4){
  set.seed(1)
  whichsamp = i
  sampall <- sample(1:4,size = length(sim_noise$infected),replace = T)
  sim_noise_train <- sim_noise[sampall!=whichsamp,]
  sim_noise_train <- sim_noise_train %>% arrange(infected)
  sim_noise_test <- sim_noise[sampall==whichsamp,]

  N=length(sim_noise_train$IgM_t2)
  Ntest=length(log10(sim_noise_test$IgM_t2+1))
  Nc2=sum(sim_noise_train$infected)
  Nc1= N - Nc2
  nP = 3

  y1 = array(as.vector(c(log10(sim_noise_train$IgM_t2+1),log10(sim_noise_train$IgG_t2+1),(sim_noise_train$HAI_t2))), dim=c(N,nP))
  y2 = array(as.vector(c(log10(sim_noise_train$IgM_t3+1),log10(sim_noise_train$IgG_t3+1),(sim_noise_train$HAI_t3))), dim=c(N,nP))

  y1_test = array(as.vector(c(log10(sim_noise_test$IgM_t2+1),(log10(sim_noise_test$IgG_t2+1)),(sim_noise_test$HAI_t2))), dim=c(Ntest,nP))
  y2_test = array(as.vector(c(log10(sim_noise_test$IgM_t3+1),(log10(sim_noise_test$IgG_t3+1)),(sim_noise_test$HAI_t3))), dim=c(Ntest,nP))

  mixture_data=list(N=N,
                    Ntest=Ntest,
                    Nc1=Nc1,
                    Nc2=Nc2,
                    y1=y1,
                    y2=y2,
                    y1_test=y1_test,
                    y2_test=y2_test,
                    C = 2,
                    K = 2,
                    nP = nP)

  fit <- mod$sample(
    data = mixture_data,
    seed = 1,
    chains = 4,
    parallel_chains = ncores,
    adapt_delta=0.9,
    iter_warmup = 2000,
    iter_sampling = 1000,
    refresh = 500 # print update every 500 iters
  )

  fit$save_object(file = paste0(getwd(),"/stan_output/fit_",stanfile,"_sim_noise_more_nonLog_",whichsamp,".RDS"))

}

# 
# # set up stan model#####
# stanfile <- "multivariate_normal_2mix_semisupervised_testing_nor22prior"
# file <- file.path(getwd(),paste0("models/stan_files/",stanfile,".stan"))
# 
# 
# mod <- cmdstan_model(file)

# HAI
for(i in 1:4){
  set.seed(1)
  whichsamp = i
  sampall <- sample(1:4,size = length(sim_noise$infected),replace = T)
  sim_noise_train <- sim_noise[sampall!=whichsamp,]
  sim_noise_train <- sim_noise_train %>% arrange(infected)
  sim_noise_test <- sim_noise[sampall==whichsamp,]
  
  N=length(sim_noise_train$IgM_t2)
  Ntest=length(log10(sim_noise_test$IgM_t2+1))
  Nc2=sum(sim_noise_train$infected)
  Nc1= N - Nc2
  nP = 1
  
  y1 = array(as.vector(c((sim_noise_train$HAI_t2))), dim=c(N,nP))
  y2 = array(as.vector(c((sim_noise_train$HAI_t3))), dim=c(N,nP))
  
  y1_test = array(as.vector(c((sim_noise_test$HAI_t2))), dim=c(Ntest,nP))
  y2_test = array(as.vector(c((sim_noise_test$HAI_t3))), dim=c(Ntest,nP))
  
  mixture_data=list(N=N,
                    Ntest=Ntest,
                    Nc1=Nc1,
                    Nc2=Nc2,
                    y1=y1,
                    y2=y2,
                    y1_test=y1_test,
                    y2_test=y2_test,
                    C = 2,
                    K = 2,
                    nP = nP)
  
  fit <- mod$sample(
    data = mixture_data,
    seed = 1,
    chains = 4,
    parallel_chains = ncores,
    adapt_delta=0.9,
    iter_warmup = 2000,
    iter_sampling = 1000,
    refresh = 500 # print update every 500 iters
  )
  
  fit$save_object(file = paste0(getwd(),"/stan_output/fit_",stanfile,"_sim_noise_more_nonLog_HAI_",whichsamp,".RDS"))
  
}

#IgM
for(i in 1:4){
  set.seed(1)
  whichsamp = i
  sampall <- sample(1:4,size = length(sim_noise$infected),replace = T)
  sim_noise_train <- sim_noise[sampall!=whichsamp,]
  sim_noise_train <- sim_noise_train %>% arrange(infected)
  sim_noise_test <- sim_noise[sampall==whichsamp,]

  N=length(sim_noise_train$IgM_t2)
  Ntest=length(log10(sim_noise_test$IgM_t2+1))
  Nc2=sum(sim_noise_train$infected)
  Nc1= N - Nc2
  nP = 1

  y1 = array(as.vector(c(log10(sim_noise_train$IgM_t2+1))), dim=c(N,nP))
  y2 = array(as.vector(c(log10(sim_noise_train$IgM_t3+1))), dim=c(N,nP))

  y1_test = array(as.vector(c((log10(sim_noise_test$IgM_t2+1)))), dim=c(Ntest,nP))
  y2_test = array(as.vector(c((log10(sim_noise_test$IgM_t3+1)))), dim=c(Ntest,nP))

  mixture_data=list(N=N,
                    Ntest=Ntest,
                    Nc1=Nc1,
                    Nc2=Nc2,
                    y1=y1,
                    y2=y2,
                    y1_test=y1_test,
                    y2_test=y2_test,
                    C = 2,
                    K = 2,
                    nP = nP)

  fit <- mod$sample(
    data = mixture_data,
    seed = 1,
    chains = 4,
    parallel_chains = ncores,
    adapt_delta=0.9,
    iter_warmup = 2000,
    iter_sampling = 1000,
    refresh = 500 # print update every 500 iters
  )

  fit$save_object(file = paste0(getwd(),"/stan_output/fit_",stanfile,"_sim_noise_more_nonLog_IgM_",whichsamp,".RDS"))

}
