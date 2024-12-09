library(cmdstanr)
library(parallelly)
library(dplyr)

check_cmdstan_toolchain()

ncores <- 4

# set up stan model##### 
stanfile <- "multivariate_normal_2mix_semisupervised_testing"
file <- file.path(getwd(),paste0("models/stan_files/",stanfile,".stan"))


mod <- cmdstan_model(file)

# run stan model #####
stan_df <- readRDS(paste0(getwd(),"/data/stan_df.RDS"))

for(i in 1:4){
  set.seed(1)
  whichquantile = i
  sampall <- sample(1:4,size = length(stan_df$DEN_M),replace = T)
  stan_df_train <- stan_df[sampall!=whichquantile,]
  stan_df_train <- stan_df_train %>% arrange(PCRconf)
  stan_df_test <- stan_df[sampall==whichquantile,]
  
  N=length(stan_df_train$DEN_M)
  Ntest=length(stan_df_test$DEN_M)
  Nc2=sum(stan_df_train$PCRconf) 
  Nc1= N - Nc2
  
  y1 = array(as.vector(c(log10(stan_df_train$DEN_M+1),
    log10(stan_df_train$DEN_G+1),
    log2(stan_df_train$avg_D_HAI))), 
  dim=c(N,3))
  y2 = array(as.vector(c(log10(stan_df_train$IgM_C_DEN+1),
    log10(stan_df_train$IgG_C_DEN+1),
    log2(stan_df_train$avg_C_best))), 
  dim=c(N,3))
  
  y1_test = array(as.vector(c(log10(stan_df_test$DEN_M+1),
    log10(stan_df_test$DEN_G+1),
    log2(stan_df_test$avg_D_HAI))), 
  dim=c(Ntest,3))
  y2_test = array(as.vector(c(log10(stan_df_test$IgM_C_DEN+1),
    log10(stan_df_test$IgG_C_DEN+1),
    log2(stan_df_test$avg_C_best))), 
  dim=c(Ntest,3))
  
  mixture_data=list(N=N,
                    Ntest=Ntest,
                    Nc1=Nc1,
                    Nc2=Nc2,
                    y1=y1,
                    y2=y2,
                    y1_test=y1_test,
                    y2_test=y2_test,
                    C=2,
                    K= 2,
                    nP= 3)
  
  # Set the number of chains
  chains <- 4
  set.seed(i)
  randUnifr <- runif(2*chains, 0.75, 1)
  randUniftheta <- runif(chains, 0.5, .7)
  
  randUnifsigma01 <- runif(chains, .1, 1)
  randUnifsigma02 <- runif(chains, randUnifsigma01, 1)
  randUnifsigma03 <- runif(chains, randUnifsigma02, 2)
  
  randUnifsigma11 <- runif(chains, .1, 1)
  randUnifsigma12 <- runif(chains, randUnifsigma11, 1)
  randUnifsigma13 <- runif(chains, randUnifsigma12, 2)
  
  randUnifbeta01 <- runif(chains, 1, 2)
  randUnifbeta02 <- runif(chains, randUnifbeta01, 2)
  
  init_values <- list(
    chain1 = list(
      r = matrix(c(randUnifr[1],1-randUnifr[1], 1-randUnifr[2],randUnifr[2]), nrow = 2, byrow = TRUE),
      theta = array(c(randUniftheta[1],1-randUniftheta[1]), dim=c(2)),
      sigma0 = array(c(randUnifsigma01[1],randUnifsigma02[1],randUnifsigma03[1]), dim=c(3)),
      sigma1 = array(c(randUnifsigma11[1],randUnifsigma12[1],randUnifsigma13[1]), dim=c(3)),
      beta0 = array(c(randUnifbeta01[1],randUnifbeta02[1], runif(1, 6.5, 8.5)), dim=c(3)),
      beta1 = array(c(runif(1, -1, 1), runif(1, -1, 1), runif(1, -1, 1)), dim=c(3)),
      L = array(c(1, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, 0,runif(1, -.5, .5),
                  1, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, 0,runif(1, -.5, .5)),dim = c(2,3,3))
    ),
    chain2 = list(
      r = matrix(c(randUnifr[3],1-randUnifr[3], 1-randUnifr[4],randUnifr[4]), nrow = 2, byrow = TRUE),
      theta = array(c(randUniftheta[2],1-randUniftheta[2]), dim=c(2)),
      sigma0 = array(c(randUnifsigma01[2],randUnifsigma02[2],randUnifsigma03[2]), dim=c(3)),
      sigma1 = array(c(randUnifsigma11[2],randUnifsigma12[2],randUnifsigma13[2]), dim=c(3)),
      beta0 = array(c(randUnifbeta01[2],randUnifbeta02[2], runif(1, 6.5, 8.5)), dim=c(3)),
      beta1 = array(c(runif(1, -1, 1), runif(1, -1, 1), runif(1, -1, 1)), dim=c(3)),
      L = array(c(1, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, 0,runif(1, -.5, .5),
                  1, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, 0,runif(1, -.5, .5)),dim = c(2,3,3))
    ),
    chain3 = list(
      r = matrix(c(randUnifr[5],1-randUnifr[5], 1-randUnifr[6],randUnifr[6]), nrow = 2, byrow = TRUE),
      theta = array(c(randUniftheta[3],1-randUniftheta[3]), dim=c(2)),
      sigma0 = array(c(randUnifsigma01[3],randUnifsigma02[3],randUnifsigma03[3]), dim=c(3)),
      sigma1 = array(c(randUnifsigma11[3],randUnifsigma12[3],randUnifsigma13[3]), dim=c(3)),
      beta0 = array(c(randUnifbeta01[3],randUnifbeta02[3], runif(1, 6.5, 8.5)), dim=c(3)),
      beta1 = array(c(runif(1, -1, 1), runif(1, -1, 1), runif(1, -1, 1)), dim=c(3)),
      L = array(c(1, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, 0,runif(1, -.5, .5),
                  1, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, 0,runif(1, -.5, .5)),dim = c(2,3,3))
    ),
    chain4 = list(
      r = matrix(c(randUnifr[7],1-randUnifr[7], 1-randUnifr[8],randUnifr[8]), nrow = 2, byrow = TRUE),
      theta = array(c(randUniftheta[4],1-randUniftheta[4]), dim=c(2)),
      sigma0 = array(c(randUnifsigma01[4],randUnifsigma02[4],randUnifsigma03[4]), dim=c(3)),
      sigma1 = array(c(randUnifsigma11[4],randUnifsigma12[4],randUnifsigma13[4]), dim=c(3)),
      beta0 = array(c(randUnifbeta01[4],randUnifbeta02[4], runif(1, 6.5, 8.5)), dim=c(3)),
      beta1 = array(c(runif(1, -1, 1), runif(1, -1, 1), runif(1, -1, 1)), dim=c(3)),
      L = array(c(1, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, 0,runif(1, -.5, .5),
                  1, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, runif(1, -.5, .5),runif(1, -.5, .5),
                  0, 0,runif(1, -.5, .5)),dim = c(2,3,3))
    )
    # Add similar sections for other chains if needed
  )
  
  fit <- mod$sample(
    data = mixture_data, 
    seed = 123, 
    chains = chains,
    parallel_chains = ncores,
    adapt_delta=0.99,
    iter_warmup = 2000,
    iter_sampling = 1000,
    refresh = 500 # print update every 500 iters
  )
  
  fit$save_object(file = paste0(getwd(),"/stan_output/fit_",stanfile,"_ACall_",whichquantile,".RDS")) 
  
}
