## Install and load serosim along with other packages #####
library(serosim)
library(tidyverse)
library(patchwork)
library(reshape2)
library(tibble)

# source an updated function for simulation
source("~/Documents/GitHub/DENV_modeling/Illness_inv_prob/Manuscript/data/runserosim_updated.R", echo=TRUE)

##### simulate #####

# create demography
times <- seq(1,365,by=1) # simulate monthly data for 1 year
N=500 # nubmer of indviduals simulated
demography <- generate_pop_demography(N=N, times=times, prob_removal=0,birth_times = 2)

## Create biomarker map
biomarker_map_original <- tibble(exposure_id=c("ifxn","ifxn","ifxn"),biomarker_id=c("HAI","IgG","IgM"))

## Reformat biomarker_map for runserosim
biomarker_map <-reformat_biomarker_map(biomarker_map_original)

## Create an empty array to store the force of exposure for all exposure types
foe_pars <- array(0, dim=c(1,max(times),n_distinct(biomarker_map$exposure_id)))
## Specify the force of exposure for exposure ID 1 which represents natural infection
foe_pars[,1:200,1] <- (.007+.001*cos(2*pi*c(1:200)/400))
foe_pars[,331:335,1] <- .09 

## Bring in the antibody parameters needed for the antibody model
## Note that the observation error parameter needed for the observation model (Section 1.7) is defined here too.
model_pars_path <- system.file("extdata", "model_pars_cs2.csv", package = "serosim")
model_pars_original <- read.csv(file = model_pars_path, header = TRUE)
model_pars_original <- model_pars_original[1,]
model_pars_original[1,] <- c("ifxn","HAI","boost_long",4,1,"log-normal")
model_pars_original[2,] <- c("ifxn","HAI","boost_short",0.75,.1,"log-normal")
model_pars_original[3,] <- c("ifxn","HAI","wane_long",.0001,.0001,"log-normal")
model_pars_original[4,] <- c("ifxn","HAI","wane_short",.005,.001,"log-normal")
model_pars_original[5,] <- c("ifxn","HAI","biomarker_prot_midpoint",8,NA,"")
model_pars_original[6,] <- c("ifxn","HAI","biomarker_prot_width",.8,NA,"")
model_pars_original[7,] <- c(NA,"HAI","obs_sd",NA,.05,"normal")
model_pars_original[8,] <- c("ifxn","IgG","boost_long",5,30,"log-normal")
model_pars_original[9,] <- c("ifxn","IgG","boost_short",80,30,"log-normal")
model_pars_original[10,] <- c("ifxn","IgG","wane_long",.005,.005,"log-normal")
model_pars_original[11,] <- c("ifxn","IgG","wane_short",.02,.05,"log-normal")
model_pars_original[12,] <- c(NA,"IgG","obs_sd",NA,1.5,"normal")
model_pars_original[13,] <- c("ifxn","IgM","boost_long",5,30,"log-normal")
model_pars_original[14,] <- c("ifxn","IgM","boost_short",50,70,"log-normal")
model_pars_original[15,] <- c("ifxn","IgM","wane_long",.005,.005,"log-normal")
model_pars_original[16,] <- c("ifxn","IgM","wane_short",.0075,.05,"log-normal")
model_pars_original[17,] <- c(NA,"IgM","obs_sd",NA,1.5,"normal")
model_pars_original[18,] <- c("ifxn","HAI","biomarker_ceiling_threshold",10,NA,"")
model_pars_original[19,] <- c("ifxn","HAI","biomarker_ceiling_gradient",.15,NA,"")
model_pars_original[20,] <- c("ifxn","IgG","biomarker_ceiling_threshold",80,NA,"")
model_pars_original[21,] <- c("ifxn","IgG","biomarker_ceiling_gradient",.01,NA,"")
model_pars_original[22,] <- c("ifxn","IgM","biomarker_ceiling_threshold",75,NA,"")
model_pars_original[23,] <- c("ifxn","IgM","biomarker_ceiling_gradient",.01,NA,"")
model_pars<-reformat_biomarker_map(model_pars_original)

model_pars$mean <- as.numeric(model_pars$mean)
model_pars$sd <- as.numeric(model_pars$sd)

bounds<-tibble(biomarker_id=c(1,1,2,2,3,3),name=rep(c("lower_bound","upper_bound"),3),
               value=c(log2(10),log2(20480),0,(500+1),0,(500+1)))

set.seed(1)
res<- runserosim_update(
  simulation_settings=list("t_start"=1,"t_end"=max(times)),
  demography,
  observation_times=tibble(i=rep(1:max(demography$i),each=3),t=rep(c(250,330,365),N)),
  foe_pars,
  biomarker_map,
  model_pars,
  exposure_model=exposure_model_fixed,
  immunity_model=immunity_model_ifxn_biomarker_prot,
  antibody_model=antibody_model_biphasic,
  observation_model=observation_model_continuous_bounded_noise,
  draw_parameters=draw_parameters_random_fx_biomarker_dep,
  VERBOSE = T,
  
  ## Other arguments needed
  max_events=c(4),
  bounds=bounds
)

## turn into the format we are using in the mixture model approach
observed_states <- res$observed_biomarker_states
observed_states <- reshape(observed_states[,-4],idvar=c("i","b"),timevar = c("t"),direction="wide")
observed_states <- reshape(observed_states,idvar=c("i"),timevar = c("b"),direction="wide")
names(observed_states) <- c("i","HAI_t1","HAI_t2","HAI_t3","IgG_t1","IgG_t2","IgG_t3","IgM_t1","IgM_t2","IgM_t3")
infectedHosts <- res$immune_histories_long[(res$immune_histories_long$t>300)&(res$immune_histories_long$value==1),]$i
observed_states <- observed_states %>% mutate(infected = if_else(i %in% infectedHosts,1,0))
observed_states <- observed_states %>% mutate(numInfections = rowSums(data.frame(res$immune_histories),na.rm = T))

## save data 
saveRDS(observed_states,"sim_noise_nonLog.RDS")

##### re-run with more noise #####
# in particular we are going to bump up the variation in all standard deviations related to boosting and waning by 50%
model_pars_original[1,] <- c("ifxn","HAI","boost_long",4,1,"log-normal")
model_pars_original[2,] <- c("ifxn","HAI","boost_short",0.75,.1,"log-normal")
model_pars_original[3,] <- c("ifxn","HAI","wane_long",.0001,.0001,"log-normal")
model_pars_original[4,] <- c("ifxn","HAI","wane_short",.005,.001,"log-normal")
model_pars_original[5,] <- c("ifxn","HAI","biomarker_prot_midpoint",8,NA,"")
model_pars_original[6,] <- c("ifxn","HAI","biomarker_prot_width",.8,NA,"")
model_pars_original[7,] <- c(NA,"HAI","obs_sd",NA,.15,"normal")
model_pars_original[8,] <- c("ifxn","IgG","boost_long",5,30,"log-normal")
model_pars_original[9,] <- c("ifxn","IgG","boost_short",80,30,"log-normal")
model_pars_original[10,] <- c("ifxn","IgG","wane_long",.005,.005,"log-normal")
model_pars_original[11,] <- c("ifxn","IgG","wane_short",.02,.05,"log-normal")
model_pars_original[12,] <- c(NA,"IgG","obs_sd",NA,2.5,"normal")
model_pars_original[13,] <- c("ifxn","IgM","boost_long",5,30,"log-normal")
model_pars_original[14,] <- c("ifxn","IgM","boost_short",50,70,"log-normal")
model_pars_original[15,] <- c("ifxn","IgM","wane_long",.005,.005,"log-normal")
model_pars_original[16,] <- c("ifxn","IgM","wane_short",.0075,.05,"log-normal")
model_pars_original[17,] <- c(NA,"IgM","obs_sd",NA,2.5,"normal")
model_pars_original[18,] <- c("ifxn","HAI","biomarker_ceiling_threshold",10,NA,"")
model_pars_original[19,] <- c("ifxn","HAI","biomarker_ceiling_gradient",.15,NA,"")
model_pars_original[20,] <- c("ifxn","IgG","biomarker_ceiling_threshold",80,NA,"")
model_pars_original[21,] <- c("ifxn","IgG","biomarker_ceiling_gradient",.01,NA,"")
model_pars_original[22,] <- c("ifxn","IgM","biomarker_ceiling_threshold",75,NA,"")
model_pars_original[23,] <- c("ifxn","IgM","biomarker_ceiling_gradient",.01,NA,"")
model_pars<-reformat_biomarker_map(model_pars_original)
model_pars$mean <- as.numeric(model_pars$mean)
model_pars$sd <- as.numeric(model_pars$sd)

set.seed(1)
res<- runserosim_update(
  simulation_settings=list("t_start"=1,"t_end"=max(times)),
  demography,
  observation_times=tibble(i=rep(1:max(demography$i),each=3),t=rep(c(250,330,365),N)),
  foe_pars,
  biomarker_map,
  model_pars,
  exposure_model=exposure_model_fixed,
  immunity_model=immunity_model_ifxn_biomarker_prot,
  antibody_model=antibody_model_biphasic,
  observation_model=observation_model_continuous_bounded_noise,
  draw_parameters=draw_parameters_random_fx_biomarker_dep,
  VERBOSE = T,
  
  ## Other arguments needed
  max_events=c(4),
  bounds=bounds
)

# get into right format and save as RDS
observed_states <- res$observed_biomarker_states
observed_states <- reshape(observed_states[,-4],idvar=c("i","b"),timevar = c("t"),direction="wide")
observed_states <- reshape(observed_states,idvar=c("i"),timevar = c("b"),direction="wide")
names(observed_states) <- c("i","HAI_t1","HAI_t2","HAI_t3","IgG_t1","IgG_t2","IgG_t3","IgM_t1","IgM_t2","IgM_t3")
infectedHosts <- res$immune_histories_long[(res$immune_histories_long$t>300)&(res$immune_histories_long$value==1),]$i
observed_states <- observed_states %>% mutate(infected = if_else(i %in% infectedHosts,1,0))
observed_states <- observed_states %>% mutate(numInfections = rowSums(data.frame(res$immune_histories),na.rm = T))

saveRDS(observed_states,"sim_noise_more_nonLog.RDS")

##### re-run with even more noise #####
# in particular we are going to bump up the variation in all standard deviations related to boosting and waning by 50%
model_pars_original[1,] <- c("ifxn","HAI","boost_long",4,1,"log-normal")
model_pars_original[2,] <- c("ifxn","HAI","boost_short",0.75,.1,"log-normal")
model_pars_original[3,] <- c("ifxn","HAI","wane_long",.0001,.0001,"log-normal")
model_pars_original[4,] <- c("ifxn","HAI","wane_short",.005,.001,"log-normal")
model_pars_original[5,] <- c("ifxn","HAI","biomarker_prot_midpoint",8,NA,"")
model_pars_original[6,] <- c("ifxn","HAI","biomarker_prot_width",.8,NA,"")
model_pars_original[7,] <- c(NA,"HAI","obs_sd",NA,.3,"normal")
model_pars_original[8,] <- c("ifxn","IgG","boost_long",5,30,"log-normal")
model_pars_original[9,] <- c("ifxn","IgG","boost_short",80,30,"log-normal")
model_pars_original[10,] <- c("ifxn","IgG","wane_long",.005,.005,"log-normal")
model_pars_original[11,] <- c("ifxn","IgG","wane_short",.02,.05,"log-normal")
model_pars_original[12,] <- c(NA,"IgG","obs_sd",NA,5,"normal")
model_pars_original[13,] <- c("ifxn","IgM","boost_long",5,30,"log-normal")
model_pars_original[14,] <- c("ifxn","IgM","boost_short",50,70,"log-normal")
model_pars_original[15,] <- c("ifxn","IgM","wane_long",.005,.005,"log-normal")
model_pars_original[16,] <- c("ifxn","IgM","wane_short",.005,.05,"log-normal")
model_pars_original[17,] <- c(NA,"IgM","obs_sd",NA,5,"normal")
model_pars_original[18,] <- c("ifxn","HAI","biomarker_ceiling_threshold",10,NA,"")
model_pars_original[19,] <- c("ifxn","HAI","biomarker_ceiling_gradient",.15,NA,"")
model_pars_original[20,] <- c("ifxn","IgG","biomarker_ceiling_threshold",80,NA,"")
model_pars_original[21,] <- c("ifxn","IgG","biomarker_ceiling_gradient",.01,NA,"")
model_pars_original[22,] <- c("ifxn","IgM","biomarker_ceiling_threshold",75,NA,"")
model_pars_original[23,] <- c("ifxn","IgM","biomarker_ceiling_gradient",.01,NA,"")
model_pars<-reformat_biomarker_map(model_pars_original)
model_pars$mean <- as.numeric(model_pars$mean)
model_pars$sd <- as.numeric(model_pars$sd)

set.seed(1)
res<- runserosim_update(
  simulation_settings=list("t_start"=1,"t_end"=max(times)),
  demography,
  observation_times=tibble(i=rep(1:max(demography$i),each=3),t=rep(c(250,330,365),N)),
  foe_pars,
  biomarker_map,
  model_pars,
  exposure_model=exposure_model_fixed,
  immunity_model=immunity_model_ifxn_biomarker_prot,
  antibody_model=antibody_model_biphasic,
  observation_model=observation_model_continuous_bounded_noise,
  draw_parameters=draw_parameters_random_fx_biomarker_dep,
  VERBOSE = T,
  
  ## Other arguments needed
  max_events=c(4),
  bounds=bounds
)

# get into right format and save as RDS
observed_states <- res$observed_biomarker_states
observed_states <- reshape(observed_states[,-4],idvar=c("i","b"),timevar = c("t"),direction="wide")
observed_states <- reshape(observed_states,idvar=c("i"),timevar = c("b"),direction="wide")
names(observed_states) <- c("i","HAI_t1","HAI_t2","HAI_t3","IgG_t1","IgG_t2","IgG_t3","IgM_t1","IgM_t2","IgM_t3")
infectedHosts <- res$immune_histories_long[(res$immune_histories_long$t>300)&(res$immune_histories_long$value==1),]$i
observed_states <- observed_states %>% mutate(infected = if_else(i %in% infectedHosts,1,0))
observed_states <- observed_states %>% mutate(numInfections = rowSums(data.frame(res$immune_histories),na.rm = T))

saveRDS(observed_states,"sim_noise_even_more_nonLog.RDS")
