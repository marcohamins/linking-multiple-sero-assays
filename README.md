# Code for *Linking multiple serological assays to infer dengue virus infections from paired samples using mixture models*

An outline of the organization of this repo:
1. /simulate contains the sim_antibody_response.R file which runs forward simulations to create simulated antibody data
2. /models/run_stan contains a few examples of how the mixture models are run on simulated data and Thai cohort data, with results getting saved to stan_output
3. /models/stan_files contain the .stan files used to run the mixture models
4. /stan_output contains an example output from Model 8 where all acute and convalescent IgG, IgM, and HAI data are incorporated and is where all output from the run_stan folder are saved
5. /data folder contains simulated data from code in /simulate folder, and antibody data from Thai cohorts