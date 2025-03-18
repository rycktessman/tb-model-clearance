library(deSolve)
library(dampack)
library(mvtnorm)
library(flexsurv)
library(lhs)
library(matrixStats)
library(dplyr)
library(stringr)

#specify options
part <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) #array
print(part)
scenario <- Sys.getenv('scenario') #base or clearance
print(scenario)
B <- as.numeric(Sys.getenv('n_samples')) #number of samples per round (initial = 10*B)
print(B)
version <- Sys.getenv('version') #whether to induce incidence declines/the method. no decline="", transmission rate="beta_decline", tx increases="omega_increase"
print(version)

source('code/model_setup.R')
source('code/model_functions.R')

pars.baseline <- pars
par.range <- read.csv(paste0('data/par_range_', scenario, '_', version, '.csv'), row.names=1)[1:3,]

n_extra_states <- 16 #number of extra states/outcomes tracked in the model (like incidence, etc.)

tic <- Sys.time()

outputs <- samples_like(B, par.range, pars.baseline)

toc <- Sys.time()
print(toc-tic)

path <- paste0("samples_", version, "/")
fsave <- paste(path, 'samples_',part, '_', scenario, '.rda',sep='')
save(outputs, file=fsave)
