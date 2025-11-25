library(dplyr)
library(stringr)
library(matrixStats)
library(mvtnorm)
library(tmvtnorm) #version that includes truncated multivariate normal
library(lhs)
library(data.table)
library(deSolve)
library(dampack)

age_groups <- 2
n_extra_states <- 20

scenario <- Sys.getenv('scenario') #"base" or "progression"
version <- Sys.getenv('version') #"beta_decline" or "omega_increase"
targets_version <- Sys.getenv('targets_set') #5, 8, or 9
print(targets_version)
B <- as.numeric(Sys.getenv('B')) #200
print(B)
number_k <- as.numeric(Sys.getenv('number_k'))
print(number_k)
p_clear_react <- as.numeric(Sys.getenv('p_clear_react')) #% of cleared infections remaining IGRA+, only used in progression scenario
print(p_clear_react)
file_label <- Sys.getenv('file_lab')
print(file_label)
B.re <- 1000
D <- 0

source('code/model_setup.R')
source('code/model_functions.R')

pars.baseline <- pars
pars.baseline[["p_clear_react"]] <- p_clear_react
if(scenario=="base") {
  par_range <- read.csv(paste0('data/par_range_', version, '.csv'), row.names=1)[1:3,]
} else if(scenario=="progression") {
  par_range <- read.csv(paste0('data/par_range_clearance_', version, '.csv'), row.names=1)[1:3,]
}
#make sure the sampling ranges are all calculated correctly
for(i in colnames(par_range)) {
  par_range["diff", i] <- as.numeric(par_range["hi", i]) - as.numeric(par_range["low", i])
}

#remove params that don't vary in calibration from par_range
par_range <- par_range[,par_range["diff",]!=0]

#"prior" and "sample.prior" are already defined as functions
#need to define likelihood function - wrapper on model_output_like
gen_out_like <- function(pars) {
  out_round <- data.frame()
  ystart_round <- data.frame()
  yend_round <- data.frame()
  for(i in 1:nrow(pars)) {
    print(i)
    out_i_all <- model_output_like(pars[i,], pars.baseline)
    out_i <- unlist(out_i_all$outputs)
    ystart_i <- out_i_all$ystart
    yend_i <- tail(out_i_all$y, 2)

    out_round <- rbind(out_round, out_i)
    ystart_round <- rbind(ystart_round, ystart_i)
    yend_round <- rbind(yend_round, yend_i) 
  }
  #fix names
  names(out_round) <- names(out_i)
  names(ystart_round) <- paste0("X", 1:ncol(ystart_round))
  names(yend_round) <- c("time", names(ystart_round))
  return(list("out"=out_round, "ystart"=ystart_round, "yend"=yend_round))
}

#run IMIS and save
post <- IMIS(B=B, B.re=B.re, number_k=number_k, D=0, targets_set=targets_version)
save(post, file=paste0("output/imis", "_", file_label, "_", scenario, ".Rda"))
