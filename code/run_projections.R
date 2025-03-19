library(deSolve)
library(mvtnorm)
library(flexsurv)
library(lhs)
library(matrixStats)
library(dplyr)

n_extra_states <- 16 #number of extra states/outcomes tracked in the model (like incidence, etc.)
time_proj <- 11 #run for 11 years
samples_per_array <- 100

#specify additional options in the bash script
part <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) #array; 1 through 20
print(part)
scenario <- Sys.getenv('scenario') #base or clearance
print(scenario)
path <- Sys.getenv('path') #directory where files stored - determines which set of targets and whether beta_decline, omega_increase, or neither
print(path)

repeat_acf <- 2 #years between ACF rounds - set to 999 for 1-time ACF
acf_years <- seq(from=0, to=time_proj, by=repeat_acf)
start <- (part-1)*samples_per_array+1
end <- part*samples_per_array
print(start)
print(end)
print(acf_years)

source('code/model_setup.R')
source('code/model_functions.R')

acf_ages <- dimnames(matind)[[3]] #all ages receive case finding/TPT

## read posterior parameters for simulation
post_all <- read.csv(paste0(path, "/params_post_", scenario, ".csv")) 
ystart_all <- read.csv(paste0(path, "/ystart_post_", scenario, ".csv")) 

post_unique <- post_all %>% group_by_all() %>% summarise(count=n()) #take only unique rows, but count them to duplicate later
post_unique <- post_unique %>% ungroup() %>% arrange(array, array2, round, sample)
ystart_unique <- ystart_all %>% group_by_all() %>% summarise(count=n()) #take only unique rows, but count them to duplicate later
ystart_unique <- ystart_unique %>% ungroup() %>% arrange(array, array2, round, sample)
ystart_unique <- ystart_unique %>% select(-c(array, sample, round, array2, count))
ystart_unique[,((max(matind)+1):ncol(ystart_unique))] <- 0

if(end>nrow(post_unique)) {
  end <- nrow(post_unique)
}
if(start>=end) {
  stop("the script ends")
} else {
  post_sub <- post_unique[start:end,]
  ystart_sub <- ystart_unique[start:end,]
  
  pars.baseline <- pars
  
  outputs_all <- data.frame()
  out_time_all <- data.frame()
  for(js in 1:nrow(post_sub)){
    print(js)
    pars_use <- post_sub[js,]
    
    #BURN IN 20 YEARS OF TX IMPROVEMENTS (need to do this again to get outputs, like incidence, over time)
    ystart.burn <- unlist(ystart_sub[js, ])
    out.burn <- model_output_proj(pars_use, pars.baseline, ystart.burn)
    
    output_names <- names(out.burn$outputs)
    
    # NO INTERVENTION
    ystart.noint <- unlist(out.burn$y[nrow(out.burn$y), 2:ncol(out.burn$y)])
    y.lastyears <- out.burn$y[(nrow(out.burn$y)-120):nrow(out.burn$y), ]
    y.lastyears[,"time"] <- seq(from=-10, to=0, by=1/12)
    #don't continue to increase omega after 20 year end of burn-in
    if("increase.omega" %in% names(pars_use)) {
      pars_use$omega <- pars_use$omega*(1+20*pars_use$increase.omega)
      pars_use$increase.omega <- 0
    }
    #do keep declining beta, but incorporate 20 years burn-in to base beta value for projections
    if("decline.beta" %in% names(pars_use)) {
      pars_use$beta <- pars_use$beta*exp(-20*pars_use$decline.beta)
    }
    out.noint <- model_output_proj(pars_use, pars.baseline, ystart.noint, tmax=time_proj)
    #also calculate outputs over time
    out.noint$out_time <- calc_outputs_time(rbind(head(y.lastyears, -1), out.noint$y), 12)
    out.noint$out_time$time <- out.noint$out_time$time - 121
    
    # ACTIVE CASE FINDING
    names.nat.hist <- rownames(matind)
    tb <- c(which(names.nat.hist %in% paste0('AA', 1:n.active)), 
            which(names.nat.hist %in% paste0('SA', 1:n.active)))
    #define additional parameters - not used if % cured by ACF is specified in the SLURM script
    acf_covg <- 0.7
    test_sens <- 0.65
    p_cure <- pars_use$k
    cures_acf <- acf_covg*test_sens*p_cure
    ystart.acf <- ystart.noint
    for(i in 1:length(acf_years)) {
      cur_acf_time <- acf_years[[i]]
      next_acf_time <- ifelse(i==length(acf_years), time_proj, acf_years[[i+1]])
      if(cur_acf_time==next_acf_time) {
        break
      }
      #ACF happens up front: AA and SA -> R 
      for(a in acf_ages) {
        ystart.acf[matind["R",1,a]] <- sum(cures_acf*ystart.acf[matind[tb,1,a]]) + ystart.acf[matind["R",1,a]] #sum over symptomatic and asymptomatic TB for susceptible pop
        ystart.acf[matind["R",2,a]] <- sum(cures_acf*ystart.acf[matind[tb,2,a]]) + ystart.acf[matind["R",2,a]] #sum over symptomatic and asymptomatic TB for resistor pop
        ystart.acf[matind[tb, , a]] <- (1-cures_acf)* ystart.acf[matind[tb, , a]]
      }
      #then run model for 20 years with updated starting state
      out.acf_tmp <- model_output_proj(pars_use, pars.baseline, ystart.acf, tmax=(next_acf_time - cur_acf_time))
      
      #bind to out.acf and update ystart
      if(i==1) {
        out.acf <- out.acf_tmp
      } else {
        out.acf$y <- rbind(out.acf$y %>% filter(time!=cur_acf_time), out.acf_tmp$y %>% mutate(time=time+cur_acf_time))
      }
      out.acf$outputs <- unlist(out.acf_tmp$outputs)
      names(out.acf$outputs) <- output_names
      ystart.acf <- unlist(tail(out.acf_tmp$y, 1)[2:ncol(out.acf_tmp$y)])
    }
    #also calculate outputs over time
    out.acf$out_time <- calc_outputs_time(rbind(head(y.lastyears, -1), out.acf$y), 12)
    out.acf$out_time$time <- out.acf$out_time$time - 121
    
    #TPT PROVISION
    ltbi <- c(which(names.nat.hist %in% c("L1", "L2", "L3"))) #needs editing if >3 LTBI states
    #define additional parameters
    tpt_covg <- acf_covg
    tst_read <- 0.865
    tst_sens <- 0.9 
    tpt_uptake <- 0.85
    tpt_completion <- 0.85
    tpt_efficacy <- 0.69 
    cures_tpt <- tpt_covg*tst_read*tst_sens*tpt_uptake*tpt_completion*tpt_efficacy
    #add TPT to ACF
    ystart.tpt <- ystart.noint
    for(i in 1:length(acf_years)) {
      cur_acf_time <- acf_years[[i]]
      next_acf_time <- ifelse(i==length(acf_years), time_proj, acf_years[[i+1]])
      if(cur_acf_time==next_acf_time) {
        break
      }
      #ACF and TPT happen up front: AA, SA, L1, L2 -> R 
      for(a in acf_ages) {
        ystart.tpt[matind["R",1,a]] <- sum(cures_acf*ystart.tpt[matind[tb,1,a]]) + ystart.tpt[matind["R",1,a]] #sum over symptomatic and asymptomatic TB for susceptible pop
        ystart.tpt[matind["R",2,a]] <- sum(cures_acf*ystart.tpt[matind[tb,2,a]]) + ystart.tpt[matind["R",2,a]] #sum over symptomatic and asymptomatic TB for resistor pop
        ystart.tpt[matind[tb, , a]] <- (1-cures_acf)* ystart.tpt[matind[tb, , a]] 
        ystart.tpt[matind["R",1,a]] <- sum(cures_tpt*ystart.tpt[matind[ltbi,1,a]]) + ystart.tpt[matind["R",1,a]] #sum over LTBI states for susceptible pop
        ystart.tpt[matind["R",2,a]] <- sum(cures_tpt*ystart.tpt[matind[ltbi,2,a]]) + ystart.tpt[matind["R",2,a]] #sum over LTBI states for resistor pop
        ystart.tpt[matind[ltbi, , a]] <- (1-cures_tpt)*ystart.tpt[matind[ltbi, , a]] 
      }
      #then run model for 20 years with updated starting state
      out.tpt_tmp <- model_output_proj(pars_use, pars.baseline, ystart.tpt, tmax=(next_acf_time - cur_acf_time))
      #bind to out.tpt and update ystart
      if(i==1) {
        out.tpt <- out.tpt_tmp
      } else {
        out.tpt$y <- rbind(out.tpt$y %>% filter(time!=cur_acf_time), out.tpt_tmp$y %>% mutate(time=time+cur_acf_time))
      }
      out.tpt$outputs <- unlist(out.tpt_tmp$outputs)
      names(out.tpt$outputs) <- output_names
      ystart.tpt <- unlist(tail(out.tpt_tmp$y, 1)[2:ncol(out.tpt_tmp$y)])
    }
    #also calculate outputs over time
    out.tpt$out_time <- calc_outputs_time(rbind(head(y.lastyears, -1), out.tpt$y), 12)
    out.tpt$out_time$time <- out.tpt$out_time$time - 121
    
    #combine everything
    out.noint$outputs <- unlist(out.noint$outputs)
    names(out.noint$outputs) <- output_names
    outputs <- rbind(as.data.frame(t(unlist(out.noint$outputs))) %>% mutate(scenario="No Intervention"),
                     as.data.frame(t(out.acf$outputs)) %>% mutate(scenario="ACF"),
                     as.data.frame(t(out.tpt$outputs)) %>% mutate(scenario="ACF & TPT")) %>%
      mutate(array=pars_use$array, sample=pars_use$sample, round=pars_use$round, array2=pars_use$array2, count=pars_use$count)
    
    out_time <- rbind(out.noint$out_time %>% mutate(scenario="No Intervention"),
                      out.acf$out_time %>% mutate(scenario="ACF"),
                      out.tpt$out_time %>% mutate(scenario="ACF & TPT")) %>%
      mutate(array=pars_use$array, sample=pars_use$sample, round=pars_use$round, array2=pars_use$array2, count=pars_use$count)
    
    outputs_all <- bind_rows(outputs_all, outputs)
    out_time_all <- bind_rows(out_time_all, out_time)
  }
  
  #fix names
  names(outputs_all) <- gsub("\\..*$", "", names(outputs_all))
  
  file_ending <- paste0(scenario, "_", repeat_acf, "yrs_", part)
  write.csv(outputs_all, file=paste0(path, "/projections_", file_ending, ".csv"), row.names=F)
  save(out_time_all, file=paste0(path, "/projections_time_", file_ending, ".Rda"))
  
  #to check after 1 iteration
  if(FALSE) {
    library(ggplot2)
    ggplot(out_time, aes(x=time, y=inc, color=scenario)) +
      geom_line() + 
      theme_bw()
    
    ggplot(out_time, aes(x=time, y=prev, color=scenario)) +
      geom_line() + 
      theme_bw()
    
    ggplot(out_time, aes(x=time, y=mort, color=scenario)) +
      geom_line() + 
      theme_bw()
    
    ggplot(out_time, aes(x=time, y=notif, color=scenario)) +
      geom_line() + 
      theme_bw()
  }
  
}
