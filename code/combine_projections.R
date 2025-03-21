library(deSolve)
library(mvtnorm)
library(flexsurv)
library(lhs)
library(matrixStats)
library(dplyr)
library(stringr)

#specify options
scenario <- Sys.getenv('scenario') #base or clearance
print(scenario)
path <- Sys.getenv('path') #directory where files stored - determines which set of targets and whether beta_decline, omega_increase, or neither
print(path)
repeat_acf <- 2

files <- list.files(path=path, pattern=paste0("projections_time_", scenario, "_", repeat_acf, "yrs_"))
print(files)
file_numbers <- as.numeric(str_extract(files, "(?<=_)[0-9]+(?=\\.Rda)"))
print(file_numbers)
arrays <- max(file_numbers)
print(arrays)

out_time_comb <- data.frame()
outputs_comb <- data.frame()
for(i in 1:arrays) {
  print(i)
  load(paste0(path, "/projections_time_", scenario, "_", repeat_acf, "yrs_", i, ".Rda"))
  outputs <- read.csv(paste0(path, "/projections_", scenario, "_", repeat_acf, "yrs_", i, ".csv"))
  
  out_time_comb <- bind_rows(out_time_comb, out_time_all)
  outputs_comb <- bind_rows(outputs_comb, outputs)
}

#keep a selection of outcomes at all month time steps
out_time_sub <- out_time_comb %>% select(time, pop, inc, prev, mort, prev_ltbi, ari_all,
                                        scenario, array, sample, round, array2, count)
#calculate cumulative and incremental outcomes
out_time_sub <- out_time_sub %>% 
  group_by(time, array, sample, round, array2) %>%
  mutate(inc_decline=(inc[scenario=="No Intervention"] - inc)/inc[scenario=="No Intervention"],
         mort_decline=(mort[scenario=="No Intervention"] - mort)/mort[scenario=="No Intervention"],
         inc_decline_prop_acf=if_else(inc_decline[scenario=="ACF"]<0, 0, inc_decline[scenario=="ACF"]/inc_decline[scenario=="ACF & TPT"]),
         mort_decline_prop_acf=if_else(mort_decline[scenario=="ACF"]<0, 0, mort_decline[scenario=="ACF"]/mort_decline[scenario=="ACF & TPT"]))
out_time_sub <- out_time_sub %>% ungroup() %>%
  group_by(scenario, array, sample, round, array2) %>%
  mutate(inc_decline0=(inc[time==0] - inc)/inc[time==0],
         mort_decline0=(mort[time==0] - mort)/mort[time==0])

#duplicate as needed to reintroduce weights
out_time_sub <- out_time_sub[rep(seq(nrow(out_time_sub)), out_time_sub$count), ]

#take means and CIs
out_time_sum <- out_time_sub %>% group_by(time, scenario) %>% 
  select(-c(array, sample, round, array2, count)) %>%
  summarise_all(.funs=list(mean=~mean(., na.rm=T),
                           lb=~quantile(., p=0.025, na.rm=T),
                           ub=~quantile(., p=0.975, na.rm=T)))

#save full distributions for selected outputs at 10 yrs
out_all <- out_time_sub %>% filter(time==120)


write.csv(out_time_sum, file=paste0(path, "/projections_time_sum_", scenario, "_", repeat_acf, "yrs.csv"), row.names=F)
write.csv(out_all, file=paste0(path, "/projections_all_10_", scenario, "_", repeat_acf, "yrs.csv"), row.names=F)
