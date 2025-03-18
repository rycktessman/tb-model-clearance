library(dplyr)
library(stringr)

scenario <- Sys.getenv('scenario') #base or clearance
path_out <- Sys.getenv('path') #file path where IMIS samples are stored

pars_all <- data.frame()
out_all <- data.frame()
ystart_all <- data.frame()
yend_all <- data.frame()

setwd(path_out)
completed_files <- list.files(pattern=paste0(scenario, ".rda"))
completed_files <- completed_files[!(str_detect(completed_files, "temp"))]
print(completed_files)
completed <- lapply(1:length(completed_files), function(x)
  str_split(completed_files[[x]], "_")[[1]][[2]])
completed <- as.numeric(unlist(completed))

for(i in completed) {
  print(i)
  
  fname <- paste0("imis_", i, "_", scenario, ".rda")
  load(fname)
  
  pars <- post$resample
  out <- post$out
  
  #only take ystart and yend that correspond to the posteriors
  ystart_everything <- post$ystart_all
  ystart <- left_join(out %>% select(array, sample, round), ystart_everything)
  
  yend_everything <- post$yend_all
  if("samples" %in% names(yend_everything)) {
    yend_everything <- yend_everything %>% rename(sample=samples)
  }
  out_temp <- rbind(out %>% mutate(t=20), out %>% mutate(t=21))
  yend <- left_join(out_temp %>% select(array, sample, round, t),
                    yend_everything, by=c("array", "sample", "round", "t"))
  
  pars_all <- bind_rows(pars_all, pars %>% mutate(array2=i))
  out_all <- bind_rows(out_all, out %>% mutate(array2=i))
  ystart_all <- bind_rows(ystart_all, ystart %>% mutate(array2=i))
  yend_all <- bind_rows(yend_all, yend %>% mutate(array2=i))
}

write.csv(pars_all, paste0("params_post_", scenario, ".csv"), row.names=F)
write.csv(out_all, paste0("out_post_", scenario, ".csv"), row.names=F)
write.csv(ystart_all, paste0("ystart_post_", scenario, ".csv"), row.names=F)
write.csv(yend_all, paste0("yend_post_", scenario, ".csv"), row.names=F)

